# -*- coding: utf-8 -*-
"""
Reconstruction of a lattice from a 3D image (ringspy.reconstruct), checked against synthetic
volumes of known geometry (ringspy.voxelize)
"""

import os
import subprocess
import sys

import numpy as np
import pytest
from scipy.spatial import cKDTree

import ringspy
from ringspy import reconstruction as R
from ringspy.synthetic import voxelize
from helpers import OUTPUT_SUFFIXES, assert_same_files

VOXEL = 0.002
# a small wood cube: 0.3 mm, cells of 20 to 40 um with walls of 6 and 10 um, 150^3 voxels of 2 um
SMALL = ringspy.Config(geoName="small", r_max=1.0, nrings=2, box_size=0.3, box_center=(0.55, 0.0), z_max=0.3,
                       seed=3, iter_max=300, stlFlag=False, inpFlag=False)
TILT = (0.06, -0.03)


@pytest.fixture(scope="module")
def synthetic():
    """A synthetic scan: tilted, twisted and wobbling fibers"""
    return voxelize(SMALL, VOXEL, tilt=TILT, twist=0.05, wobble_amplitude=0.004, wobble_wavelength=0.2, seed=3)


@pytest.fixture(scope="module")
def reconstructed(synthetic, tmp_path_factory):
    """The reconstruction of the synthetic scan (run once), with the output folder made absolute"""
    workdir = tmp_path_factory.mktemp("reconstructed")
    cwd = os.getcwd()
    os.chdir(workdir)
    try:
        result = ringspy.reconstruct(synthetic, geoName="rec", nsegments=2)
    finally:
        os.chdir(cwd)
    result.output_dir = workdir / result.output_dir
    return result


def true_vertices(synthetic, z, margin=0.0):
    z_min, z_max = synthetic.box[4], synthetic.box[5]
    return synthetic.vertices_at(min(max(z, z_min + 1e-9), z_max - 1e-9), margin)


# ---------------------------------------------------------------- synthetic volumes

def test_voxelize_volume(synthetic):
    assert synthetic.volume.dtype == bool and synthetic.shape == (150, 150, 150)
    assert synthetic.voxel_size == (VOXEL, VOXEL, VOXEL)
    np.testing.assert_allclose(synthetic.box, [0.4, 0.7, -0.15, 0.15, 0.0, 0.3])
    assert 0.5 < synthetic.porosity < 0.9
    assert len(synthetic.sites) == len(synthetic.wall_thickness) > 50
    assert set(np.round(synthetic.wall_thickness, 6)) <= {0.006, 0.01}
    vertices = synthetic.vertices_at(0.15)
    assert len(vertices) > 20
    assert (vertices[:, 0] > 0.4).all() and (vertices[:, 0] < 0.7).all()
    # the fibers lean: the cells have moved from the bottom to the top of the volume
    moved = synthetic.sites_at(0.3) - synthetic.sites_at(0.0)
    assert np.linalg.norm(moved, axis=1).mean() > 0.3*np.linalg.norm(TILT)*0.3


def test_voxelize_from_arrays_and_grayscale():
    sites = np.array([[0.0, 0.0], [0.1, 0.0], [0.0, 0.1], [0.1, 0.1], [0.05, 0.05]])
    geometry = {"sites": sites, "wall_thickness": 0.01, "box": (0.0, 0.1, 0.0, 0.1, 0.0, 0.02)}
    plain = voxelize(geometry, 0.002)
    assert plain.shape == (50, 50, 10) and plain.volume.dtype == bool
    assert 0.1 < plain.volume.mean() < 0.6
    gray = voxelize(geometry, 0.002, noise=0.1, blur=0.5, seed=1)
    assert gray.volume.dtype == np.float32 and gray.shape == plain.shape
    assert abs(gray.porosity - plain.porosity) < 0.1
    with pytest.raises(ValueError):
        voxelize(geometry, 0.0)


# ---------------------------------------------------------------- image processing

def test_otsu_threshold_separates_two_modes():
    rng = np.random.default_rng(0)
    values = np.concatenate((rng.normal(0.2, 0.05, 5000), rng.normal(0.8, 0.05, 2000)))
    threshold = R.otsu_threshold(values)
    assert 0.35 < threshold < 0.65


def test_binarize_bright_and_dark():
    image = np.zeros((6, 6, 6), dtype=np.float32)
    image[2:4] = 1.0
    bright, threshold = R.binarize(image, "bright")
    dark, _ = R.binarize(image, "dark", threshold=0.5)
    assert bright.sum() == 72 and dark.sum() == 6*6*6 - 72 and 0 < threshold < 1
    walls, threshold = R.binarize(image > 0.5, "bright")
    assert threshold is None and walls.sum() == 72
    with pytest.raises(ValueError):
        R.binarize(image, "shiny")


def test_detect_longitudinal_axis(synthetic):
    for axis in range(3):
        moved = np.moveaxis(synthetic.volume, 2, axis)
        assert R.detect_longitudinal_axis(moved) == axis


SITES = np.array([[8.0, 8.0], [32.0, 8.0], [8.0, 32.0], [32.0, 32.0], [20.0, 20.0]])


def voronoi_labels(sites=SITES, n=40, thickness=3.0):
    """A 2D cell image: the Voronoi cells of `sites` (pixels), walls `thickness` pixels thick between them.
    The default is a diamond-shaped cell with corners at (20, 8), (32, 20), (20, 32), (8, 20), between four cells cut
    by the sides of the image."""
    i, j = np.meshgrid(np.arange(n) + 0.5, np.arange(n) + 0.5, indexing="ij")
    pixels = np.column_stack((i.reshape(-1), j.reshape(-1)))
    distance, index = cKDTree(sites).query(pixels, k=len(sites))
    walls = np.zeros(len(pixels), dtype=bool)
    for k in range(1, len(sites)):
        separation = np.linalg.norm(sites[index[:, k]] - sites[index[:, 0]], axis=1)
        walls |= (distance[:, k]**2 - distance[:, 0]**2)/(2*separation) < thickness/2
    labels = (index[:, 0] + 1).reshape(n, n)
    walls = walls.reshape(n, n)
    labels[walls] = 0
    return labels, walls


def test_slice_network_of_a_cell_between_four_cells():
    labels, walls = voronoi_labels()
    net = R.slice_network(labels, walls, (0.1, 0.1), (0.0, 0.0))
    inside = net.kind == 0
    assert inside.sum() == 4 and (net.kind == 1).sum() == 4 and (net.kind == 2).sum() == 4
    expected = np.array([[2.0, 0.8], [3.2, 2.0], [2.0, 3.2], [0.8, 2.0]])
    distance = cKDTree(expected).query(net.points[inside])[0]
    assert distance.max() < 0.1  # within a pixel of the true corners
    assert all(5 in labels and len(labels) == 3 for labels in [net.label_sets[i] for i in np.nonzero(inside)[0]])
    # the 4 walls around the center cell, 4 walls from its corners to the sides, 8 pieces of the perimeter
    assert len(net.ridges) == 16
    finite = [r for r, (a, b) in enumerate(net.ridges) if inside[a] and inside[b]]
    cut = [r for r, (a, b) in enumerate(net.ridges) if inside[a] != inside[b]]
    assert len(finite) == 4 and len(cut) == 4
    assert sorted(net.pairs[r] for r in finite) == [(1, 5), (2, 5), (3, 5), (4, 5)]
    np.testing.assert_allclose(net.thickness[finite + cut], 0.3, atol=0.1)
    np.testing.assert_allclose(net.median_cell_diameter, 2*np.sqrt(np.median(list(net.cell_area.values()))/np.pi))
    ordered = R.order_network(net, boundary_on=True)
    assert ordered.nvertices_in == 4 and ordered.nboundary_pts == 8 and ordered.nboundary_pts_featured == 4
    assert len(ordered.finite_ridges) == 4 and len(ordered.boundary_ridges) == 12
    assert (ordered.kind[:4] == 0).all() and (ordered.kind[4:8] == 1).all() and (ordered.kind[8:] == 2).all()
    np.testing.assert_array_equal(ordered.ridges[:4], ordered.finite_ridges)
    ordered = R.order_network(net, boundary_on=False)
    assert ordered.nvertex == 8 and ordered.nridge == 8 and ordered.nboundary_pts_featured == 0
    with pytest.raises(ValueError, match="no complete cell"):
        R.order_network(R.slice_network(*voronoi_labels(SITES[:2]), (0.1, 0.1), (0.0, 0.0)))


def test_merge_short_ridges():
    labels, walls = voronoi_labels()
    net = R.slice_network(labels, walls, (0.1, 0.1), (0.0, 0.0))
    # add a second inside point close to the center, joined by a short ridge
    i = int(np.nonzero(net.kind == 0)[0][0])
    points = np.vstack((net.points, net.points[i] + [0.05, 0.0]))
    label_sets = net.label_sets + [frozenset({1, 5, 9})]
    kind = np.append(net.kind, 0)
    ridges = np.vstack((net.ridges, [[i, len(points) - 1]]))
    pairs = net.pairs + [(1, 5)]
    thickness = np.append(net.thickness, 0.3)
    grown = R.SliceNetwork(points, label_sets, kind, ridges, pairs, thickness, net.cell_area, net.box)
    merged = R.merge_short_ridges(grown, tol=0.1)
    assert len(merged.points) == len(net.points) and len(merged.ridges) == len(net.ridges)
    assert merged.label_sets[i] == net.label_sets[i] | {9}
    np.testing.assert_allclose(merged.points[i], net.points[i] + [0.025, 0.0])
    untouched = R.merge_short_ridges(grown, tol=0.01)
    assert len(untouched.points) == len(grown.points)


def test_direction_field_recovers_tilt():
    tilt = (0.08, -0.04)
    straight = voxelize(SMALL, VOXEL, tilt=tilt, seed=3)
    field = R.DirectionField(straight.volume, straight.voxel_size, straight.origin, scale=0.05)
    expected = np.degrees(np.arctan(np.linalg.norm(tilt)))
    assert abs(field.grain_angle_deg - expected) < 0.5
    assert abs((field.azimuth_deg - np.degrees(np.arctan2(tilt[1], tilt[0])) + 180) % 360 - 180) < 5
    assert field.dispersion_deg < 1.0
    # the field moves points as the fibers lean
    moved = field.advect(np.array([[0.55, 0.0]]), 0.05, 0.25)
    np.testing.assert_allclose(moved[0], [0.55 + 0.2*tilt[0], 0.2*tilt[1]], atol=0.002)


def test_fit_midpoints_interpolates_quadratic_beams():
    # a beam through (0,0), (1,1), (2,0) at the 5 layers of one segment: layers 1 and 3 are moved so that
    # the two quadratic Bezier elements pass through the followed positions
    layers = np.array([[[0.0, 0.0]], [[0.5, 0.75]], [[1.0, 1.0]], [[1.5, 0.75]], [[2.0, 0.0]]])
    fitted = R.fit_midpoints(layers, nsegments=1)
    np.testing.assert_allclose(fitted[[0, 2, 4]], layers[[0, 2, 4]])
    np.testing.assert_allclose(fitted[1, 0], 2*layers[1, 0] - (layers[0, 0] + layers[2, 0])/2)
    bezier_mid = (fitted[0, 0] + 2*fitted[1, 0] + fitted[2, 0])/4
    np.testing.assert_allclose(bezier_mid, layers[1, 0])


# ---------------------------------------------------------------- end to end

def test_reconstruct_writes_all_files(reconstructed):
    out = reconstructed.output_dir
    expected = ["rec" + suffix for suffix in OUTPUT_SUFFIXES] + ["rec.stl", "rec.inp", "rec-grains.npz"]
    assert sorted(p.name for p in out.iterdir()) == sorted(expected)
    assert sorted(p.name for p in reconstructed.files) == sorted(expected)
    for name in expected:
        assert (out / name).stat().st_size > 0, name
    log = (out / "rec-report.log").read_text()
    assert "grain angle" in log and "matched by its cells" in log
    assert reconstructed.shape == (150, 150, 150) and reconstructed.longitudinal_axis == 2
    assert reconstructed.threshold is None
    assert reconstructed.nParticles > 20 and reconstructed.nbeamElem > 100
    assert 0 < reconstructed.porosity < 1 and abs(reconstructed.image_porosity - 0.734) < 0.02


def test_reconstruct_follows_the_cell_corners(reconstructed, synthetic):
    grains = np.load(reconstructed.output_dir / "rec-grains.npz")
    layer_vertices, status, kind, z = grains["layer_vertices"], grains["status"], grains["kind"], grains["z_coord"]
    assert layer_vertices.shape == (9, len(kind), 2) and status.shape == (9, len(kind))  # 2*nsegments elements: 9 half-layers
    np.testing.assert_array_equal(layer_vertices, reconstructed.layer_vertices)
    inside = kind == 0
    counts = reconstructed.tracking_counts
    total = sum(counts.values())
    assert counts["matched by its cells"] > 0.8*total
    for layer in range(len(z)):
        truth = true_vertices(synthetic, z[layer])
        distance = cKDTree(truth).query(layer_vertices[layer][inside])[0]
        exact = status[layer][inside] == 0
        assert np.median(distance[exact]) < 1.0*VOXEL            # found corners: within a voxel
        assert np.percentile(distance[exact], 90) < 2.0*VOXEL
        assert distance.max() < 0.5*reconstructed.median_cell_diameter   # all corners: far less than a cell
    # the corners of the box stay put, boundary points stay on the sides
    corners = kind == 2
    assert np.all(layer_vertices[:, corners] == layer_vertices[0, corners])
    x0, x1, y0, y1 = 0.4, 0.7, -0.15, 0.15
    on_side = np.isclose(layer_vertices[..., 0], x0) | np.isclose(layer_vertices[..., 0], x1) | \
        np.isclose(layer_vertices[..., 1], y0) | np.isclose(layer_vertices[..., 1], y1)
    assert on_side[:, kind == 1].all()
    assert on_side[:, inside].mean() < 0.02  # a corner reaching a side as the fibers lean is rare


def test_reconstruct_grain_angles(reconstructed, synthetic):
    # the mean fiber direction and the inclination of the grains agree with the true motion of the cells in the box
    sites = synthetic.sites_at(0.15)
    in_box = (sites[:, 0] > 0.4) & (sites[:, 0] < 0.7) & (np.abs(sites[:, 1]) < 0.15)
    moved = (synthetic.sites_at(0.3) - synthetic.sites_at(0.0))[in_box]
    expected = np.degrees(np.arctan(np.linalg.norm(moved.mean(axis=0))/0.3))
    assert abs(reconstructed.grain_angle - expected) < 1.0
    assert abs(reconstructed.grain_inclination_mean - expected) < 1.5
    assert reconstructed.grain_dispersion < 4.0
    assert reconstructed.grain_inclination_max > reconstructed.grain_inclination_mean


def test_reconstruct_measures_wall_thickness(reconstructed, synthetic):
    grains = np.load(reconstructed.output_dir / "rec-grains.npz")
    thickness = grains["wall_thickness"]
    assert 0.005 < np.percentile(thickness, 10) < 0.0075
    assert 0.0085 < np.percentile(thickness, 90) < 0.0115
    assert 0.006 < reconstructed.median_wall_thickness < 0.011
    # the connector widths vary from layer to layer (measured in every layer)
    mesh = np.loadtxt(reconstructed.output_dir / "rec-mesh.txt", skiprows=10)
    assert not np.isnan(mesh).any()
    assert mesh[:, 11].min() > 0.003 and mesh[:, 11].max() < 0.014
    assert np.linalg.norm(mesh[:, 5:8] - mesh[:, 8:11], axis=1).min() > 0


def test_reconstruct_fiber_ends_from_image(tmp_path, monkeypatch):
    # cells ending inside the volume become fiber ends (longitudinal connectors) of the grains at their corners
    monkeypatch.chdir(tmp_path)
    scan = voxelize(SMALL, VOXEL, tilt=TILT, wobble_amplitude=0.003, seed=3, fiber_length=0.4, fiber_length_cov=0.4)
    assert scan.cell_ends is not None and sum(len(e) for e in scan.cell_ends) > 0
    result = ringspy.reconstruct(scan, geoName="ends", nsegments=4, stlFlag=False, inpFlag=False)
    grains = np.load(result.output_dir / "ends-grains.npz")
    ends, z_grid = grains["fiber_ends"], grains["z_grid"]
    h = z_grid[1] - z_grid[0]
    assert len(ends) > 0 and result.nfibers > len(result.vertices)
    assert np.all((ends[:, 1] >= 1) & (ends[:, 1] <= 7))
    # every true cell end inside the box has a detected fiber end within a grid spacing (the ends are rounded to the grid)
    truth = [z for i, e in enumerate(scan.cell_ends) for z in e
             if 0.42 < scan.sites[i, 0] < 0.68 and abs(scan.sites[i, 1]) < 0.13]
    detected = np.unique(z_grid[ends[:, 1]])
    found = sum(1 for z in truth if np.any(np.abs(detected - z) <= h))
    assert len(truth) > 5 and found >= 0.8*len(truth)
    log = (result.output_dir / "ends-report.log").read_text()
    assert "fiber_ends:                              image" in log and "nfiber_ends:" in log
    mesh = np.loadtxt(result.output_dir / "ends-mesh.txt", skiprows=10)
    header = (result.output_dir / "ends-mesh.txt").read_text().splitlines()
    nbot, nreg, ntop, nlong = [int(header[i]) for i in (2, 4, 6, 8)]
    assert nlong > 0 and len(mesh) == nbot + nreg + ntop + nlong
    assert nbot == ntop and nbot > len(result.ridges)  # bot/top connectors at the faces and at every fiber end
    assert not np.isnan(mesh).any()
    # the fibers of the grains without a cell end run through the whole height: one patch
    iga = (result.output_dir / "endsIGA.txt").read_text()
    assert "Number of Patches \n{:d}".format(result.nfibers) in iga
    # other modes
    none = ringspy.reconstruct(scan, geoName="none", nsegments=2, fiber_ends="none", stlFlag=False, inpFlag=False)
    assert none.fiber_ends is None and none.layer_vertices.shape[0] == 10
    rnd = ringspy.reconstruct(scan, geoName="rnd", nsegments=3, fiber_ends="random", fiber_length=0.15, seed=4,
                              stlFlag=False, inpFlag=False)
    assert rnd.nfibers > 2*len(rnd.vertices)
    with pytest.raises(ValueError, match="fiber_length"):
        ringspy.reconstruct(scan, geoName="bad", fiber_ends="random")


def test_reconstruct_is_deterministic(synthetic, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    first = ringspy.reconstruct(synthetic, geoName="a", stlFlag=False, inpFlag=False, nsegments=1)
    (tmp_path / "again").mkdir()
    monkeypatch.chdir(tmp_path / "again")
    second = ringspy.reconstruct(synthetic, geoName="a", stlFlag=False, inpFlag=False, nsegments=1)
    assert_same_files(tmp_path / first.output_dir, tmp_path / "again" / second.output_dir)


def test_reconstruct_tracking_modes(synthetic, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    straight = ringspy.reconstruct(synthetic, geoName="straight", tracking="straight", stlFlag=False, inpFlag=False)
    assert np.all(straight.layer_vertices == straight.layer_vertices[0])
    assert straight.grain_inclination_max == 0 and straight.grain_angle == 0
    direction = ringspy.reconstruct(synthetic, geoName="direction", tracking="direction", stlFlag=False, inpFlag=False)
    assert direction.tracking_counts["matched by its cells"] == 0
    assert direction.tracking_counts["moved along the fiber direction"] > 0
    assert abs(direction.grain_inclination_mean - direction.grain_angle) < 1.5
    constant = ringspy.reconstruct(synthetic, geoName="constant", wall_thickness=0.007, stlFlag=False, inpFlag=False,
                                   boundaryFlag=False, merge_operation=False)
    mesh = np.loadtxt(constant.output_dir / "constant-mesh.txt", skiprows=10)
    np.testing.assert_allclose(mesh[:, 11], 0.007)
    vertices = np.loadtxt(constant.output_dir / "constant-vertex.mesh", skiprows=6, usecols=(0, 1))
    assert not (np.isclose(vertices[:, 0], 0.4) & np.isclose(vertices[:, 1], -0.15)).any()  # no box corners


def test_reconstruct_array_with_other_axis_order(synthetic, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    stack = np.moveaxis(synthetic.volume, 2, 0)  # a stack of cross-sections, as a TIFF stack is read
    result = ringspy.reconstruct(stack, voxel_size=VOXEL, origin=[0.4, -0.15, 0.0], geoName="stack", nsegments=1,
                                 stlFlag=False, inpFlag=False)
    assert result.longitudinal_axis == 0 and result.shape == (150, 150, 150)
    direct = ringspy.reconstruct(synthetic, geoName="direct", nsegments=1, stlFlag=False, inpFlag=False)
    assert len(result.vertices) == len(direct.vertices)
    np.testing.assert_allclose(result.vertices, direct.vertices)
    cropped = ringspy.reconstruct(synthetic, geoName="cropped", crop=[20, 120, 10, 140, 30, 120], nsegments=1,
                                  stlFlag=False, inpFlag=False)
    assert cropped.shape == (100, 130, 90)
    np.testing.assert_allclose(cropped.z_coord[0], 0.0)
    assert len(cropped.vertices) < len(direct.vertices)


@pytest.mark.slow
def test_reconstruct_grayscale_scan(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    scan = voxelize(SMALL, VOXEL, tilt=TILT, wobble_amplitude=0.003, noise=0.25, blur=0.6, seed=3)
    clean = voxelize(SMALL, VOXEL, tilt=TILT, wobble_amplitude=0.003, seed=3)
    noisy = ringspy.reconstruct(scan, geoName="noisy", smooth_sigma=1.0, stlFlag=False, inpFlag=False)
    exact = ringspy.reconstruct(clean, geoName="clean", stlFlag=False, inpFlag=False)
    assert 0.3 < noisy.threshold < 0.7
    assert abs(len(noisy.vertices) - len(exact.vertices)) <= 0.15*len(exact.vertices)
    assert abs(noisy.grain_angle - exact.grain_angle) < 0.5


@pytest.mark.slow
def test_cli_reconstruct(synthetic, tmp_path):
    np.save(tmp_path / "scan.npy", synthetic.volume)
    (tmp_path / "scan.toml").write_text('geoName = "scan"\nvolume = "scan.npy"\nvoxel_size = 0.002\n'
                                        'origin = [0.4, -0.15, 0.0]\nnsegments = 1\nstlFlag = false\n')
    subprocess.run([sys.executable, "-m", "ringspy", "reconstruct", str(tmp_path / "scan.toml"), "--geoName", "cli"],
                   cwd=tmp_path, check=True, env=dict(os.environ, MPLBACKEND="Agg"))
    out = tmp_path / "meshes" / "cli"
    assert (out / "cli-report.log").is_file() and (out / "cli.inp").is_file()
    assert "volume:" in (out / "cli-report.log").read_text()
    # the TOML file names the volume relative to its own folder
    config = ringspy.ReconstructionConfig.from_toml(tmp_path / "scan.toml")
    assert config.volume == str(tmp_path / "scan.npy")
    failed = subprocess.run([sys.executable, "-m", "ringspy", "reconstruct", str(tmp_path / "missing.toml")],
                            cwd=tmp_path, capture_output=True, text=True)
    assert failed.returncode == 2 and "error" in failed.stderr


def test_reconstruction_config_checks():
    config = ringspy.ReconstructionConfig(voxel_size=0.003, reference_slice=12, longitudinal_axis=1, crop=[0, 5, 0, 5, 0, 5])
    assert config.voxel_size == [0.003, 0.003, 0.003] and config.reference_slice == "12" and config.longitudinal_axis == "1"
    assert config.replace(nsegments=3).nsegments == 3
    assert config.to_dict()["crop"] == [0, 5, 0, 5, 0, 5]
    for bad in [dict(solid="shiny"), dict(tracking="magic"), dict(crop=[0, 1]), dict(voxel_size=[1, 2]),
                dict(voxel_size=-1), dict(longitudinal_axis="w"), dict(reference_slice="somewhere"),
                dict(nsegments=0), dict(merge_tol=0), dict(inpType="Ansys"), dict(vtk_format="xml")]:
        with pytest.raises((ValueError, TypeError)):
            ringspy.ReconstructionConfig(**bad)
    with pytest.raises(ValueError, match="did you mean"):
        ringspy.ReconstructionConfig.from_dict({"nsegment": 2})
    with pytest.raises(ValueError, match="no volume"):
        ringspy.reconstruct(ringspy.ReconstructionConfig())
    with pytest.raises(TypeError):
        ringspy.reconstruct(42)
