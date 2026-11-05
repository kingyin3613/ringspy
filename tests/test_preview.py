# -*- coding: utf-8 -*-
"""
The quick 2D preview (previewFlag, ringspy preview, ringspy run --preview): the same cell mesh as the
full run, no extrusion, no mesh files; and the figures of plot_mesh
"""

import csv

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import ringspy  # noqa: E402
from ringspy.cli import main  # noqa: E402
from ringspy.plot import COLORINGS, plot_mesh  # noqa: E402
from helpers import EXAMPLES_DIR  # noqa: E402

SMALL = dict(r_max=1.0, nrings=2, box_size=0.8, box_center=[0.4, 0], iter_max=200, seed=4, stlFlag=False, inpFlag=False)
MESH_FIELDS = ["sites", "radii", "vertices", "ridges", "vertex_ring", "vertex_dense", "vertex_boundary",
               "vertex_thickness", "ridge_ring", "ridge_dense", "ridge_boundary", "site_type", "vertex_vessel",
               "vertex_ray", "ridge_vessel", "ridge_ray", "boundary_points"]


def test_preview_same_mesh_as_full_run_without_the_3d_model(tmp_path):
    full = ringspy.generate(dict(SMALL, geoName="full", path=str(tmp_path), statsFlag=True))
    quick = ringspy.preview(dict(SMALL, geoName="quick", path=str(tmp_path), statsFlag=True))
    assert quick.preview and not full.preview
    assert quick.config.previewFlag and not full.config.previewFlag
    assert quick.nParticles == full.nParticles
    for name in MESH_FIELDS:
        np.testing.assert_array_equal(getattr(quick, name), getattr(full, name), err_msg=name)
    for name in ["nbeamElem", "mass", "bulk_volume", "bulk_density", "porosity"]:
        assert getattr(quick, name) is None and getattr(full, name) is not None
    assert quick.summary() == full.summary()
    for stats in ["cell_stats", "ring_stats", "wall_stats"]:
        for key, value in getattr(full, stats).items():
            np.testing.assert_array_equal(getattr(quick, stats)[key], value, err_msg=stats + " " + key)
    # the preview writes the 2D mesh files, its figure and the statistics, nothing of the 3D model
    names = sorted(f.name for f in quick.files)
    assert names == ["quick-cells.csv", "quick-ridge.mesh", "quick-rings.csv", "quick-vertex.mesh",
                     "quick-walls.csv", "quick.png"]
    assert (tmp_path / "quick" / "quick.png").stat().st_size > 10000
    for suffix in ["-vertex.mesh", "-ridge.mesh"]:
        assert (tmp_path / "quick" / ("quick" + suffix)).read_bytes() == \
            (tmp_path / "full" / ("full" + suffix)).read_bytes()
    assert not any(f.suffix in [".vtu", ".txt", ".log", ".stl", ".inp"] for f in quick.files)


def test_preview_other_files_and_paths(tmp_path):
    # the sites and the graph can be saved from a preview; the graph has no 3D properties
    res = ringspy.preview(dict(SMALL, geoName="p", path=str(tmp_path), save_sites=True, graphFlag=True))
    names = sorted(f.name for f in res.files)
    assert names == ["p-graph.npz", "p-ridge.mesh", "p-vertex.mesh", "p.png", "p_radii.npy", "p_sites.npy"]
    graph = np.load(tmp_path / "p" / "p-graph.npz")
    assert "nParticles" in graph.files and "pos" in graph.files
    assert not any(k in graph.files for k in ["nbeamElem", "mass", "bulk_volume", "bulk_density", "porosity"])
    assert int(graph["nParticles"]) == res.nParticles
    # fibers of various lengths, hardwood, a shaped ring, a polygon: the preview stops before any of the extrusions
    fibers = ringspy.preview(dict(SMALL, geoName="f", path=str(tmp_path), fiber_length=0.3, nsegments=2))
    assert fibers.preview and fibers.nvertex > 10
    hardwood = ringspy.preview(dict(SMALL, geoName="h", path=str(tmp_path), vessel_fraction=0.1, vessel_size=0.1,
                                    ray_spacing=0.3, ring_eccentricity=0.3))
    assert hardwood.ridge_vessel.any() and hardwood.ridge_ray.any()
    assert (tmp_path / "h" / "h.png").is_file()
    polygon = ringspy.preview(dict(SMALL, geoName="d", path=str(tmp_path), box_shape="polygon",
                                   box_polygon=[[0, -0.4], [0.8, -0.4], [0.8, 0.4], [0.4, 0.1], [0, 0.4]]))
    assert len(polygon.boundary_points) == 5
    # the figure of the preview is closed afterwards
    assert plt.get_fignums() == []


def toml_text(params):
    return "\n".join("{} = {}".format(k, str(v).lower() if isinstance(v, bool) else
                                      '"{}"'.format(v) if isinstance(v, str) else v) for k, v in params.items())


def test_preview_from_toml_and_cli(tmp_path):
    toml = tmp_path / "p.toml"
    toml.write_text(toml_text(dict(SMALL, geoName="p", previewFlag=True)))
    config = ringspy.Config.from_toml(toml)
    assert config.previewFlag
    res = ringspy.generate(config, path=str(tmp_path / "m"))
    assert res.preview and (tmp_path / "m" / "p" / "p.png").is_file()
    # --preview on the command line, with a file that does not ask for one
    toml.write_text(toml_text(dict(SMALL, geoName="c")))
    assert main(["run", str(toml), "--preview", "--path", str(tmp_path / "cli"), "--seed", "2"]) == 0
    files = sorted(f.name for f in (tmp_path / "cli" / "c").iterdir())
    assert files == ["c-ridge.mesh", "c-vertex.mesh", "c.png"]
    # a batch of previews: the table has the statistics but no mass
    assert main(["batch", str(toml), "--preview", "--seeds", "1-2", "--path", str(tmp_path / "b"),
                 "--geoName", "pb"]) == 0
    with open(tmp_path / "b" / "pb-batch.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    assert [r["geoName"] for r in rows] == ["pb_000", "pb_001"]
    assert all(r["mass"] == "" and r["nbeamElem"] == "" and float(r["porosity_2d"]) > 0.5 for r in rows)
    assert all(int(r["nvertex"]) > 10 for r in rows)
    assert not (tmp_path / "b" / "pb_000" / "pb_000-report.log").exists()


@pytest.mark.parametrize("color", COLORINGS)
def test_plot_mesh_colorings(tmp_path, color):
    res = ringspy.preview(dict(SMALL, geoName="h", path=str(tmp_path), vessel_fraction=0.1, vessel_size=0.1,
                               ray_spacing=0.3))
    fig, ax = plt.subplots()
    out = plot_mesh(res, ax=ax, color=color, sites=True, linewidth=0.8 if color == "dense" else "thickness")
    assert out is ax
    collections = ax.collections
    assert len(collections) >= 2  # the walls and the sites
    walls = collections[0]
    assert len(walls.get_segments()) == res.nridge
    widths = np.asarray(walls.get_linewidths(), dtype=float)
    if color == "dense":
        assert np.allclose(widths, 0.8)
    else:
        assert len(widths) == res.nridge and widths.max() > widths.min()
    # the vessel walls are drawn in their colour in every colouring (types=True)
    colors = walls.get_colors()
    if len(colors) == res.nridge:
        vessel = np.asarray(res.ridge_vessel, dtype=bool)
        assert np.all(colors[vessel] == colors[vessel][0])
        assert not np.all(colors[~vessel] == colors[vessel][0])
    assert ax.get_aspect() == 1.0
    assert ax.get_title().startswith("h: ")
    plt.close(fig)


def test_plot_mesh_from_a_full_result_and_options(tmp_path):
    res = ringspy.generate(dict(SMALL, geoName="full", path=str(tmp_path)))
    ax = res.plot(rings=False, boundary=False, legend=False, title=False)
    assert ax.get_title() == "" and ax.get_legend() is None
    assert len(ax.patches) == 0 and len(ax.lines) == 0  # no outline, no ring boundaries
    plt.close(ax.figure)
    ax = res.plot()
    assert len(ax.patches) == 1 and len(ax.lines) == len(res.radii) - 1 and ax.get_legend() is not None
    plt.close(ax.figure)
    with pytest.raises(ValueError):
        plot_mesh(res, color="rings")


def test_preview_of_the_examples_is_the_full_mesh(tmp_path):
    for name in ["wood_cube", "dogbone"]:
        full = ringspy.generate(EXAMPLES_DIR / (name + ".toml"), path=str(tmp_path / "full"), stlFlag=False,
                                inpFlag=False)
        quick = ringspy.preview(EXAMPLES_DIR / (name + ".toml"), path=str(tmp_path / "quick"))
        np.testing.assert_array_equal(quick.vertices, full.vertices)
        np.testing.assert_array_equal(quick.ridges, full.ridges)
        np.testing.assert_array_equal(quick.ridge_ring, full.ridge_ring)
