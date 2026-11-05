# -*- coding: utf-8 -*-
"""
Hardwood features: vessels (ring-porous or diffuse-porous) and rays placed before the fibers, their
labels in the Result, the ParaView and Abaqus files, the graph and the log
"""

import re
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial import cKDTree

import ringspy
import ringspy.MeshGenTools as rpgen
from helpers import EXAMPLES_DIR
from test_rings import ring_sets, vtu_cell_data

SMALL = dict(r_max=1.5, nrings=2, box_size=1.2, box_center=[0.9, 0], iter_max=300, seed=5, stlFlag=False,
             cellsize_sparse=0.02, cellsize_dense=0.012)
RING_POROUS = dict(vessel_fraction=0.25, vessel_size=0.1, vessel_arrangement="ring", vessel_zone=0.5, ray_spacing=0.25)
DIFFUSE = dict(vessel_fraction=0.12, vessel_size=0.08, vessel_arrangement="diffuse", ray_spacing=0.2, ray_width=0.02,
               ray_cell_length=0.06)


# ----------------------------------------------------------------------------------- placement

def test_ring_radii_and_preplaced_cells_keep_the_softwood_identical():
    # the hardwood draws the same ring radii as the softwood, and preplaced cells come first in the sites
    rng = np.random.default_rng(1)
    radii = rpgen.RingRadii(0.3, 0.7, 0.3, 2, rng)
    sites, radii_b = rpgen.CellPlacement_Binary((0, 0), 2.0, 0, 2, 0.3, 0.7, 0.3, 0.02, 0.01, 200, 10**6, rng=1)
    np.testing.assert_array_equal(radii, radii_b)
    rng = np.random.default_rng(1)
    radii_drawn = rpgen.RingRadii(0.3, 0.7, 0.3, 2, rng)  # the draw CellPlacement_Binary skips when radii are given
    sites2, radii2 = rpgen.CellPlacement_Binary((0, 0), 2.0, 0, 2, 0.3, 0.7, 0.3, 0.02, 0.01, 200, 10**6,
                                                rng=rng, radii=radii_drawn)
    np.testing.assert_array_equal(sites, sites2)
    np.testing.assert_array_equal(radii2, radii)
    pre = [[0.5, 0.0, 0.05], [0.0, 0.6, 0.05]]
    sites3, _ = rpgen.CellPlacement_Binary((0, 0), 2.0, 0, 2, 0.3, 0.7, 0.3, 0.02, 0.01, 200, 10**6, rng=1, preplaced=pre)
    np.testing.assert_array_equal(sites3[:2], pre)
    d = np.linalg.norm(sites3[2:, :2, None] - np.array(pre)[:, :2].T[None], axis=1)
    assert np.all(d >= 0.05 + 0.01 - 1e-12)  # the fibers keep clear of the preplaced cells


def test_place_rays():
    radii = np.array([0, 0.3, 1.0, 1.3, 2.0, 2.3])
    cellsize = rpgen.RingProfile(radii, 0.02, 0.01)
    circles, nrays, segments = rpgen.PlaceRays((0, 0), radii, 0.25, 0.02, 0.08, cellsize, rng=0)
    circles = np.array(circles)
    segments = np.array(segments)
    assert len(segments) == len(circles) - nrays and segments.shape[1] == 5  # one between consecutive cells of a ray
    assert np.all(segments[:, 4] > 0)
    r = np.linalg.norm(circles[:, :2], axis=1)
    assert nrays >= round(2*np.pi*0.3/0.25) and len(circles) > nrays
    assert r.min() > 0.3 and r.max() < 2.3
    # the rays are radial lines: the sites of each lie on one angle
    angles = np.round(np.arctan2(circles[:, 1], circles[:, 0]), 9)
    assert len(np.unique(angles)) == nrays
    # the spacing around the circumference stays between 0.75 and 1.5 times ray_spacing (checked at 3 radii)
    for rr in [0.5, 1.2, 2.2]:
        present = np.unique(angles[np.abs(r - rr) < 0.08])
        if len(present) > 2:
            gaps = np.diff(np.concatenate((present, [present[0] + 2*np.pi])))*rr
            assert 0.6*0.25 < np.median(gaps) < 1.6*0.25
    # ray cells about ray_cell_length apart along a ray
    one = np.sort(r[angles == angles[0]])
    assert 0.06 < np.median(np.diff(one)) < 0.1
    # the fibers (and vessels) keep out of the corridors between the consecutive cells of a ray
    sites, _, site_type, _ = rpgen.CellPlacement_Hardwood((0, 0), 2.3, 0, 2, 0.3, 0.7, 0.3, 0.02, 0.01, 300, 10**6, rng=0,
                                                          ray_spacing=0.25, ray_width=0.02, ray_cell_length=0.08,
                                                          vessel_fraction=0.1, vessel_size=0.08)
    rays = sites[site_type == 2]
    angle_of = np.round(np.arctan2(rays[:, 1], rays[:, 0]), 9)
    segments = []
    for a in np.unique(angle_of):
        cells = rays[angle_of == a]
        cells = cells[np.argsort(np.linalg.norm(cells[:, :2], axis=1))]
        segments += [[*cells[k, :2], *cells[k + 1, :2], min(cells[k, 2], cells[k + 1, 2])] for k in range(len(cells) - 1)]
    index = rpgen._SegmentIndex(segments)
    assert len(segments) > 50
    assert all(index.clear(c) for c in sites[site_type != 2].tolist())
    assert not index.clear([*rays[1, :2], 0.01]) or len(rays) < 2
    with pytest.raises(ValueError):
        rpgen.PlaceRays((0, 0), radii, 0.0, 0.02, 0.08, cellsize)


@pytest.mark.parametrize("arrangement", ["ring", "diffuse"])
def test_place_vessels(arrangement):
    radii = np.array([0, 0.3, 1.0, 1.3, 2.0, 2.3])
    cellsize = rpgen.RingProfile(radii, 0.02, 0.01)
    grid = rpgen._CircleGrid(0.1)
    circles = np.array(rpgen.PlaceVessels((0, 0), radii, 0.1, 0.2, arrangement, 0.5, cellsize, grid, rng=0))
    r = np.linalg.norm(circles[:, :2], axis=1)
    band = np.searchsorted(radii, r, side="right") - 1
    assert len(circles) > 10 and np.all(band >= 1)  # no vessels in the heart
    if arrangement == "ring":
        assert np.all(band % 2 == 1)  # the sparse parts only ...
        for b in [1, 3]:
            sel = band == b
            assert sel.any() and np.all(r[sel] <= radii[b] + 0.5*(radii[b+1] - radii[b]) + 1e-12)  # ... their inner half
    else:
        assert np.any(band % 2 == 0) and np.any(band % 2 == 1)
    # the vessels keep clear of each other: at least their exclusion radii apart
    d = cKDTree(circles[:, :2]).query(circles[:, :2], k=2)[0][:, 1]
    assert d.min() >= 2*min(circles[:, 2]) - 1e-12
    # about the asked area fraction (ring: of the inner half of the sparse parts)
    if arrangement == "diffuse":
        area = np.pi*(radii[-1]**2 - radii[1]**2)
        assert 0.6 < len(circles)*np.pi*0.05**2/(0.2*area) < 1.1
    assert rpgen.PlaceVessels((0, 0), radii, 0.1, 0.0, arrangement, 0.5, cellsize, grid, rng=0) == []
    with pytest.raises(ValueError):
        rpgen.PlaceVessels((0, 0), radii, 0.1, 0.2, "clustered", 0.5, cellsize, grid)


def test_hardwood_placement_types_and_counts():
    sites, radii, site_type, info = rpgen.CellPlacement_Hardwood((0, 0), 2.0, 0, 2, 0.3, 0.7, 0.3, 0.02, 0.01, 300, 10**6,
                                                                 rng=2, vessel_fraction=0.2, vessel_size=0.1,
                                                                 vessel_arrangement="ring", ray_spacing=0.3)
    assert len(site_type) == len(sites) and set(np.unique(site_type)) == {0, 1, 2}
    assert info["nvessels"] == np.sum(site_type == 1) > 0
    assert info["nray_cells"] == np.sum(site_type == 2) > info["nrays"] > 0
    # the rays and vessels come first, then the fibers
    assert np.all(site_type[:info["nvessels"] + info["nray_cells"]] != 0) and np.all(site_type[info["nvessels"] + info["nray_cells"]:] == 0)
    # the same ring radii as the softwood with the same seed
    _, radii_soft = rpgen.CellPlacement_Binary((0, 0), 2.0, 0, 2, 0.3, 0.7, 0.3, 0.02, 0.01, 300, 10**6, rng=2)
    np.testing.assert_array_equal(radii, radii_soft)
    # vessel sites keep the fibers at least half a vessel away (the walls halfway: vessel_size wide cells)
    vessels, fibers = sites[site_type == 1], sites[site_type == 0]
    d = cKDTree(fibers[:, :2]).query(vessels[:, :2])[0]
    assert d.min() >= 0.05 + 0.01 - 1e-12


# ------------------------------------------------------------------------------------- Config

def test_config_hardwood_parameters():
    c = ringspy.Config(vessel_fraction=0.2, vessel_size=0.1, vessel_arrangement="ring", ray_spacing=0.3)
    assert c.ray_width == 0.02 and c.ray_cell_length == 0.1 and ringspy.Config().vessel_fraction == 0
    assert ringspy.Config().ray_spacing is None


@pytest.mark.parametrize("params,match", [
    ({"vessel_fraction": 1.0}, "vessel_fraction"),
    ({"vessel_fraction": -0.1}, "vessel_fraction"),
    ({"vessel_size": 0}, "vessel_size"),
    ({"vessel_arrangement": "clustered"}, "vessel_arrangement"),
    ({"vessel_zone": 0}, "vessel_zone"),
    ({"ray_spacing": 0}, "ray_spacing"),
    ({"ray_width": -1}, "ray_width"),
    ({"ray_cell_length": 0}, "ray_cell_length"),
    ({"ray_spacing": 0.3, "radial_growth_rule": "binary_lloyd"}, "binary"),
    ({"vessel_fraction": 0.1, "radial_growth_rule": "regular_hexagonal"}, "binary"),
])
def test_config_rejects_bad_hardwood_values(params, match):
    with pytest.raises(ValueError, match=match):
        ringspy.Config(**params)


# ------------------------------------------------------------------------------ generation

@pytest.fixture(scope="module")
def ring_porous(tmp_path_factory):
    import os
    workdir = tmp_path_factory.mktemp("ring_porous")
    cwd = os.getcwd()
    os.chdir(workdir)
    try:
        result = ringspy.generate(dict(SMALL, geoName="oak", graphFlag=True, save_sites=True, **RING_POROUS))
    finally:
        os.chdir(cwd)
    result.output_dir = workdir / result.output_dir
    return result


def test_result_labels(ring_porous):
    res = ring_porous
    assert len(res.site_type) == res.nParticles and set(np.unique(res.site_type)) == {0, 1, 2}
    for name in ["vertex_vessel", "vertex_ray"]:
        assert getattr(res, name).dtype == bool and len(getattr(res, name)) == res.nvertex
    for name in ["ridge_vessel", "ridge_ray"]:
        assert getattr(res, name).dtype == bool and len(getattr(res, name)) == res.nridge
    assert res.ridge_vessel.sum() > 0 and res.ridge_ray.sum() > 0 and not np.any(res.ridge_vessel & res.ridge_ray)
    # a vessel wall has a vessel site among the two cells nearest to its midpoint
    mid = res.vertices[res.ridges].mean(axis=1)
    nearest = cKDTree(res.sites[:, :2]).query(mid, k=2)[1]
    types = res.site_type[nearest]
    np.testing.assert_array_equal(res.ridge_vessel, np.any(types == 1, axis=1))
    np.testing.assert_array_equal(res.ridge_ray, np.any(types == 2, axis=1) & ~res.ridge_vessel)
    # a vertex of a vessel has a vessel wall
    has_vessel_wall = np.zeros(res.nvertex, dtype=bool)
    has_vessel_wall[res.ridges[res.ridge_vessel].ravel()] = True
    np.testing.assert_array_equal(res.vertex_vessel, has_vessel_wall)
    # vessels in the inner half of the earlywood only
    r = np.linalg.norm(res.sites[res.site_type == 1, :2], axis=1)
    band = np.searchsorted(res.radii, r, side="right") - 1
    assert np.all(band % 2 == 1)
    # the walls of a vessel enclose a cell about vessel_size wide: the vessel site is far from its walls
    d = cKDTree(mid[res.ridge_vessel]).query(res.sites[res.site_type == 1, :2])[0]
    assert np.median(d) > 0.3*RING_POROUS["vessel_size"]
    # the cell types are saved with the sites
    assert Path(res.output_dir, "oak_types.npy").is_file()
    np.testing.assert_array_equal(np.load(Path(res.output_dir, "oak_types.npy")), res.site_type)


def test_files_carry_vessel_and_ray_labels(ring_porous):
    res = ring_porous
    out = res.output_dir
    # ParaView: Vessel and Ray cell data after Ring and Dense, consistent with the Result
    for suffix in ["_vertices", "_beams", "_conns", "_conns_vol"]:
        text = (out / ("oak" + suffix + ".vtu")).read_text()
        names = re.findall(r'<DataArray type="Int32" Name="(\w+)" format="ascii">', text)
        assert names[-4:] == ["Ring", "Dense", "Vessel", "Ray"], suffix
    _, vessel = vtu_cell_data(out / "oak_conns.vtu", "Vessel")
    _, ray = vtu_cell_data(out / "oak_conns.vtu", "Ray")
    nt = res.nridge*3*2  # transverse connectors: 3 per segment, 2 segments
    np.testing.assert_array_equal(vessel[:nt], np.tile(res.ridge_vessel, 6))
    np.testing.assert_array_equal(ray[:nt], np.tile(res.ridge_ray, 6))
    _, vessel_beams = vtu_cell_data(out / "oak_beams.vtu", "Vessel")
    np.testing.assert_array_equal(vessel_beams[:res.nvertex], res.vertex_vessel)
    # Abaqus: VesselBeams/TConns/LConns and RayBeams/TConns/LConns, next to the ring sets
    inp = out / "oak.inp"
    text = inp.read_text()
    sets = {m.group(1): [int(x) for x in m.group(2).replace(",", " ").split()]
            for m in re.finditer(r"\*Elset, elset=((?:Vessel|Ray)(?:Beams|TConns|LConns))\n((?:[\d, ]+\n)*)", text)}
    assert set(sets) == {"VesselBeams", "VesselTConns", "VesselLConns", "RayBeams", "RayTConns", "RayLConns"}
    assert len(sets["VesselBeams"]) == res.vertex_vessel.sum()*2*2  # 2 elements per beam, 2 segments
    assert len(sets["VesselTConns"]) == res.ridge_vessel.sum()*6
    assert not set(sets["VesselTConns"]) & set(sets["RayTConns"])
    assert ring_sets(inp)["Ring1TConns"]  # the ring sets are still there
    assert text.index("** Element sets of the vessel cells") < text.index("** Section: Section-1")
    # the log
    log = (out / "oak-report.log").read_text()
    assert "HARDWOOD:" in log
    for key in ["vessel_fraction", "vessel_arrangement", "ray_spacing", "nvessels", "nrays", "nray_cells", "nvessel_walls"]:
        assert re.search(r"^{}:\s+\S".format(key), log, re.M), key
    assert int(re.search(r"^nvessels:\s+(\d+)", log, re.M).group(1)) == np.sum(res.site_type == 1)
    # the graph: vessel and ray features
    g = np.load(out / "oak-graph.npz")
    assert list(g["node_features"][-2:]) == ["vessel", "ray"] and list(g["edge_features"][-2:]) == ["vessel", "ray"]
    np.testing.assert_array_equal(g["x"][:, -2], res.vertex_vessel)
    np.testing.assert_array_equal(g["edge_attr"][:res.nridge, -1], res.ridge_ray)


def test_softwood_has_no_hardwood_labels_or_sections():
    res = ringspy.generate(dict(SMALL, geoName="soft", inpFlag=False, graphFlag=True, save_sites=True))
    assert not res.site_type.any() and not res.ridge_vessel.any() and not res.vertex_ray.any()
    text = (res.output_dir / "soft_conns.vtu").read_text()
    assert 'Name="Vessel"' not in text and 'Name="Ray"' not in text
    assert "HARDWOOD" not in (res.output_dir / "soft-report.log").read_text()
    assert not (res.output_dir / "soft_types.npy").exists()
    assert res.graph["x"].shape[1] == 8 and not res.graph["x"][:, 6:].any()


def test_diffuse_porous_and_saved_types(tmp_path):
    first = ringspy.generate(dict(SMALL, geoName="maple", inpFlag=False, save_sites=True, **DIFFUSE))
    r = np.linalg.norm(first.sites[first.site_type == 1, :2], axis=1)
    band = np.searchsorted(first.radii, r, side="right") - 1
    assert np.any(band % 2 == 0) and np.any(band % 2 == 1) and np.all(band >= 1)
    # the saved sites and types mesh the same hardwood again
    second = ringspy.generate(dict(SMALL, geoName="maple2", inpFlag=False, radial_growth_rule="maple_sites.npy",
                                   sites_dir=str(first.output_dir)))
    np.testing.assert_array_equal(second.sites, first.sites)
    np.testing.assert_array_equal(second.site_type, first.site_type)
    np.testing.assert_array_equal(second.ridge_vessel, first.ridge_vessel)
    assert "HARDWOOD:" in (second.output_dir / "maple2-report.log").read_text()


def test_hardwood_with_fibers_of_various_lengths():
    res = ringspy.generate(dict(SMALL, geoName="oak_fibers", fiber_length=0.6, nsegments=3, z_max=1.2, **RING_POROUS))
    assert res.ridge_vessel.any() and res.ridge_ray.any()
    text = (res.output_dir / "oak_fibers_conns.vtu").read_text()
    assert 'Name="Vessel"' in text and 'Name="Ray"' in text
    inp = (res.output_dir / "oak_fibers.inp").read_text()
    assert "*Elset, elset=VesselTConns" in inp and "*Elset, elset=RayBeams" in inp
    _, vessel = vtu_cell_data(res.output_dir / "oak_fibers_conns.vtu", "Vessel")
    assert vessel.sum() > 0


def test_hardwood_examples_load():
    for name in ["hardwood_ring_porous", "hardwood_diffuse"]:
        c = ringspy.load_config(EXAMPLES_DIR / (name + ".toml"))
        assert c.vessel_fraction > 0 and c.ray_spacing is not None
    assert ringspy.load_config(EXAMPLES_DIR / "hardwood_ring_porous.toml").vessel_arrangement == "ring"
    assert ringspy.load_config(EXAMPLES_DIR / "hardwood_diffuse.toml").vessel_arrangement == "diffuse"
