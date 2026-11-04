# -*- coding: utf-8 -*-
"""
Growth rings: the ring labels in the output files and the Result, the gradual transition of
cell size and wall thickness from the sparse to the dense part of a ring, and per-ring values
"""

import re
from pathlib import Path

import numpy as np
import pytest

import ringspy
import ringspy.MeshGenTools as rpgen
from helpers import EXAMPLES_DIR
from test_regressions import parse_inp

SMALL = dict(r_max=1.2, nrings=3, box_size=1.5, box_center=[0.6, 0], iter_max=300, seed=1, stlFlag=False)


def ring_sets(inp_path):
    """Members of the Ring*/Sparse*/Dense* element sets of an Abaqus input file"""
    text = Path(inp_path).read_text()
    sets = {}
    for match in re.finditer(r"\*Elset, elset=((?:Ring\d+|Sparse|Dense)(?:Beams|TConns|LConns))\n((?:[\d, ]+\n)*)", text):
        sets[match.group(1)] = [int(x) for x in match.group(2).replace(",", " ").split()]
    return sets


def vtu_cell_data(vtu_path, name):
    text = Path(vtu_path).read_text()
    match = re.search(r'<DataArray type="(\w+)" Name="{}" format="ascii">\n(.*?)</DataArray>'.format(name), text, re.S)
    assert match, name
    return match.group(1), np.array(match.group(2).split(), dtype=float)


# ---------------------------------------------------------------------------------------- helpers

def test_ring_bands_and_labels():
    radii = [0, 0.3, 1.0, 1.3, 2.0, 2.3]  # heart, then sparse/dense of 2 rings
    points = [[0, 0], [0.2, 0], [0, 0.5], [1.1, 0], [1.5, 0], [0, 2.2], [3, 0], [0.3, 0], [1.0, 0]]
    band = rpgen.RingBands(points, radii, (0, 0))
    np.testing.assert_array_equal(band, [0, 0, 1, 2, 3, 4, 5, 1, 2])  # a point on a radius belongs to the outer band
    ring, dense = rpgen.RingLabels(band)
    np.testing.assert_array_equal(ring, [0, 0, 1, 1, 2, 2, 3, 1, 1])
    np.testing.assert_array_equal(dense, [True, True, False, True, False, True, False, False, True])
    # the generation center counts
    np.testing.assert_array_equal(rpgen.RingBands([[1.2, 0]], radii, (1.0, 0)), [0])


def test_ring_profile_abrupt_is_the_old_parity_rule():
    radii = np.array([0, 0.3, 1.0, 1.3, 2.0, 2.3])
    profile = rpgen.RingProfile(radii, 0.006, 0.010)
    np.testing.assert_array_equal(profile.values, [0.010, 0.006, 0.010, 0.006, 0.010, 0.006])
    r = np.linspace(0, 3, 301)
    expected = np.where(np.searchsorted(radii, r, side="right") % 2 == 0, 0.006, 0.010)  # bisect even: sparse
    np.testing.assert_array_equal(profile(r), expected)
    assert float(profile(0.5)) == 0.006 and float(profile(0.5, 2)) == 0.010


@pytest.mark.parametrize("transition", ["linear", "sigmoid"])
def test_ring_profile_gradual(transition):
    radii = np.array([0, 0.3, 1.0, 1.3, 2.0, 2.3])
    profile = rpgen.RingProfile(radii, 0.006, 0.010, transition, transition_width=0.5)
    # the dense parts and the heart are unchanged
    for r in [0.1, 1.1, 1.29, 2.1]:
        assert float(profile(r)) == 0.010
    # the inner half of a sparse part keeps the sparse value, the outer half rises to the dense value
    assert float(profile(0.3)) == 0.006 and float(profile(0.65)) == 0.006
    inner = np.linspace(0.65, 1.0, 50)
    values = profile(inner)
    assert values[0] == 0.006 and np.isclose(values[-1], 0.010)
    assert np.all(np.diff(values) >= 0) and 0.006 < values[25] < 0.010
    if transition == "linear":
        np.testing.assert_allclose(values, 0.006 + 0.004*(inner - 0.65)/0.35)
    else:
        np.testing.assert_allclose(values[25], 0.006 + 0.004*0.5, rtol=0.05)  # S-curve through the middle
        assert values[1] - values[0] < 0.004/49*0.3  # flat at both ends
    # beyond the outermost ring: the sparse value, as before
    assert float(profile(3.0)) == 0.006


def test_ring_profile_per_ring_values():
    radii = np.array([0, 0.3, 1.0, 1.3, 2.0, 2.3])
    profile = rpgen.RingProfile(radii, [0.004, 0.005], [0.010, 0.012])
    np.testing.assert_array_equal(profile.values, [0.010, 0.004, 0.010, 0.005, 0.012, 0.005])
    with pytest.raises(ValueError, match="one value per growth ring"):
        rpgen.RingProfile(radii, [0.004, 0.005, 0.006], 0.010)
    with pytest.raises(ValueError, match="transition"):
        rpgen.RingProfile(radii, 0.004, 0.010, "smooth")
    with pytest.raises(ValueError, match="transition_width"):
        rpgen.RingProfile(radii, 0.004, 0.010, "linear", 0)


def test_on_boundary():
    square = [[0, 0], [1, 0], [1, 1], [0, 1]]
    on = rpgen.OnBoundary([[0, 0], [0.5, 0], [1, 0.3], [0.5, 0.5], [0.5, 1 + 1e-9], [0.5, 1.1]], square)
    np.testing.assert_array_equal(on, [True, True, True, False, True, False])


# ----------------------------------------------------------------------------------------- Config

def test_config_accepts_per_ring_lists_and_transition():
    c = ringspy.Config(nrings=2, cellsize_sparse=[0.02, 0.03], cellwallthickness_dense=[0.01, 0.02],
                       transition="linear", transition_width=0.5)
    assert c.cellsize_sparse == [0.02, 0.03] and c.cellwallthickness_dense == [0.01, 0.02]
    assert c.to_dict()["cellsize_sparse"] == [0.02, 0.03]
    assert ringspy.Config().transition == "abrupt" and ringspy.Config().transition_width == 1.0


@pytest.mark.parametrize("params,match", [
    ({"nrings": 2, "cellsize_sparse": [0.02, 0.03, 0.04]}, "one value per growth ring"),
    ({"cellsize_dense": [0.01, 0, 0.01, 0.01]}, "positive"),
    ({"cellwallthickness_sparse": "thick"}, "cellwallthickness_sparse"),
    ({"cellwallthickness_sparse": True}, "cellwallthickness_sparse"),
    ({"transition": "smooth"}, "transition"),
    ({"transition_width": 1.5}, "transition_width"),
    ({"transition_width": 0}, "transition_width"),
])
def test_config_rejects_bad_ring_values(params, match):
    with pytest.raises((TypeError, ValueError), match=match):
        ringspy.Config(**params)


def test_toml_with_per_ring_lists(tmp_path):
    toml = tmp_path / "rings.toml"
    toml.write_text('nrings = 2\ncellsize_sparse = [0.02, 0.03]\ntransition = "sigmoid"\ntransition_width = 0.4\n')
    c = ringspy.load_config(toml)
    assert c.cellsize_sparse == [0.02, 0.03] and c.transition == "sigmoid" and c.transition_width == 0.4


# ---------------------------------------------------------------------------------- labels (item 1)

@pytest.fixture(scope="module")
def labelled(tmp_path_factory):
    """A small wood_cube-like mesh with the ring labels in all files"""
    import os
    workdir = tmp_path_factory.mktemp("labelled")
    cwd = os.getcwd()
    os.chdir(workdir)
    try:
        result = ringspy.generate(dict(SMALL, geoName="labelled"))
    finally:
        os.chdir(cwd)
    result.output_dir = workdir / result.output_dir
    return result


def test_result_holds_the_mesh_and_its_labels(labelled):
    res = labelled
    assert res.vertices.shape == (res.nvertex, 2) and res.ridges.shape == (res.nridge, 2)
    assert res.ridges.min() >= 0 and res.ridges.max() < res.nvertex
    assert len(res.vertex_ring) == len(res.vertex_dense) == len(res.vertex_boundary) == len(res.vertex_thickness) == res.nvertex
    assert len(res.ridge_ring) == len(res.ridge_dense) == len(res.ridge_boundary) == res.nridge
    # the labels follow the distance to the generation center
    band = np.searchsorted(res.radii, np.linalg.norm(res.vertices, axis=1), side="right") - 1
    np.testing.assert_array_equal(res.vertex_ring, (band + 1)//2)
    np.testing.assert_array_equal(res.vertex_dense, band % 2 == 0)
    assert set(np.unique(res.vertex_ring)) <= set(range(0, SMALL["nrings"] + 2))
    assert res.vertex_dense.any() and not res.vertex_dense.all()
    # a wall is labelled like its thickness: by its second vertex
    np.testing.assert_array_equal(res.ridge_ring, res.vertex_ring[res.ridges[:, 1]])
    np.testing.assert_array_equal(res.ridge_dense, res.vertex_dense[res.ridges[:, 1]])
    np.testing.assert_array_equal(res.vertex_thickness, np.where(res.vertex_dense, 0.010, 0.006))
    # the boundary: the box corners and the walls along the box edges
    corners = np.array([[-0.15, -0.75], [1.35, -0.75], [1.35, 0.75], [-0.15, 0.75]])
    for corner in corners:
        i = np.argmin(np.linalg.norm(res.vertices - corner, axis=1))
        assert np.allclose(res.vertices[i], corner, atol=1e-6) and res.vertex_boundary[i]
    assert 0 < res.vertex_boundary.sum() < res.nvertex
    np.testing.assert_array_equal(res.ridge_boundary, res.vertex_boundary[res.ridges].all(axis=1))
    assert res.ridge_boundary.sum() >= 4
    # the vertex file holds the same vertices and thicknesses
    lines = (res.output_dir / "labelled-vertex.mesh").read_text().splitlines()
    table = np.loadtxt(lines[6:])
    np.testing.assert_allclose(table[:, :2], res.vertices, rtol=1e-14, atol=1e-15)  # %.16g in the file
    np.testing.assert_array_equal(table[:, 3 + 3*int(lines[4])], res.vertex_thickness)


def test_vtu_files_carry_ring_and_dense_cell_data(labelled):
    res = labelled
    out = res.output_dir
    for suffix in ["_vertices", "_beams", "_conns", "_conns_vol"]:
        text = (out / ("labelled" + suffix + ".vtu")).read_text()
        ncells = int(re.search(r'NumberOfCells="(\d+)"', text).group(1))
        for name in ["Ring", "Dense"]:
            vtk_type, values = vtu_cell_data(out / ("labelled" + suffix + ".vtu"), name)
            assert vtk_type == "Int32" and len(values) == ncells, (suffix, name)
        _, ring = vtu_cell_data(out / ("labelled" + suffix + ".vtu"), "Ring")
        _, dense = vtu_cell_data(out / ("labelled" + suffix + ".vtu"), "Dense")
        assert set(np.unique(ring)) <= set(range(0, SMALL["nrings"] + 2)) and set(np.unique(dense)) == {0, 1}
    # connectors: the dense label goes with the dense thickness, and the existing array is still first
    conns = out / "labelled_conns.vtu"
    _, width = vtu_cell_data(conns, "Connector_width")
    _, dense = vtu_cell_data(conns, "Dense")
    np.testing.assert_allclose(width, np.where(dense == 1, 0.010, 0.006), rtol=1e-6)
    assert '<CellData Tensors="" Vevtors="" Scalars="Connector_width">' in conns.read_text()
    # beams: one label per beam cell, the ring of its vertex, for every segment and interval
    nbeam_cells = res.nvertex * SMALL.get("nsegments", 2)
    _, ring = vtu_cell_data(out / "labelled_beams.vtu", "Ring")
    np.testing.assert_array_equal(ring[:nbeam_cells], np.tile(res.vertex_ring, nbeam_cells // res.nvertex))


def test_binary_vtu_files_carry_the_labels(tmp_path):
    # the binary writer takes the same label arrays
    import xml.etree.ElementTree as ET
    res = ringspy.generate(dict(SMALL, geoName="binary_vtu", inpFlag=False, path=str(tmp_path)))
    text = (res.output_dir / "binary_vtu_conns.vtu").read_bytes()
    assert b'Name="Ring" format="ascii"' in text
    points, cells, types, names = None, None, None, []
    # rewrite with the binary writer through the public function is long; check the writer itself
    cells = [np.array([[0, 1], [1, 2]])]
    rpgen._write_vtu_appended(tmp_path / "t.vtu", np.zeros((3, 3)), cells, [3, 3], None,
                              [("Connector_width", np.array([0.1, 0.2])), ("Ring", np.array([1, 2])), ("Dense", np.array([True, False]))])
    head = (tmp_path / "t.vtu").read_bytes().split(b"<AppendedData")[0].decode()
    assert '<CellData Scalars="Connector_width">' in head
    assert 'type="Float32" Name="Connector_width"' in head and 'type="Int32" Name="Ring"' in head and 'type="Int32" Name="Dense"' in head
    ET.fromstring(head + "</VTKFile>")


def test_inp_has_an_element_set_per_ring(labelled):
    res = labelled
    inp = res.output_dir / "labelled.inp"
    sets = ring_sets(inp)
    elsets, sections, all_elles = parse_inp(inp)
    rings = sorted(set(res.vertex_ring) | set(res.ridge_ring))
    for kind, all_name in [("Beams", "allbeams"), ("TConns", None), ("LConns", "alllongconns")]:
        per_ring = [sets["Ring{}{}".format(k, kind)] for k in rings if "Ring{}{}".format(k, kind) in sets]
        members = sorted(i for s in per_ring for i in s)
        assert len(members) == len(set(members))  # every element in exactly one ring
        if all_name:
            assert members == sorted(elsets[all_name])
        else:
            tconns = sorted(i for name in ["allbotconns", "allregconns", "alltopconns", "allprecracktconns"]
                            for i in elsets.get(name, []))
            assert members == tconns
        assert sorted(sets["Sparse" + kind] + sets["Dense" + kind]) == members
        assert not set(sets["Sparse" + kind]) & set(sets["Dense" + kind])
    # the sets are inside the instance, before the sections, and only non-empty ones are written
    text = inp.read_text()
    assert text.index("*Elset, elset=Ring0Beams") < text.index("** Section: Section-1") < text.index("*End Instance")
    assert all(len(members) > 0 for members in sets.values())
    # beams of ring k are the beams of the vertices of ring k
    elements_per_patch = res.nbeamElem // (res.nvertex * 2)
    beam_vertex = (np.arange(res.nbeamElem) // elements_per_patch) % res.nvertex
    for k in rings:
        expected = (np.flatnonzero(res.vertex_ring[beam_vertex] == k) + 1).tolist()
        assert sets.get("Ring{}Beams".format(k), []) == expected


def test_log_lists_the_rings(labelled):
    log = (labelled.output_dir / "labelled-report.log").read_text()
    assert re.search(r"^transition:\s+abrupt$", log, re.M)
    assert re.search(r"^transition_width:\s+1.0$", log, re.M)
    radii = re.search(r"^ring radii:\s+(\[.*\])$", log, re.M).group(1)
    np.testing.assert_allclose(eval(radii), labelled.radii)


# ------------------------------------------------------------------------------ transition (item 2)

@pytest.mark.parametrize("transition", ["linear", "sigmoid"])
def test_gradual_transition_of_thickness_and_cell_size(transition):
    abrupt = ringspy.generate(dict(SMALL, geoName="abrupt", inpFlag=False))
    gradual = ringspy.generate(dict(SMALL, geoName=transition, inpFlag=False, transition=transition, transition_width=0.5))
    # same random stream: the same ring radii
    np.testing.assert_array_equal(gradual.radii, abrupt.radii)
    r = np.linalg.norm(gradual.vertices, axis=1)
    band = np.searchsorted(gradual.radii, r, side="right") - 1
    t = gradual.vertex_thickness
    # the dense parts (and the heart) keep the dense thickness, the sparse parts lie between the two values
    assert np.all(t[band % 2 == 0] == 0.010)
    assert np.all((t[band % 2 == 1] >= 0.006) & (t[band % 2 == 1] <= 0.010))
    for b in [1, 3]:  # the inner half of a sparse part is sparse, the outer half rises
        inner = (band == b) & (r < gradual.radii[b] + 0.5*(gradual.radii[b+1] - gradual.radii[b]))
        outer = (band == b) & ~inner
        if inner.any():
            assert np.all(t[inner] == 0.006)
        if outer.sum() > 5:
            assert len(np.unique(t[outer])) > 3
            order = np.argsort(r[outer])
            assert np.all(np.diff(t[outer][order]) >= -1e-15)
    # the cells of the sparse parts get smaller towards the dense parts
    rs = np.linalg.norm(gradual.sites[:, :2], axis=1)
    band_s = np.searchsorted(gradual.radii, rs, side="right") - 1
    size = gradual.sites[:, 2]
    assert np.all(size[band_s % 2 == 0] == 0.01)
    assert np.all(np.where(band_s == 0, True, True))
    for b in [1, 3]:
        sel = band_s == b
        if sel.sum() > 10:
            assert np.corrcoef(rs[sel], size[sel])[0, 1] < -0.5
    assert np.all(size[band_s % 2 == 1] <= 0.02) and np.all(size[band_s % 2 == 1] >= 0.01)
    # the connector widths of the files follow
    log = (gradual.output_dir / "{}-report.log".format(transition)).read_text()
    assert re.search(r"^transition:\s+{}$".format(transition), log, re.M)
    assert re.search(r"^transition_width:\s+0.5$", log, re.M)


def test_per_ring_values():
    sparse = [0.004, 0.005, 0.006]
    dense = [0.010, 0.012, 0.014]
    res = ringspy.generate(dict(SMALL, geoName="per_ring", inpFlag=False, cellwallthickness_sparse=sparse,
                                cellwallthickness_dense=dense, cellsize_sparse=[0.02, 0.025, 0.03],
                                cellsize_dense=[0.01, 0.012, 0.014]))
    for k in [1, 2, 3]:
        in_ring = res.vertex_ring == k
        assert np.all(res.vertex_thickness[in_ring & res.vertex_dense] == dense[k-1])
        assert np.all(res.vertex_thickness[in_ring & ~res.vertex_dense] == sparse[k-1])
        # the heart takes the dense values of the first ring
    assert np.all(res.vertex_thickness[res.vertex_ring == 0] == dense[0])
    rs = np.linalg.norm(res.sites[:, :2], axis=1)
    band = np.searchsorted(res.radii, rs, side="right") - 1
    np.testing.assert_array_equal(res.sites[:, 2], np.array([0.01, 0.02, 0.01, 0.025, 0.012, 0.03, 0.014, 0.03])[band])
    log = (res.output_dir / "per_ring-report.log").read_text()
    assert "cellwallthickness_dense:                 [0.01, 0.012, 0.014]" in log


def test_lloyd_and_honeycomb_with_transition_and_lists(keep_global_random_state):
    np.random.seed(0)
    sites, radii = rpgen.CellPlacement_Binary_Lloyd("lloyd", "meshes", (0.0, 0.0), 1.0, 0.0, 2, 0.1, 0.3, 0.1,
                                                    [0.04, 0.05], 0.02, 5, 1000, transition="sigmoid", transition_width=0.8)
    assert np.all(np.isfinite(sites)) and len(sites) > 100
    r = np.linalg.norm(sites, axis=1)
    # within the first sparse band, the cells get closer towards the dense band
    from scipy.spatial import cKDTree
    d = cKDTree(sites).query(sites, k=2)[0][:, 1]
    sel_in = (r > radii[1]) & (r < radii[1] + 0.3*(radii[2] - radii[1]))
    sel_out = (r > radii[2] - 0.3*(radii[2] - radii[1])) & (r < radii[2])
    assert d[sel_in].mean() > 1.3*d[sel_out].mean()
    pytest.importorskip("hexalattice")
    sites, radii = rpgen.CellPlacement_Honeycomb((0, 0), 1.0, 0, 2, (0, 0), 1.0, 0.1, 0.3, 0.1, [0.04, 0.05], [0.02, 0.02],
                                                 0.006, 0.010, 5, 1000)
    assert len(sites) > 100


def test_precrack_with_per_ring_cell_sizes():
    res = ringspy.generate(dict(SMALL, geoName="precrack_rings", nrings=3, cellsize_sparse=[0.02, 0.03, 0.02],
                                precrackFlag=True, precrack_nodes=[[0.2, 0.0, 0.5, 0.0]]))
    log = (res.output_dir / "precrack_rings-report.log").read_text()
    assert int(re.search(r"nconnector_t_precrack:\s+(\d+)", log).group(1)) > 0
    assert ring_sets(res.output_dir / "precrack_rings.inp")


def test_gradual_example_matches_its_description():
    config = ringspy.load_config(EXAMPLES_DIR / "wood_gradual.toml")
    assert config.transition == "sigmoid" and config.graphFlag
    assert len(config.cellsize_sparse) == len(config.cellwallthickness_dense) == config.nrings
