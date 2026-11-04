# -*- coding: utf-8 -*-
"""
The graph file (graphFlag) and the batch generation (ringspy.batch, ringspy batch)
"""

import csv
import json
from pathlib import Path

import numpy as np
import pytest

import ringspy
from ringspy.batch import parse_seeds, parse_set
from ringspy.cli import main
from ringspy.pipeline import EDGE_FEATURES, NODE_FEATURES, graph_arrays
from helpers import EXAMPLES_DIR

SMALL = dict(r_max=1.0, nrings=2, box_size=0.8, box_center=[0.4, 0], iter_max=200, seed=4, stlFlag=False, inpFlag=False)


def test_graph_file(tmp_path):
    res = ringspy.generate(dict(SMALL, geoName="graph", graphFlag=True, path=str(tmp_path)))
    file = res.output_dir / "graph-graph.npz"
    assert file.is_file() and file in res.files
    g = np.load(file)
    assert set(g.files) >= {"pos", "x", "edge_index", "edge_attr", "node_features", "edge_features", "sites", "radii",
                            "generation_center", "nParticles", "nbeamElem", "mass", "bulk_volume", "bulk_density",
                            "porosity", "config"}
    n, e = res.nvertex, res.nridge
    assert g["pos"].shape == (n, 2) and g["x"].shape == (n, len(NODE_FEATURES))
    assert g["edge_index"].shape == (2, 2*e) and g["edge_index"].dtype == np.int64
    assert g["edge_attr"].shape == (2*e, len(EDGE_FEATURES))
    assert list(g["node_features"]) == NODE_FEATURES and list(g["edge_features"]) == EDGE_FEATURES
    # the same arrays as the Result
    for key, value in res.graph.items():
        np.testing.assert_array_equal(g[key], value)
    np.testing.assert_array_equal(g["pos"], res.vertices)
    # every wall in both directions, and the node features
    np.testing.assert_array_equal(g["edge_index"][:, :e], res.ridges.T)
    np.testing.assert_array_equal(g["edge_index"][:, e:], res.ridges[:, ::-1].T)
    x = g["x"]
    np.testing.assert_array_equal(x[:, 0], res.vertex_ring)
    np.testing.assert_array_equal(x[:, 1], res.vertex_dense)
    np.testing.assert_array_equal(x[:, 2], res.vertex_boundary)
    np.testing.assert_array_equal(x[:, 3], res.vertex_thickness)
    np.testing.assert_array_equal(x[:, 4], np.bincount(res.ridges.ravel(), minlength=n))
    np.testing.assert_allclose(x[:, 5], np.linalg.norm(res.vertices, axis=1))
    assert np.all(x[:, 4] >= 2)
    # edge features: length and angle of each direction, thickness and labels of the wall
    a = g["edge_attr"]
    vec = res.vertices[res.ridges[:, 1]] - res.vertices[res.ridges[:, 0]]
    np.testing.assert_allclose(a[:e, 0], np.linalg.norm(vec, axis=1))
    np.testing.assert_allclose(a[e:, 0], a[:e, 0])
    np.testing.assert_allclose(a[:e, 1], np.arctan2(vec[:, 1], vec[:, 0]))
    np.testing.assert_allclose(np.cos(a[e:, 1]), -np.cos(a[:e, 1]), atol=1e-12)
    np.testing.assert_array_equal(a[:e, 2], res.vertex_thickness[res.ridges[:, 1]])
    np.testing.assert_array_equal(a[:e, 3], res.ridge_ring)
    np.testing.assert_array_equal(a[:e, 4], res.ridge_dense)
    np.testing.assert_array_equal(a[:e, 5], res.ridge_boundary)
    # model properties and the parameters
    assert float(g["porosity"]) == res.porosity and int(g["nParticles"]) == res.nParticles
    assert json.loads(str(g["config"]))["geoName"] == "graph"
    np.testing.assert_array_equal(g["sites"], res.sites)
    log = (res.output_dir / "graph-report.log").read_text()
    assert "Graph File:" in log and "graph-graph.npz" in log


def test_graph_arrays_without_the_flag():
    res = ringspy.generate(dict(SMALL, geoName="nograph"))
    assert res.graph is None and not (res.output_dir / "nograph-graph.npz").exists()
    g = graph_arrays(res)
    assert g["x"].shape == (res.nvertex, len(NODE_FEATURES))


def test_batch_sweeps_and_seeds(tmp_path):
    base = ringspy.Config(**dict(SMALL, geoName="sweep", path=str(tmp_path / "out")))
    results = ringspy.batch(base, seeds=[1, 2], sweep={"cellsize_sparse": [0.02, 0.03], "theta_max": [0.0]},
                            graphFlag=True)
    assert len(results) == 4 and all(r is not None for r in results)
    assert [r.config.geoName for r in results] == ["sweep_000", "sweep_001", "sweep_002", "sweep_003"]
    assert [r.config.seed for r in results] == [1, 2, 1, 2]
    assert [r.config.cellsize_sparse for r in results] == [0.02, 0.02, 0.03, 0.03]
    assert all(r.config.theta_max == 0.0 and r.config.graphFlag for r in results)
    assert all((r.output_dir / (r.config.geoName + "-graph.npz")).is_file() for r in results)
    # the same seed and parameters give the same cells, in a batch as in generate()
    np.testing.assert_array_equal(results[0].sites, results[2].sites[:0].shape and
                                  ringspy.generate(base, seed=1, cellsize_sparse=0.02, theta_max=0.0, geoName="single").sites)
    # the table
    with open(tmp_path / "out" / "sweep-batch.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 4
    assert list(rows[0]) == ["run", "geoName", "seed", "cellsize_sparse", "theta_max", "nParticles", "nvertex", "nridge",
                             "nbeamElem", "mass", "bulk_volume", "bulk_density", "porosity", "output_dir", "error"]
    assert [int(r["run"]) for r in rows] == [0, 1, 2, 3]
    assert [r["seed"] for r in rows] == ["1", "2", "1", "2"]
    assert [float(r["cellsize_sparse"]) for r in rows] == [0.02, 0.02, 0.03, 0.03]
    for row, res in zip(rows, results):
        assert int(row["nParticles"]) == res.nParticles and int(row["nvertex"]) == res.nvertex
        assert float(row["porosity"]) == res.porosity and row["error"] == ""
        assert Path(row["output_dir"]) == res.output_dir


def test_batch_records_failed_runs_and_list_values(tmp_path):
    base = dict(SMALL, geoName="robust", path=str(tmp_path / "out"))
    results = ringspy.batch(base, sweep={"radial_growth_rule": ["binary", "missing_sites.npy"],
                                         "cellwallthickness_dense": [[0.01, 0.02]]})
    assert results[0] is not None and results[1] is None
    with open(tmp_path / "out" / "robust-batch.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    assert rows[0]["error"] == "" and "missing_sites.npy" in rows[1]["error"] and rows[1]["nParticles"] == ""
    assert json.loads(rows[0]["cellwallthickness_dense"]) == [0.01, 0.02]
    assert results[0].config.cellwallthickness_dense == [0.01, 0.02]


def test_batch_parallel(tmp_path):
    base = dict(SMALL, geoName="par", path=str(tmp_path / "out"))
    results = ringspy.batch(base, seeds=range(3), jobs=2)
    assert all(r is not None for r in results) and len(results) == 3
    assert all((r.output_dir / (r.config.geoName + "-report.log")).is_file() for r in results)


def test_batch_rejects_bad_input():
    with pytest.raises(ValueError, match="unknown parameter"):
        ringspy.batch(SMALL, sweep={"cell_size_sparse": [0.02]})
    with pytest.raises(ValueError, match="non-empty list"):
        ringspy.batch(SMALL, sweep={"cellsize_sparse": []})
    with pytest.raises(ValueError, match="jobs"):
        ringspy.batch(SMALL, jobs=0)
    with pytest.raises(TypeError):
        ringspy.batch(42)


def test_parse_seeds_and_set():
    assert parse_seeds(["0-2", "7", "10,11"]) == [0, 1, 2, 7, 10, 11]
    with pytest.raises(ValueError, match="seeds"):
        parse_seeds(["a"])
    assert parse_set(["cellsize_sparse=0.02,0.03", "stlFlag=false", "box_shape=square,hexagon",
                      'inpType="Abaqus"', "cellsize_dense=[0.01,0.02],[0.02,0.02]"]) == {
        "cellsize_sparse": [0.02, 0.03], "stlFlag": [False], "box_shape": ["square", "hexagon"],
        "inpType": ["Abaqus"], "cellsize_dense": [[0.01, 0.02], [0.02, 0.02]]}
    with pytest.raises(ValueError, match="--set"):
        parse_set(["cellsize_sparse"])
    with pytest.raises(ValueError, match="at least one value"):
        parse_set(["cellsize_sparse="])


def test_cli_batch(tmp_path):
    toml = tmp_path / "small.toml"
    toml.write_text("".join("{} = {}\n".format(k, json.dumps(v) if not isinstance(v, bool) else str(v).lower())
                            for k, v in dict(SMALL, geoName="cli").items()))
    assert main(["batch", str(toml), "--seeds", "1-2", "--set", "theta_max=0,0.05", "--path", str(tmp_path / "m"),
                 "--geoName", "clibatch"]) == 0
    with open(tmp_path / "m" / "clibatch-batch.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    assert [r["geoName"] for r in rows] == ["clibatch_000", "clibatch_001", "clibatch_002", "clibatch_003"]
    assert [r["theta_max"] for r in rows] == ["0", "0", "0.05", "0.05"]
    assert (tmp_path / "m" / "clibatch_003" / "clibatch_003-report.log").is_file()
    # a failing run makes the command exit with 1, after the other runs
    assert main(["batch", str(toml), "--set", "radial_growth_rule=binary,nothing.npy", "--path", str(tmp_path / "f")]) == 1
    assert (tmp_path / "f" / "cli_000" / "cli_000-report.log").is_file()
    with pytest.raises(SystemExit):
        main(["batch", str(toml), "--seeds", "x"])
    with pytest.raises(SystemExit):
        main(["batch", str(toml), "--set", "no_such_parameter=1"])
