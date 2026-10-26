# -*- coding: utf-8 -*-
"""
End-to-end runs of the example configurations in examples/
"""

import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import ringspy
import ringspy.MeshGenTools as rpgen
from helpers import EXAMPLES, OUTPUT_SUFFIXES, ROOT, assert_same_files, inside_polygon

pytestmark = pytest.mark.slow
examples = pytest.mark.parametrize("toml", EXAMPLES, ids=[p.stem for p in EXAMPLES])


def test_examples_found():
    assert {p.stem for p in EXAMPLES} >= {"wood_cube", "wood_cube_lloyd", "honeycomb_hex", "honeycomb_tri",
                                          "saved_sites", "dogbone"}


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    """generated(toml) -> Result of ringspy.generate(toml), run once per example in a temporary folder"""
    results = {}

    def run(toml):
        if toml not in results:
            workdir = tmp_path_factory.mktemp(toml.stem)
            cwd = os.getcwd()
            os.chdir(workdir)
            try:
                results[toml] = ringspy.generate(toml)
            finally:
                os.chdir(cwd)
            results[toml].output_dir = workdir / results[toml].output_dir
        return results[toml]

    return run


@examples
def test_example_writes_all_files(toml, generated):
    result = generated(toml)
    config = result.config
    out = result.output_dir
    assert out.name == config.geoName and out.parent.name == config.path
    expected = [config.geoName + suffix for suffix in OUTPUT_SUFFIXES]
    if config.stlFlag:
        expected.append(config.geoName + ".stl")
    if config.inpFlag:
        expected.append(config.geoName + ".inp")
    assert sorted(p.name for p in out.iterdir()) == sorted(expected)
    for name in expected:
        assert (out / name).stat().st_size > 0, name

    if config.stlFlag:
        # binary STL: 80-byte header, triangle count, 50 bytes per triangle
        stl = (out / (config.geoName + ".stl")).read_bytes()
        ntriangles = int.from_bytes(stl[80:84], "little")
        assert ntriangles > 0 and len(stl) == 84 + 50 * ntriangles

    log = (out / (config.geoName + "-report.log")).read_text()
    assert re.search(r"^seed:\s+{}$".format(config.seed), log, re.M)
    assert re.search(r"^nParticles:\s+{}$".format(result.nParticles), log, re.M)
    assert result.nParticles == len(result.sites) > 0
    assert 0 < result.porosity < 1


@examples
def test_example_vertices_inside_box(toml, generated):
    result = generated(toml)
    config = result.config
    vertices = np.loadtxt(result.output_dir / (config.geoName + "-vertex.mesh"), skiprows=6, usecols=(0, 1))
    assert len(vertices) > 0
    *_, polygon, _ = rpgen.Clipping_Box(config.box_shape, config.box_center, config.box_size, config.boundaryFlag,
                                        polygon=config.box_polygon)
    assert inside_polygon(vertices, polygon, tol=1e-6).all()
    if not config.boundaryFlag:
        return
    # the mesh fills the box: its vertices include the box corners
    x_min, x_max, y_min, y_max = polygon[:, 0].min(), polygon[:, 0].max(), polygon[:, 1].min(), polygon[:, 1].max()
    np.testing.assert_allclose(vertices.min(axis=0), [x_min, y_min], atol=1e-6)
    np.testing.assert_allclose(vertices.max(axis=0), [x_max, y_max], atol=1e-6)


@examples
def test_example_same_seed_same_files(toml, generated, tmp_path):
    # a second run, through the command line, gives the same files
    result = generated(toml)
    assert result.config.seed is not None or result.config.radial_growth_rule.endswith(".npy")
    subprocess.run([sys.executable, "-m", "ringspy", "run", str(toml)], cwd=tmp_path, check=True,
                   env=dict(os.environ, MPLBACKEND="Agg"))
    assert_same_files(result.output_dir, tmp_path / result.config.path / result.config.geoName)


def test_different_seed_different_cells():
    config = ringspy.load_config(ROOT / "examples" / "wood_cube.toml").replace(
        r_max=1.0, nrings=1, box_size=0.5, stlFlag=False, inpFlag=False)
    results = [ringspy.generate(config, seed=seed, geoName="seed{}_{}".format(seed, i))
               for i, seed in enumerate([1, 1, 2])]
    np.testing.assert_array_equal(results[0].sites, results[1].sites)
    assert results[0].sites.shape != results[2].sites.shape or \
        not np.array_equal(results[0].sites, results[2].sites)


def test_example_scripts_run_their_toml(tmp_path):
    # the thin scripts in examples/ load the TOML file next to them
    script = ROOT / "examples" / "honeycomb_tri.py"
    subprocess.run([sys.executable, str(script)], cwd=tmp_path, check=True, env=dict(os.environ, MPLBACKEND="Agg"))
    assert (tmp_path / "meshes" / "honeycomb_tri" / "honeycomb_tri-report.log").is_file()


def test_save_sites_then_reuse(tmp_path):
    # save_sites = true writes the cells; a *_sites.npy growth rule meshes them again
    base = ringspy.load_config(ROOT / "examples" / "honeycomb_hex.toml").replace(stlFlag=False, inpFlag=False)
    first = ringspy.generate(base, geoName="first", save_sites=True)
    saved = Path(first.output_dir, "first_sites.npy")
    assert saved.is_file() and Path(first.output_dir, "first_radii.npy").is_file()
    second = ringspy.generate(base, geoName="second", radial_growth_rule="first_sites.npy",
                              sites_dir=str(first.output_dir))
    np.testing.assert_array_equal(second.sites, first.sites)
    for suffix in ["-vertex.mesh", "-ridge.mesh", "-mesh.txt", "IGA.txt"]:
        assert Path(first.output_dir, "first" + suffix).read_bytes() == \
            Path(second.output_dir, "second" + suffix).read_bytes(), suffix
