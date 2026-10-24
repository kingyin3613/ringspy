# -*- coding: utf-8 -*-
"""
Template scripts of RingsPy <= 1.1 keep working, and ringspy.generate() gives the same files
"""

import importlib.util
from pathlib import Path

import numpy as np

import ringspy
from helpers import assert_same_files

TEMPLATE = Path(__file__).with_name("legacy_template.py")


def run_legacy_template():
    spec = importlib.util.spec_from_file_location("legacy_template", TEMPLATE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.legacy_template()


def test_legacy_template_matches_generate(tmp_path, monkeypatch):
    # the parameters of legacy_template.py; generate() fills in the ring widths and z_max as the script does
    config = ringspy.Config(geoName="legacy", r_max=1.5, nrings=2, cellsize_sparse=0.04, cellsize_dense=0.02,
                            box_size=1.0, precrackFlag="on")
    seed = 2023
    for folder, run in [("legacy", run_legacy_template), ("new", lambda: ringspy.generate(config))]:
        (tmp_path / folder).mkdir()
        monkeypatch.chdir(tmp_path / folder)
        np.random.seed(seed)  # seed=None: the cells come from the global random state, as in the template
        run()
    out = Path("meshes", "legacy")
    assert (tmp_path / "legacy" / out / "legacy.stl").is_file()
    assert "nconnector_t_precrack" in (tmp_path / "legacy" / out / "legacy-report.log").read_text()
    assert_same_files(tmp_path / "legacy" / out, tmp_path / "new" / out)
