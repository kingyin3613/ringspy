# -*- coding: utf-8 -*-
"""
Config (TOML input files), the seed of the cell placement, and the command line
"""

from pathlib import Path

import numpy as np
import pytest

import ringspy
import ringspy.MeshGenTools as rpgen
from ringspy.cli import main
from ringspy.config import Config
from ringspy.pipeline import default_precrack_nodes
from helpers import EXAMPLES, EXAMPLES_DIR


def write(tmp_path, text, name="input.toml"):
    path = tmp_path / name
    path.write_text(text)
    return path


@pytest.mark.parametrize("toml", EXAMPLES, ids=[p.stem for p in EXAMPLES])
def test_examples_load(toml):
    config = ringspy.load_config(toml)
    assert config.geoName == toml.stem
    assert Path(config.sites_dir) == EXAMPLES_DIR


def test_wood_cube_example_matches_defaults():
    # the defaults are the values of examples/wood_cube.toml
    config = ringspy.load_config(EXAMPLES_DIR / "wood_cube.toml")
    defaults = Config().resolved()
    assert config.replace(seed=None, sites_dir=None) == defaults


def test_defaults_are_the_template_expressions():
    c = Config(r_min=0.5, r_max=3, nrings=3, box_size=2.5).resolved()
    assert c.width_heart == 0.3*(3-0.5)/3
    assert c.width_sparse == 0.7*(3-0.5)/3
    assert c.width_dense == 0.3*(3-0.5)/3
    assert c.z_max == 2.5
    # explicit values are kept, and replace() keeps the defaults following the other parameters
    assert Config(width_heart=0.1, z_max=1).resolved().width_heart == 0.1
    assert Config().replace(box_size=2.0).resolved().z_max == 2.0


@pytest.mark.parametrize("value,expected", [(True, True), (False, False), ("on", True), ("off", False),
                                            ("On", True), ("yes", True), ("No", False), ("TRUE", True),
                                            (np.bool_(True), True)])
def test_flags_accept_booleans_and_on_off_strings(value, expected):
    assert Config(stlFlag=value).stlFlag is expected


@pytest.mark.parametrize("value", ["onn", 1, None])
def test_bad_flag_raises(value):
    with pytest.raises(ValueError, match="stlFlag"):
        Config(stlFlag=value)


@pytest.mark.parametrize("params,match", [
    ({"box_shape": "circle"}, "box_shape"),
    ({"radial_growth_rule": "ternary"}, "radial_growth_rule"),
    ({"inpType": "LS-DYNA"}, "inpType"),
    ({"seed": -1}, "seed"),
    ({"seed": 1.5}, "seed"),
    ({"nrings": 2.0}, "nrings"),
    ({"r_max": "4"}, "r_max"),
    ({"box_center": [0, 0, 0]}, "box_center"),
    ({"cellsize_dense": 0}, "cellsize_dense"),
])
def test_bad_values_raise(params, match):
    with pytest.raises((TypeError, ValueError), match=match):
        Config(**params)


def test_unsupported_inp_type_is_fine_without_inp_file():
    assert Config(inpFlag=False, inpType="LS-DYNA").inpType == "LS-DYNA"


def test_sequences_become_tuples():
    c = Config(generation_center=[1, 2], box_center=np.array([0.5, 0.25]), precrack_nodes=[[0, 0, 1, 0]])
    assert c.generation_center == (1, 2) and c.box_center == (0.5, 0.25)
    assert c.precrack_nodes == [(0, 0, 1, 0)]
    assert c.to_dict()["generation_center"] == [1, 2]


def test_toml(tmp_path):
    path = write(tmp_path, 'geoName = "x"\nseed = 3\nstlFlag = "off"\nbox_center = [1, 0.5]\n'
                           'radial_growth_rule = "my_sites.npy"\n')
    c = ringspy.load_config(path)
    assert (c.geoName, c.seed, c.stlFlag, c.box_center) == ("x", 3, False, (1, 0.5))
    assert Path(c.sites_dir) == tmp_path  # saved sites are looked for next to the TOML file
    assert c.merge_operation is True  # default
    write(tmp_path, 'sites_dir = "data"\n', "other.toml")
    assert Path(ringspy.load_config(tmp_path / "other.toml").sites_dir) == tmp_path / "data"


def test_toml_unknown_parameter_suggests_name(tmp_path):
    path = write(tmp_path, 'geoName = "x"\nstlflag = false\n')
    with pytest.raises(ValueError, match=r"unknown parameter 'stlflag' \(did you mean 'stlFlag'\?\)"):
        ringspy.load_config(path)


def test_unknown_override_suggests_name():
    with pytest.raises(ValueError, match=r"did you mean 'box_size'"):
        Config().replace(boxsize=1.0)
    with pytest.raises(ValueError, match=r"did you mean 'seed'"):
        ringspy.generate(Config(), seeds=1)


def test_toml_absolute_sites_dir(tmp_path):
    data = (tmp_path / "data").resolve()
    path = write(tmp_path, "sites_dir = '{}'\n".format(data.as_posix()))
    assert Path(ringspy.load_config(path).sites_dir) == data


def test_toml_tables_rejected(tmp_path):
    with pytest.raises(ValueError, match=r"\[box\] tables"):
        ringspy.load_config(write(tmp_path, '[box]\nbox_size = 1.0\n'))


def test_toml_errors_name_the_file(tmp_path):
    with pytest.raises(ValueError, match="input.toml.*box_shape"):
        ringspy.load_config(write(tmp_path, 'box_shape = "circle"\n'))


def test_generate_rejects_other_types():
    with pytest.raises(TypeError):
        ringspy.generate(42)


def test_default_precrack_is_the_template_precrack():
    c = Config(box_size=4.0, box_center=(0, 0.5))
    x_min = -2.0
    # the template scripts: x_indent = x_min + box_size*0.120, x_precrack = x_indent + box_size*0.1
    assert default_precrack_nodes(c, x_min) == [[x_min + 4.0*0.120, 0.5, x_min + 4.0*0.120 + 4.0*0.1, 0.5]]


# ---------------------------------------------------------------------------
# seed / rng of the cell placement functions

PLACEMENTS = {
    "binary": lambda rng: rpgen.CellPlacement_Binary((0, 0), 0.7, 0, 2, 0.1, 0.2, 0.1, 0.04, 0.02, 200, 10**6,
                                                     rng=rng),
    "lloyd": lambda rng: rpgen.CellPlacement_Binary_Lloyd("x", "meshes", (0, 0), 0.7, 0, 2, 0.1, 0.2, 0.1,
                                                          0.04, 0.02, 3, 10**6, rng=rng),
}


@pytest.mark.parametrize("placement", PLACEMENTS.values(), ids=PLACEMENTS.keys())
def test_same_seed_same_cells(placement):
    sites1, radii1 = placement(np.random.default_rng(5))
    sites2, radii2 = placement(5)  # an int seed works too
    np.testing.assert_array_equal(sites1, sites2)
    np.testing.assert_array_equal(radii1, radii2)
    sites3, radii3 = placement(6)
    assert not np.array_equal(radii1, radii3)


@pytest.mark.parametrize("placement", PLACEMENTS.values(), ids=PLACEMENTS.keys())
def test_no_seed_uses_global_random_state(placement):
    # rng=None keeps the old behaviour: np.random.seed makes the placement repeatable
    np.random.seed(11)
    sites1, radii1 = placement(None)
    np.random.seed(11)
    sites2, radii2 = placement(None)
    np.testing.assert_array_equal(sites1, sites2)
    np.testing.assert_array_equal(radii1, radii2)
    # and a seeded placement leaves the global state alone
    np.random.seed(11)
    placement(3)
    assert np.random.random() == np.random.RandomState(11).random()


def test_flag_helpers():
    assert rpgen._is_on(True) and rpgen._is_on("on") and rpgen._is_on("Yes")
    assert not rpgen._is_on(False) and not rpgen._is_on("off") and not rpgen._is_on("ON")  # as before
    assert rpgen._flag_text(True) == "on" and rpgen._flag_text(False) == "off" and rpgen._flag_text("On") == "On"
    for shape in ["square", "triangle", "hexagon"]:
        on_bool = rpgen.Clipping_Box(shape, (0, 0), 1.0, True)
        on_str = rpgen.Clipping_Box(shape, (0, 0), 1.0, "on")
        np.testing.assert_array_equal(on_bool[5], on_str[5])


# ---------------------------------------------------------------------------
# command line

def test_cli_version(capsys):
    with pytest.raises(SystemExit) as e:
        main(["--version"])
    assert e.value.code == 0
    assert ringspy.__version__ in capsys.readouterr().out


def test_cli_needs_a_command(capsys):
    with pytest.raises(SystemExit) as e:
        main([])
    assert e.value.code == 2


def test_cli_missing_file(capsys):
    with pytest.raises(SystemExit) as e:
        main(["run", "missing.toml"])
    assert e.value.code == 2
    assert "missing.toml" in capsys.readouterr().err


def test_cli_bad_parameter(tmp_path, capsys):
    path = write(tmp_path, 'box_shape = "circle"\n')
    with pytest.raises(SystemExit) as e:
        main(["run", str(path)])
    assert e.value.code == 2
    assert "box_shape" in capsys.readouterr().err


def test_cli_overrides(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr("ringspy.pipeline.generate", lambda config: calls.append(config) or
                        type("R", (), {"output_dir": Path(config.path, config.geoName)})())
    path = write(tmp_path, 'geoName = "a"\nseed = 1\n')
    assert main(["run", str(path), "--seed", "7", "--path", "out", "--geoName", "b"]) == 0
    assert (calls[0].seed, calls[0].path, calls[0].geoName) == (7, "out", "b")
