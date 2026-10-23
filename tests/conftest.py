import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest


@pytest.fixture(autouse=True)
def run_in_tmp_path(tmp_path, monkeypatch):
    """Write the generated meshes/ folder to a temporary directory and start every test with a clean figure"""
    monkeypatch.chdir(tmp_path)
    plt.close("all")
    yield
    plt.close("all")


@pytest.fixture
def keep_global_random_state():
    """For tests that seed np.random: restore its state afterwards, the template tests use it unseeded"""
    state = np.random.get_state()
    yield
    np.random.set_state(state)
