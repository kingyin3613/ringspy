import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest


@pytest.fixture(autouse=True)
def run_in_tmp_path(tmp_path, monkeypatch):
    """Write the generated meshes/ folder to a temporary directory and start every test with a clean figure"""
    monkeypatch.chdir(tmp_path)
    plt.close("all")
    yield
    plt.close("all")
