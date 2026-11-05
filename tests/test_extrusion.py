# -*- coding: utf-8 -*-
"""
Taper and spiral grain of the extrusion (MeshGenTools.ExtrudeLayers, GrainAngle, ConnectorTangents)
"""

import re

import numpy as np
import pytest

import ringspy
import ringspy.MeshGenTools as rpgen
from helpers import EXAMPLES_DIR

SMALL = dict(r_max=2.0, nrings=2, box_size=1.6, box_center=[0.9, 0], iter_max=300, seed=3, nsegments=3, z_max=2.4,
             stlFlag=False)


def nodes_of(inp):
    """The nodes of an Abaqus input file (x, y, z), in order"""
    block = inp.read_text().split("*Node,nset=AllNodes\n", 1)[1].split("*Element", 1)[0]
    return np.array([[float(v) for v in line.split(",")[1:4]] for line in block.strip().splitlines()])


# ------------------------------------------------------------------------------------- helpers

def test_grain_angle_profile():
    radii = np.array([0, 0.3, 1.0, 1.3, 2.0, 2.3])
    alpha = rpgen.GrainAngle(radii, [0.05, 0.1])
    np.testing.assert_allclose(alpha.knots, [0, 0.8, 1.8])  # the pith and the middles of the rings
    np.testing.assert_allclose(alpha([0, 0.4, 0.8, 1.3, 1.8, 3.0]), [0, 0.025, 0.05, 0.075, 0.1, 0.1])
    constant = rpgen.GrainAngle(radii, 0.08)
    np.testing.assert_allclose(constant([0.8, 1.8, 5.0]), 0.08)
    assert constant(0.4) == pytest.approx(0.04)
    with pytest.raises(ValueError, match="one per growth ring"):
        rpgen.GrainAngle(radii, [0.1, 0.1, 0.1])
    with pytest.raises(ValueError, match="pi/2"):
        rpgen.GrainAngle(radii, 1.6)


def test_extrude_layers_without_taper_or_grain_is_the_rigid_twist():
    rng = np.random.default_rng(0)
    points = rng.random((50, 2))*2 - 1
    theta = np.linspace(0, 0.05, 7)
    z = np.linspace(0, 3, 7)
    layers = rpgen.ExtrudeLayers(points, theta, z, (0.1, -0.2))
    np.testing.assert_array_equal(layers, rpgen._rotate_layers(points, theta, z, (0.1, -0.2))[:, :, 0:2])


def test_extrude_layers_taper_and_spiral_grain():
    rng = np.random.default_rng(0)
    points = rng.random((200, 2))*2 - 1 + [0.1, -0.2]
    theta = np.zeros(5)
    z = np.linspace(0, 2, 5)
    center = (0.1, -0.2)
    d = points - center
    r = np.linalg.norm(d, axis=1)
    # taper: the distances to the center scale linearly with z
    layers = rpgen.ExtrudeLayers(points, theta, z, center, (0, 2), stem_taper=-0.3)
    for i, zi in enumerate(z):
        np.testing.assert_allclose(np.linalg.norm(layers[i] - center, axis=1), r*(1 - 0.3*zi/2))
    np.testing.assert_allclose(layers[0], points)
    # spiral grain: a point turns by tan(alpha)*z/r (in the sense of the twist theta), its radius unchanged
    radii = np.array([0, 0.3, 1.0, 1.3, 2.0, 2.3])
    layers = rpgen.ExtrudeLayers(points, theta, z, center, (0, 2), grain_angle=0.1, radii=radii)
    alpha = rpgen.GrainAngle(radii, 0.1)(r)
    for i, zi in enumerate(z):
        di = layers[i] - center
        np.testing.assert_allclose(np.linalg.norm(di, axis=1), r)
        turned = np.arctan2(d[:, 1], d[:, 0]) - np.arctan2(di[:, 1], di[:, 0])  # the code's twist is clockwise
        turned = (turned + np.pi) % (2*np.pi) - np.pi
        np.testing.assert_allclose(turned, np.tan(alpha)*zi/r, atol=1e-9)
    # the same sense as the twist: with theta and the grain angle, the angles add up
    both = rpgen.ExtrudeLayers(points, np.full(5, 0.05), z, center, (0, 2), grain_angle=0.1, radii=radii)
    rigid = rpgen._rotate_layers(points, np.full(5, 0.05), z, center)[:, :, 0:2]
    turned = np.arctan2(rigid[-1][:, 1] - center[1], rigid[-1][:, 0] - center[0]) - np.arctan2(both[-1][:, 1] - center[1], both[-1][:, 0] - center[0])
    np.testing.assert_allclose((turned + np.pi) % (2*np.pi) - np.pi, np.tan(alpha)*2/r, atol=1e-9)
    # the pith: a finite twist (the angle rises linearly through the heart)
    pith = rpgen.ExtrudeLayers(np.array([[0.1, -0.2], [0.1 + 1e-9, -0.2]]), theta, z, center, (0, 2), grain_angle=0.1, radii=radii)
    assert np.all(np.isfinite(pith))
    with pytest.raises(ValueError, match="radii"):
        rpgen.ExtrudeLayers(points, theta, z, center, (0, 2), grain_angle=0.1)
    # with shaped rings, the grain angle follows the ring coordinate
    shape = rpgen.RingShape(center, ellipticity=1.5)
    shaped = rpgen.ExtrudeLayers(points, theta, z, center, (0, 2), grain_angle=[0.0, 0.2], radii=radii, shape=shape)
    assert not np.allclose(shaped[-1], layers[-1])


def test_connector_tangents():
    npt, nctrl = 3, 5
    base = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]])
    nodes = np.vstack([np.column_stack((base*(1 + 0.1*i), np.full(npt, 0.5*i))) for i in range(nctrl)])  # grains leaning outward
    t = rpgen.ConnectorTangents(nodes, np.array([[1, 2], [4, 5], [13, 14]]), npt, nctrl)
    np.testing.assert_allclose(np.linalg.norm(t, axis=1), 1.0)
    assert np.all(t[:, 2] > 0.9) and t[0, 0] > 0 and t[0, 1] > 0


# ------------------------------------------------------------------------------------- Config

@pytest.mark.parametrize("params,match", [
    ({"stem_taper": -1.0}, "stem_taper"),
    ({"grain_angle": [0.1, 0.1, 0.1]}, "grain_angle"),
    ({"grain_angle": 1.6}, "grain_angle"),
    ({"grain_angle": "steep"}, "grain_angle"),
])
def test_config_rejects_bad_extrusion_values(params, match):
    with pytest.raises((TypeError, ValueError), match=match):
        ringspy.Config(**params)


def test_config_extrusion_defaults():
    c = ringspy.Config(stem_taper=-0.1, grain_angle=[0.02, 0.06, 0.08, 0.05])
    assert c.grain_angle == [0.02, 0.06, 0.08, 0.05] and ringspy.Config().stem_taper == 0 and ringspy.Config().grain_angle is None
    assert ringspy.Config(grain_angle=0.05).grain_angle == 0.05


# ---------------------------------------------------------------------------------- generation

@pytest.mark.parametrize("fibers", [False, True], ids=["segments", "fibers"])
def test_generated_lattice_is_tapered_and_twisted(fibers):
    extra = dict(fiber_length=0.8) if fibers else {}
    res = ringspy.generate(dict(SMALL, geoName="spiral", stem_taper=-0.2, grain_angle=[0.05, 0.1], **extra))
    nodes = nodes_of(res.output_dir / "spiral.inp")
    n = res.nvertex
    bottom = nodes[np.isclose(nodes[:, 2], nodes[:, 2].min())]
    top = nodes[np.isclose(nodes[:, 2], nodes[:, 2].max())]
    z_top = top[0, 2]
    scale = 1 - 0.2*(z_top - 0)/2.4
    rb, rt = np.linalg.norm(bottom[:, :2], axis=1), np.linalg.norm(top[:, :2], axis=1)
    if fibers:
        # the fibers of a grain end at different heights and the cells taper to their ends, so the nodes of the
        # bottom and top layers are not the plain vertices: check the taper on the whole layers
        assert np.median(rt)/np.median(rb) == pytest.approx(scale, rel=0.05)
    else:
        bottom, top = bottom[:n], top[:n]
        np.testing.assert_allclose(bottom[:, :2], res.vertices, atol=1e-9)
        np.testing.assert_allclose(rt[:n]/rb[:n], scale, rtol=1e-6)
        # the twist on top of the rigid 0.05 (plus the longitudinal joints' share): tan(alpha(r))*z/r
        turned = np.arctan2(bottom[:, 1], bottom[:, 0]) - np.arctan2(top[:, 1], top[:, 0])
        turned = (turned + np.pi) % (2*np.pi) - np.pi
        alpha = rpgen.GrainAngle(res.radii, [0.05, 0.1])(rb[:n])
        expected = np.tan(alpha)*z_top/rb[:n] + 0.05 + 2*(0.05/3*0.02)
        inner = rb[:n] > 0.3
        np.testing.assert_allclose(turned[inner], expected[inner], atol=2e-3)
    log = (res.output_dir / "spiral-report.log").read_text()
    assert re.search(r"^stem_taper:\s+-0.2$", log, re.M) and re.search(r"^grain_angle:\s+\[0.05, 0.1\]$", log, re.M)
    # the files are complete and the connector cross-sections follow the leaning grains
    text = (res.output_dir / "spiral_conns.vtu").read_text()
    assert 'Name="Connector_width"' in text and 'Name="Ring"' in text
    assert (res.output_dir / "spiral_conns_vol.vtu").stat().st_size > 0


def test_plain_extrusion_has_no_extrusion_lines():
    res = ringspy.generate(dict(SMALL, geoName="plain", inpFlag=False))
    assert "stem_taper" not in (res.output_dir / "plain-report.log").read_text()


def test_spiral_example_loads():
    c = ringspy.load_config(EXAMPLES_DIR / "wood_spiral.toml")
    assert c.stem_taper < 0 and isinstance(c.grain_angle, list) and len(c.grain_angle) == c.nrings
