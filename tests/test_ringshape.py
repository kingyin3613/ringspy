# -*- coding: utf-8 -*-
"""
Shape of the growth rings: ellipses, an eccentric pith (rings wider on one side) and wavy rings
(MeshGenTools.RingShape), in the placement, the wall thickness and the ring labels
"""

import re

import numpy as np
import pytest

import ringspy
import ringspy.MeshGenTools as rpgen
from helpers import EXAMPLES_DIR

SMALL = dict(r_max=1.5, nrings=2, box_size=2.4, box_center=[0, 0], iter_max=300, seed=7, stlFlag=False, inpFlag=False)
SHAPE = dict(ring_ellipticity=1.4, ring_ellipse_angle=0.4, ring_eccentricity=0.25, ring_eccentricity_angle=-1.2,
             ring_waviness=0.04, ring_wave_order=4, ring_wave_phase=0.3)


# ------------------------------------------------------------------------------------ RingShape

def test_circular_shape_is_the_distance_to_the_center():
    s = rpgen.RingShape((0.5, -0.25))
    assert s.circular and s.scale_max == 1.0 and s.scale_min == 1.0
    points = np.array([[0.5, -0.25], [1.5, -0.25], [0.5, 0.75], [3.5, 3.75]])
    np.testing.assert_array_equal(s.rho(points), np.linalg.norm(points - [0.5, -0.25], axis=1))
    np.testing.assert_array_equal(s.scale(np.array([0.0, 1.0, 2.0])), 1.0)
    assert s.circumference(2.0) == 2*np.pi*2.0
    x, y = s.point(1.0, np.pi/2)
    np.testing.assert_allclose([x, y], [0.5, 0.75])
    assert s.inside_band([1.0, -0.25, 0.1], 0.3, 0.7) and not s.inside_band([1.1, -0.25, 0.1], 0.3, 0.7)


def test_shape_roundtrip_and_boundaries():
    s = rpgen.RingShape((0.2, 0.1), **{k.replace("ring_", ""): v for k, v in SHAPE.items()})
    assert not s.circular
    rho = np.array([0.3, 1.0, 2.5])
    theta = np.array([0.1, 2.0, 5.0])
    x, y = s.point(rho, theta)
    np.testing.assert_allclose(s.rho(np.column_stack((x, y))), rho)
    # the boundary at rho is the curve rho(point) == rho, closed
    b = s.boundary(1.3, npoints=100)
    assert b.shape == (101, 2) and np.allclose(b[0], b[-1])
    np.testing.assert_allclose(s.rho(b), 1.3)
    # the circumference: an ellipse's for the plain ellipse (Ramanujan's approximation)
    e = rpgen.RingShape(ellipticity=2.0)
    a, bb = 1.0, 0.5
    np.testing.assert_allclose(e.circumference(1.0), np.pi*(3*(a + bb) - np.sqrt((3*a + bb)*(a + 3*bb))), rtol=1e-4)
    assert e.circumference(3.0) == 3*e.circumference(1.0)


def test_ellipse_eccentricity_and_waviness_geometry():
    e = rpgen.RingShape(ellipticity=2.0, ellipse_angle=np.pi/2)  # long axis along y
    x, y = e.point(1.0, np.array([0.0, np.pi/2]))
    np.testing.assert_allclose([x[0], y[1]], [0.5, 1.0])  # short semi-axis 1/2 along x, long semi-axis 1 along y
    ecc = rpgen.RingShape(eccentricity=0.3, eccentricity_angle=np.pi)  # wide towards -x
    np.testing.assert_allclose(ecc.scale(np.array([np.pi, 0.0, np.pi/2])), [1.3, 0.7, 1.0])
    w = rpgen.RingShape(waviness=0.1, wave_order=3, wave_phase=0.5)
    theta = np.linspace(0, 2*np.pi, 3601)
    s = w.scale(theta)
    crests = theta[np.flatnonzero((s[1:-1] > s[:-2]) & (s[1:-1] > s[2:])) + 1]
    assert len(crests) == 3 and np.isclose(crests[0], 0.5/3, atol=1e-3) and np.isclose(s.max(), 1.1)
    # the heart and the rings keep their ring coordinate: a band of the shape is a band of rho
    radii = np.array([0, 0.3, 1.0, 1.3, 2.0, 2.3])
    pts = np.column_stack(w.point(np.array([0.1, 0.5, 1.1, 2.5]), np.array([0.2, 1.0, 3.0, 4.0])))
    np.testing.assert_array_equal(rpgen.RingBands(pts, radii, (0, 0), shape=w), [0, 1, 2, 5])


@pytest.mark.parametrize("params", [dict(ellipticity=0.5), dict(eccentricity=-0.1), dict(waviness=0.5, eccentricity=0.6),
                                    dict(wave_order=0), dict(wave_order=2.5)])
def test_bad_shape_values(params):
    with pytest.raises(ValueError):
        rpgen.RingShape(**params)


# --------------------------------------------------------------------------------------- Config

def test_config_ring_shape():
    c = ringspy.Config(**SHAPE)
    s = c.ring_shape()
    assert isinstance(s, rpgen.RingShape) and s.ellipticity == 1.4 and s.wave_order == 4 and not s.circular
    assert ringspy.Config().ring_shape().circular
    assert ringspy.Config(generation_center=[1, 2]).ring_shape().generation_center == (1.0, 2.0)


@pytest.mark.parametrize("params,match", [
    ({"ring_ellipticity": 0.8}, "ring_ellipticity"),
    ({"ring_eccentricity": -0.1}, "ring_eccentricity"),
    ({"ring_eccentricity": 0.6, "ring_waviness": 0.5}, "less than 1"),
    ({"ring_wave_order": 0}, "ring_wave_order"),
    ({"ring_wave_order": 2.5}, "ring_wave_order"),
])
def test_config_rejects_bad_shape_values(params, match):
    with pytest.raises((TypeError, ValueError), match=match):
        ringspy.Config(**params)


# ----------------------------------------------------------------------------------- generation

def bands_of(res):
    return np.searchsorted(res.radii, res.ring_shape.rho(res.vertices), side="right") - 1


@pytest.mark.parametrize("rule", ["binary", "binary_lloyd"])
def test_cells_thickness_and_labels_follow_the_shape(rule):
    extra = dict(radial_growth_rule=rule, iter_max=20) if rule == "binary_lloyd" else {}
    res = ringspy.generate(dict(SMALL, geoName="shaped_" + rule, **SHAPE, **extra))
    s = res.ring_shape
    assert not s.circular and s.ellipticity == 1.4
    # the ring labels and the wall thickness are those of the ring coordinate, not of the distance
    band = bands_of(res)
    np.testing.assert_array_equal(res.vertex_ring, (band + 1)//2)
    np.testing.assert_array_equal(res.vertex_dense, band % 2 == 0)
    np.testing.assert_array_equal(res.vertex_thickness, np.where(res.vertex_dense, 0.010, 0.006))
    distance_band = np.searchsorted(res.radii, np.linalg.norm(res.vertices, axis=1), side="right") - 1
    assert np.mean(distance_band != band) > 0.1  # the shape matters
    # the cells sit in the bands they were placed in: the sites of the dense bands are small
    site_band = np.searchsorted(res.radii, s.rho(res.sites[:, :2]), side="right") - 1
    if rule == "binary":
        np.testing.assert_array_equal(res.sites[:, 2], np.where(site_band % 2 == 0, 0.01, 0.02))
    # the rings are wider towards the eccentricity angle than opposite: more cells there
    theta = np.arctan2(res.sites[:, 1], res.sites[:, 0])
    towards = np.cos(theta - SHAPE["ring_eccentricity_angle"])
    inner = (site_band >= 1) & (site_band <= 2)
    assert np.sum(inner & (towards > 0.7)) > 1.3*np.sum(inner & (towards < -0.7))
    log = (res.output_dir / ("shaped_" + rule + "-report.log")).read_text()
    assert "RING SHAPE:" in log and re.search(r"^ring_ellipticity:\s+1.4$", log, re.M)


def test_circular_defaults_are_unchanged():
    # same cells and files as before the ring shape: the circular shape takes the old code paths
    a = ringspy.generate(dict(SMALL, geoName="plain"))
    b = ringspy.generate(dict(SMALL, geoName="explicit", ring_ellipticity=1.0, ring_eccentricity=0.0, ring_waviness=0.0))
    np.testing.assert_array_equal(a.sites, b.sites)
    assert a.ring_shape.circular
    assert "RING SHAPE" not in (a.output_dir / "plain-report.log").read_text()
    sites, radii = rpgen.CellPlacement_Binary((0, 0), 2.0, 0, 2, 0.3, 0.7, 0.3, 0.02, 0.01, 200, 10**6, rng=1)
    sites2, _ = rpgen.CellPlacement_Binary((0, 0), 2.0, 0, 2, 0.3, 0.7, 0.3, 0.02, 0.01, 200, 10**6, rng=1,
                                           shape=rpgen.RingShape((0, 0)))
    np.testing.assert_array_equal(sites, sites2)


def test_hardwood_follows_the_shape():
    res = ringspy.generate(dict(SMALL, geoName="shaped_oak", vessel_fraction=0.2, vessel_size=0.1,
                                vessel_arrangement="ring", ray_spacing=0.3, **SHAPE))
    s = res.ring_shape
    vessels = res.sites[res.site_type == 1, :2]
    band = np.searchsorted(res.radii, s.rho(vessels), side="right") - 1
    assert len(vessels) > 5 and np.all(band % 2 == 1)  # in the earlywood of the shaped rings
    # the rays follow the directions of the shape: the ray cells of one ray share their angle around the center
    rays = res.sites[res.site_type == 2, :2]
    angles = np.round(np.arctan2(rays[:, 1], rays[:, 0]), 9)
    assert len(np.unique(angles)) < len(rays)/3
    rho = s.rho(rays)
    assert rho.min() > res.radii[1] and rho.max() < res.radii[-1]


def test_synthetic_thickness_uses_the_shape():
    from ringspy.synthetic import _resolve_geometry
    c = ringspy.Config(**dict(SMALL, geoName="syn", **SHAPE))
    sites, thickness, box, center = _resolve_geometry(c, np.random.default_rng(7))
    band = np.searchsorted(ringspy.generate(dict(SMALL, geoName="syn_ref", **SHAPE)).radii, c.ring_shape().rho(sites[:, :2]),
                           side="right") - 1
    np.testing.assert_array_equal(thickness, np.where(band % 2 == 0, 0.010, 0.006))


def test_eccentric_example_loads():
    c = ringspy.load_config(EXAMPLES_DIR / "wood_eccentric.toml")
    assert c.ring_eccentricity > 0 and c.ring_ellipticity > 1 and c.ring_waviness > 0
