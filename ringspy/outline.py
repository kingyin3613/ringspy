# -*- coding: utf-8 -*-
"""
Clipping outlines from drawing files: the polygon of box_shape = "polygon" read from an SVG or a
DXF file, flattened curves, units converted to mm, and simplified with the Douglas-Peucker
algorithm so that no edge is shorter than the merge tolerance

author: king_yin3613
email: haoyin2022@u.northwestern.edu
"""

import math
import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from .MeshGenTools import cross2D, polygon_vertices

OUTLINE_SUFFIXES = ['.svg', '.dxf']
CURVE_SEGMENTS = 24      # straight segments per curve (Bezier, arc) when no tolerance is given
SVG_UNITS_MM = {'mm': 1.0, 'cm': 10.0, 'm': 1000.0, 'in': 25.4, 'pt': 25.4/72, 'pc': 25.4/6, 'px': 25.4/96, '': None}
DXF_UNITS_MM = {0: None, 1: 25.4, 2: 304.8, 4: 1.0, 5: 10.0, 6: 1000.0, 8: 0.0254, 9: 0.0254, 10: 914.4,
                11: 1e-7, 12: 1e-6, 13: 1e-3, 14: 100.0, 15: 10000.0}


def read_outline(file, tolerance=None, scale=1.0, center=None):
    """
    The clipping polygon in a drawing file, as an (n, 2) array of vertices in mm, counter-clockwise.

    file :: an SVG file (the first of its paths, polygons, polylines, rects, circles and ellipses that
        encloses the largest area; a path's subpath of largest area; transforms of the element and
        its groups applied; curves flattened; the y axis flipped so that the drawing is not mirrored)
        or a DXF file (its closed LWPOLYLINE or POLYLINE of largest area, or a closed loop of LINE and
        ARC entities, or a CIRCLE; bulges and arcs flattened).
    tolerance :: Douglas-Peucker simplification tolerance (mm): the outline keeps within this
        distance of the drawn one, with as few vertices as that allows (None: no simplification,
        curves flattened into CURVE_SEGMENTS straight segments each). Use the merge tolerance of the
        mesh, so that no edge of the polygon is shorter than it.
    scale :: factor applied to the coordinates after the file's units are converted to mm (an SVG
        without units is taken in px at 96 dpi, a DXF without $INSUNITS as mm).
    center :: (x, y): the outline is moved so that the center of its bounding box is here.
    """
    file = Path(file)
    suffix = file.suffix.lower()
    if suffix == '.svg':
        points = read_svg(file, tolerance)
    elif suffix == '.dxf':
        points = read_dxf(file, tolerance)
    else:
        raise ValueError('outline files must be .svg or .dxf, not {!r}'.format(file.name))
    points = np.asarray(points, dtype=float)*float(scale)
    if tolerance is not None:
        points = simplify_polygon(points, tolerance*float(scale) if scale else tolerance)
    if center is not None:
        shift = np.asarray(center, dtype=float) - (points.min(axis=0) + points.max(axis=0))/2
        if np.linalg.norm(shift) > 1e-12*max(np.ptp(points, axis=0).max(), 1e-300): # an outline drawn in place stays exact
            points = points + shift
    return polygon_vertices(points)


# =============================================================================
# Douglas-Peucker simplification
# =============================================================================

def simplify_polygon(points, tolerance):
    """
    A closed polygon with fewer vertices that stays within tolerance of the given one
    (Douglas-Peucker): the polygon is split at its two vertices farthest apart, the two chains are
    simplified keeping their ends, and joined again. Repeated vertices are dropped first.
    """
    p = np.asarray(points, dtype=float).reshape(-1, 2)
    if len(p) > 1 and np.linalg.norm(p[0] - p[-1]) <= 1e-12*max(np.ptp(p, axis=0).max(), 1e-300):
        p = p[:-1]
    keep = np.ones(len(p), dtype=bool)
    keep[1:] = np.linalg.norm(np.diff(p, axis=0), axis=1) > 0
    p = p[keep]
    if len(p) <= 3 or tolerance <= 0:
        return p
    i0 = 0
    i1 = int(np.argmax(np.linalg.norm(p - p[0], axis=1)))
    chain1 = _douglas_peucker(p[i0:i1 + 1], tolerance)
    chain2 = _douglas_peucker(np.vstack((p[i1:], p[0:1])), tolerance)
    return np.vstack((chain1[:-1], chain2[:-1]))


def _douglas_peucker(chain, tolerance):
    """The open chain of points with its ends kept and the in-between points dropped where the chain
    stays within tolerance of the straight line (iterative, so that long chains do not recurse deeply)"""
    chain = np.asarray(chain, dtype=float)
    n = len(chain)
    keep = np.zeros(n, dtype=bool)
    keep[0] = keep[-1] = True
    stack = [(0, n - 1)]
    while stack:
        a, b = stack.pop()
        if b - a < 2:
            continue
        seg = chain[b] - chain[a]
        length = np.linalg.norm(seg)
        rel = chain[a + 1:b] - chain[a]
        if length > 0:
            distance = np.abs(cross2D(np.broadcast_to(seg, rel.shape), rel))/length
        else:
            distance = np.linalg.norm(rel, axis=1)
        k = int(np.argmax(distance))
        if distance[k] > tolerance:
            keep[a + 1 + k] = True
            stack.append((a, a + 1 + k))
            stack.append((a + 1 + k, b))
    return chain[keep]


# =============================================================================
# SVG
# =============================================================================

_NUMBER = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')
_PATH_TOKEN = re.compile(r'[MmLlHhVvCcSsQqTtAaZz]|[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')


def read_svg(file, tolerance=None):
    """The outline of an SVG file (see read_outline), in mm when the file gives its size in units, else in
    its user units (px); the y axis is flipped so that the drawing is not mirrored"""
    tree = ET.parse(os.fspath(file))
    root = tree.getroot()
    unit_scale = _svg_unit_scale(root)
    best, best_area = None, 0.0
    for element, transform in _svg_elements(root, np.eye(3)):
        tag = element.tag.split('}')[-1]
        loops = _svg_shape_loops(tag, element, tolerance)
        for loop in loops:
            if len(loop) < 3:
                continue
            pts = np.asarray(loop, dtype=float)
            pts = (transform @ np.vstack((pts.T, np.ones(len(pts)))))[0:2].T
            area = abs(0.5*np.sum(cross2D(pts, np.roll(pts, -1, axis=0))))
            if area > best_area:
                best, best_area = pts, area
    if best is None:
        raise ValueError('{}: no closed path, polygon, rect, circle or ellipse found'.format(file))
    best = best*unit_scale
    best[:, 1] = -best[:, 1]  # SVG's y axis points down
    return best


def _svg_unit_scale(root):
    """mm per user unit: from width (with a unit) and viewBox, else px at 96 dpi"""
    width = root.get('width', '')
    viewbox = root.get('viewBox') or root.get('viewbox')
    m = re.match(r'\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)\s*([a-z%]*)\s*$', width or '')
    if m and viewbox and m.group(2) in SVG_UNITS_MM and m.group(2) != '%':
        parts = _NUMBER.findall(viewbox)
        unit = SVG_UNITS_MM[m.group(2)]
        if len(parts) == 4 and float(parts[2]) > 0 and unit is not None:
            return float(m.group(1))*unit/float(parts[2])
    return SVG_UNITS_MM['px']


def _svg_elements(element, transform):
    """Every element with the transform of itself and its ancestors (a 3x3 matrix), depth first"""
    t = transform @ _svg_transform(element.get('transform', ''))
    yield element, t
    for child in element:
        yield from _svg_elements(child, t)


def _svg_transform(text):
    """3x3 matrix of an SVG transform attribute (matrix, translate, scale, rotate, skewX, skewY)"""
    m = np.eye(3)
    for name, args in re.findall(r'(matrix|translate|scale|rotate|skewX|skewY)\s*\(([^)]*)\)', text or ''):
        v = [float(x) for x in _NUMBER.findall(args)]
        t = np.eye(3)
        if name == 'matrix' and len(v) == 6:
            t[0, 0], t[1, 0], t[0, 1], t[1, 1], t[0, 2], t[1, 2] = v
        elif name == 'translate':
            t[0, 2] = v[0] if v else 0.0
            t[1, 2] = v[1] if len(v) > 1 else 0.0
        elif name == 'scale':
            t[0, 0] = v[0] if v else 1.0
            t[1, 1] = v[1] if len(v) > 1 else t[0, 0]
        elif name == 'rotate':
            a = math.radians(v[0]) if v else 0.0
            r = np.array([[math.cos(a), -math.sin(a), 0], [math.sin(a), math.cos(a), 0], [0, 0, 1]])
            if len(v) >= 3:
                c = np.array([[1, 0, v[1]], [0, 1, v[2]], [0, 0, 1]])
                ci = np.array([[1, 0, -v[1]], [0, 1, -v[2]], [0, 0, 1]])
                t = c @ r @ ci
            else:
                t = r
        elif name == 'skewX':
            t[0, 1] = math.tan(math.radians(v[0])) if v else 0.0
        elif name == 'skewY':
            t[1, 0] = math.tan(math.radians(v[0])) if v else 0.0
        m = m @ t
    return m


def _svg_shape_loops(tag, element, tolerance):
    """The closed loops (lists of [x, y]) drawn by one SVG element"""
    n = _curve_segments(tolerance)
    if tag == 'path':
        return [loop for loop in svg_path_subpaths(element.get('d', ''), n) if len(loop) >= 3]
    if tag in ('polygon', 'polyline'):
        v = [float(x) for x in _NUMBER.findall(element.get('points', ''))]
        return [[[v[i], v[i + 1]] for i in range(0, len(v) - 1, 2)]]
    if tag == 'rect':
        x, y = float(element.get('x', 0)), float(element.get('y', 0))
        w, h = float(element.get('width', 0)), float(element.get('height', 0))
        return [[[x, y], [x + w, y], [x + w, y + h], [x, y + h]]]
    if tag in ('circle', 'ellipse'):
        cx, cy = float(element.get('cx', 0)), float(element.get('cy', 0))
        rx = float(element.get('r', element.get('rx', 0)))
        ry = float(element.get('r', element.get('ry', 0)))
        t = np.linspace(0, 2*np.pi, max(4*n, 16), endpoint=False)
        return [np.column_stack((cx + rx*np.cos(t), cy + ry*np.sin(t))).tolist()]
    return []


def _curve_segments(tolerance):
    return CURVE_SEGMENTS if tolerance is None else max(CURVE_SEGMENTS, 64)


def svg_path_subpaths(d, segments=CURVE_SEGMENTS):
    """
    The subpaths of an SVG path (its d attribute) as lists of [x, y] points, curves (C, S, Q, T, A)
    flattened into straight segments. Every subpath is returned closed (a Z or not).
    """
    tokens = _PATH_TOKEN.findall(d)
    subpaths = []
    current = []
    pos = np.zeros(2)
    start = np.zeros(2)
    last_ctrl = None
    last_cmd = ''
    i = 0
    cmd = None

    def take(k):
        nonlocal i
        values = [float(x) for x in tokens[i:i + k]]
        if len(values) < k:
            raise ValueError('malformed SVG path near {!r}'.format(' '.join(tokens[max(i - 3, 0):i + k])))
        i += k
        return values

    def flush():
        if len(current) >= 3:
            subpaths.append([list(p) for p in current])
        current.clear()

    while i < len(tokens):
        if re.match(r'[A-Za-z]', tokens[i]):
            cmd = tokens[i]
            i += 1
            if cmd in 'Zz':
                flush()
                pos = start.copy()
                last_cmd = cmd
                continue
        elif cmd is None:
            raise ValueError('SVG path must start with a command')
        elif cmd in 'Mm':  # further pairs after a moveto are linetos
            cmd = 'L' if cmd == 'M' else 'l'
        rel = cmd.islower()
        c = cmd.upper()
        if c == 'M':
            x, y = take(2)
            flush()
            pos = pos + [x, y] if rel else np.array([x, y])
            start = pos.copy()
            current.append(pos.copy())
            last_ctrl = None
        elif c == 'L':
            x, y = take(2)
            pos = pos + [x, y] if rel else np.array([x, y])
            current.append(pos.copy())
            last_ctrl = None
        elif c == 'H':
            x, = take(1)
            pos = np.array([pos[0] + x if rel else x, pos[1]])
            current.append(pos.copy())
            last_ctrl = None
        elif c == 'V':
            y, = take(1)
            pos = np.array([pos[0], pos[1] + y if rel else y])
            current.append(pos.copy())
            last_ctrl = None
        elif c in 'CS':
            if c == 'C':
                x1, y1, x2, y2, x, y = take(6)
                p1 = pos + [x1, y1] if rel else np.array([x1, y1])
            else:
                x2, y2, x, y = take(4)
                p1 = 2*pos - last_ctrl if (last_ctrl is not None and last_cmd.upper() in 'CS') else pos.copy()
            p2 = pos + [x2, y2] if rel else np.array([x2, y2])
            p3 = pos + [x, y] if rel else np.array([x, y])
            t = np.linspace(0, 1, segments + 1)[1:, np.newaxis]
            pts = (1 - t)**3*pos + 3*(1 - t)**2*t*p1 + 3*(1 - t)*t**2*p2 + t**3*p3
            current.extend(pts.tolist())
            last_ctrl, pos = p2, p3
        elif c in 'QT':
            if c == 'Q':
                x1, y1, x, y = take(4)
                p1 = pos + [x1, y1] if rel else np.array([x1, y1])
            else:
                x, y = take(2)
                p1 = 2*pos - last_ctrl if (last_ctrl is not None and last_cmd.upper() in 'QT') else pos.copy()
            p2 = pos + [x, y] if rel else np.array([x, y])
            t = np.linspace(0, 1, segments + 1)[1:, np.newaxis]
            pts = (1 - t)**2*pos + 2*(1 - t)*t*p1 + t**2*p2
            current.extend(pts.tolist())
            last_ctrl, pos = p1, p2
        elif c == 'A':
            rx, ry, rot, large, sweep, x, y = take(7)
            end = pos + [x, y] if rel else np.array([x, y])
            current.extend(_arc_points(pos, end, rx, ry, rot, large != 0, sweep != 0, segments).tolist())
            pos = end
            last_ctrl = None
        else:
            raise ValueError('unsupported SVG path command {!r}'.format(cmd))
        last_cmd = cmd
    flush()
    return subpaths


def _arc_points(p0, p1, rx, ry, rotation, large, sweep, segments):
    """Points along an SVG elliptical arc from p0 to p1 (the end point included, the start not)"""
    p0, p1 = np.asarray(p0, dtype=float), np.asarray(p1, dtype=float)
    rx, ry = abs(rx), abs(ry)
    if rx == 0 or ry == 0 or np.allclose(p0, p1):
        return p1[np.newaxis]
    phi = math.radians(rotation)
    cos, sin = math.cos(phi), math.sin(phi)
    d = (p0 - p1)/2
    x1 = cos*d[0] + sin*d[1]
    y1 = -sin*d[0] + cos*d[1]
    lam = (x1/rx)**2 + (y1/ry)**2
    if lam > 1:
        rx, ry = rx*math.sqrt(lam), ry*math.sqrt(lam)
    num = rx*rx*ry*ry - rx*rx*y1*y1 - ry*ry*x1*x1
    den = rx*rx*y1*y1 + ry*ry*x1*x1
    coef = math.sqrt(max(num/den, 0.0)) if den > 0 else 0.0
    if large == sweep:
        coef = -coef
    cx1, cy1 = coef*rx*y1/ry, -coef*ry*x1/rx
    center = np.array([cos*cx1 - sin*cy1, sin*cx1 + cos*cy1]) + (p0 + p1)/2

    def angle(u, v):
        a = math.atan2(u[0]*v[1] - u[1]*v[0], u[0]*v[0] + u[1]*v[1])
        return a
    v0 = np.array([(x1 - cx1)/rx, (y1 - cy1)/ry])
    v1 = np.array([(-x1 - cx1)/rx, (-y1 - cy1)/ry])
    theta0 = angle(np.array([1.0, 0.0]), v0)
    dtheta = angle(v0, v1)
    if not sweep and dtheta > 0:
        dtheta -= 2*math.pi
    elif sweep and dtheta < 0:
        dtheta += 2*math.pi
    t = theta0 + dtheta*np.linspace(0, 1, segments + 1)[1:]
    pts = np.column_stack((rx*np.cos(t), ry*np.sin(t)))
    pts = pts @ np.array([[cos, sin], [-sin, cos]]) + center
    pts[-1] = p1
    return pts


# =============================================================================
# DXF
# =============================================================================

def read_dxf(file, tolerance=None):
    """The outline of a DXF file (see read_outline), in mm when the header gives $INSUNITS, else as drawn"""
    with open(os.fspath(file), 'r', errors='replace') as f:
        lines = [line.rstrip('\r\n') for line in f]
    pairs = [(lines[k].strip(), lines[k + 1].strip()) for k in range(0, len(lines) - 1, 2)]
    unit_scale = 1.0
    entities = []   # (type, {code: [values]})
    section = None
    current = None
    for code, value in pairs:
        if code == '0':
            if value == 'SECTION':
                section = None
            elif value == 'ENDSEC':
                section = None
            elif section == 'ENTITIES':
                current = {}
                entities.append((value, current))
            else:
                current = None
            continue
        if code == '2' and section is None and value in ('HEADER', 'ENTITIES', 'BLOCKS', 'TABLES', 'CLASSES', 'OBJECTS'):
            section = value
            continue
        if section == 'HEADER' and code == '9':
            current = {'name': value}
            continue
        if section == 'HEADER' and current is not None and current.get('name') == '$INSUNITS' and code == '70':
            unit_scale = DXF_UNITS_MM.get(int(float(value)), None) or 1.0
            continue
        if section == 'ENTITIES' and current is not None:
            current.setdefault(code, []).append(value)
    n = _curve_segments(tolerance)
    loops = []
    open_pieces = []  # chains of points from LINE and ARC entities, to be joined into loops
    vertex_run = None
    for kind, data in entities:
        if kind == 'LWPOLYLINE':
            xs, ys = [float(v) for v in data.get('10', [])], [float(v) for v in data.get('20', [])]
            bulges = [float(v) for v in data.get('42', [])]
            closed = int(float(data.get('70', ['0'])[0])) & 1
            loops.append(_polyline_points(list(zip(xs, ys)), bulges, closed, n))
        elif kind == 'POLYLINE':
            vertex_run = {'points': [], 'bulges': [], 'closed': int(float(data.get('70', ['0'])[0])) & 1}
        elif kind == 'VERTEX' and vertex_run is not None:
            vertex_run['points'].append((float(data.get('10', ['0'])[0]), float(data.get('20', ['0'])[0])))
            vertex_run['bulges'].append(float(data.get('42', ['0'])[0]))
        elif kind == 'SEQEND' and vertex_run is not None:
            loops.append(_polyline_points(vertex_run['points'], vertex_run['bulges'], vertex_run['closed'], n))
            vertex_run = None
        elif kind == 'LINE':
            p0 = (float(data['10'][0]), float(data['20'][0]))
            p1 = (float(data['11'][0]), float(data['21'][0]))
            open_pieces.append(np.array([p0, p1]))
        elif kind == 'ARC':
            c = np.array([float(data['10'][0]), float(data['20'][0])])
            r = float(data['40'][0])
            a0, a1 = math.radians(float(data['50'][0])), math.radians(float(data['51'][0]))
            if a1 <= a0:
                a1 += 2*math.pi
            t = np.linspace(a0, a1, max(2, int(math.ceil(n*(a1 - a0)/(2*math.pi))) + 1))
            open_pieces.append(np.column_stack((c[0] + r*np.cos(t), c[1] + r*np.sin(t))))
        elif kind == 'CIRCLE':
            c = np.array([float(data['10'][0]), float(data['20'][0])])
            r = float(data['40'][0])
            t = np.linspace(0, 2*np.pi, max(4*n, 16), endpoint=False)
            loops.append(np.column_stack((c[0] + r*np.cos(t), c[1] + r*np.sin(t))))
    loops.extend(_join_pieces(open_pieces))
    best, best_area = None, 0.0
    for loop in loops:
        pts = np.asarray(loop, dtype=float)
        if len(pts) < 3:
            continue
        area = abs(0.5*np.sum(cross2D(pts, np.roll(pts, -1, axis=0))))
        if area > best_area:
            best, best_area = pts, area
    if best is None:
        raise ValueError('{}: no closed polyline, loop of lines and arcs, or circle found'.format(file))
    return best*unit_scale


def _polyline_points(points, bulges, closed, n):
    """The vertices of a polyline with its bulge arcs flattened (closed or not: the loop is taken as closed)"""
    pts = [np.asarray(p, dtype=float) for p in points]
    out = []
    m = len(pts)
    for k in range(m):
        out.append(pts[k])
        bulge = bulges[k] if k < len(bulges) else 0.0
        nxt = pts[(k + 1) % m]
        if bulge != 0 and (closed or k < m - 1) and np.linalg.norm(nxt - pts[k]) > 0:
            theta = 4*math.atan(bulge)  # included angle of the arc
            chord = nxt - pts[k]
            d = np.linalg.norm(chord)
            r = d/(2*math.sin(abs(theta)/2))
            mid = (pts[k] + nxt)/2
            h = r*math.cos(theta/2)
            normal = np.array([-chord[1], chord[0]])/d
            center = mid + (-normal*h if theta > 0 else normal*h)
            # points along the arc between the two vertices
            a0 = math.atan2(pts[k][1] - center[1], pts[k][0] - center[0])
            steps = max(2, int(math.ceil(n*abs(theta)/(2*math.pi))))
            for s in range(1, steps):
                a = a0 + theta*s/steps
                out.append(center + r*np.array([math.cos(a), math.sin(a)]))
    return np.array(out)


def _join_pieces(pieces, tol=None):
    """Closed loops made of open chains whose ends meet"""
    pieces = [np.asarray(p, dtype=float) for p in pieces if len(p) >= 2]
    if not pieces:
        return []
    if tol is None:
        extent = max(np.ptp(np.vstack(pieces), axis=0).max(), 1e-300)
        tol = 1e-6*extent
    loops = []
    used = [False]*len(pieces)
    for s in range(len(pieces)):
        if used[s]:
            continue
        used[s] = True
        chain = pieces[s].copy()
        while True:
            if np.linalg.norm(chain[-1] - chain[0]) <= tol and len(chain) > 2:
                loops.append(chain[:-1])
                break
            found = False
            for k in range(len(pieces)):
                if used[k]:
                    continue
                piece = pieces[k]
                if np.linalg.norm(piece[0] - chain[-1]) <= tol:
                    chain = np.vstack((chain, piece[1:]))
                elif np.linalg.norm(piece[-1] - chain[-1]) <= tol:
                    chain = np.vstack((chain, piece[::-1][1:]))
                else:
                    continue
                used[k] = True
                found = True
                break
            if not found:
                break
    return loops


def outline_toml(points, decimals=6):
    """The vertices as a TOML line: box_polygon = [[x, y], ...]"""
    return 'box_polygon = [' + ', '.join('[{}, {}]'.format(round(float(x), decimals), round(float(y), decimals))
                                         for x, y in np.asarray(points)) + ']'
