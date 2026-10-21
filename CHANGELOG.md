# Changelog

## v1.1.0

Maintenance release: RingsPy installs and runs again on current Python, NumPy and setuptools.

### Compatibility
- Removed the `pkg_resources` import, which made `import ringspy.MeshGenTools` fail with setuptools >= 81.
- Replaced `np.cross` on 2D vectors (an error in NumPy 2) with `MeshGenTools.cross2D`.
- Requires Python >= 3.9. Packaging moved to `pyproject.toml` (`setup.py`, `setup.cfg`, `tox.ini` and the requirements files are gone).
- `vtk` and `hexalattice` are optional: `pip install "ringspy[all]"`, or the `[vtk]` / `[hex]` extras. A clear error explains what to install if they are missing. `hexalattice` is no longer pinned to 1.2.1.
- CI runs on current GitHub runners (Python 3.9-3.14 on Ubuntu, 3.13 on Windows and macOS). PyPI publishing supports Trusted Publishing.

### Fixes
- Hexagon clipping box: `box_center` x and y were swapped, and `RebuildVoronoi` replaced `box_center` with `generation_center`.
- Triangle clipping box with `merge_operation = 'off'` could hang forever while connecting the boundary points.
- `merge_operation = 'on'` clipped triangle and hexagon boxes with their bounding rectangle and crashed with random cell sites; merged points rounded to `merge_tol` could also end up slightly outside slanted box edges.
- `RebuildVoronoi_merge` crashed (`IndexError`) when the clipping box extended beyond the cell sites.
- Boundary points are now connected in order along the box perimeter (the nearest-neighbour search could jump across the box).
- `CellPlacement_Binary_Lloyd` ignored `generation_center`.
- `ReadSavedSites` looped for 180 s on Windows when the file was missing; it now raises `FileNotFoundError` right away.
- Precracks: the Abaqus input file referenced a missing `VisualLongConns` element set, only the first segment was cracked, longitudinal connectors were "cracked" by unrelated indices, `AllElles` counted precracked connectors twice, and cracks longer than a few cells were only cut near their midpoint.
- The `path` input parameter is now used for all output files (they always went to `meshes/`).
- Library functions raise exceptions instead of calling `exit()`.
- Template scripts call `ReadSavedSites` with the right arguments.

### Other changes
- STL files are written in binary format (about 3 times smaller).
- Removed unused material property arrays from `AbaqusFile`.
- New regression tests (`tests/test_regressions.py`); tests write their output to temporary folders.
