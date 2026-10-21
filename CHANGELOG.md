# Changelog

## Unreleased

### Performance
Mesh generation is much faster, with identical output files for the same random seed (checked byte for byte against v1.1.0 for all template scripts, with and without merging, boundary and precracks). The `wood_cube` template takes 2.5 s instead of 28 s with `iter_max = 500` (3.3k cells), and 6 s instead of 106 s with `iter_max = 2000` (10k cells), and needs less memory; `benchmarks/benchmark_pipeline.py` reports the time per stage.
- STL export: the surface quads are split into triangles with `vtkTriangleFilter` before `vtkSTLWriter`, which splits them about 100 times slower itself (same STL file). Only the needed VTK modules are imported.
- `CellPlacement_Binary` checks new cells for overlaps only against the placed cells in the neighbouring cells of a grid, instead of all placed cells (same sites).
- `RebuildVoronoi_merge` finds points by their coordinates with dictionaries instead of comparing with all points, and replaces merged vertex indices in one pass.
- `VisualizationFiles`, `AbaqusFile` and the `.mesh`/`-mesh.txt` files format their data in large chunks instead of line by line, and `VisualizationFiles`, `GenerateBeamElement`, `ConnectorMeshFile` and `VertexandRidgeinfo` build their arrays without Python loops over vertices, ridges or connectors.
- The preview plot draws the Voronoi ridges as one `LineCollection` instead of one line per ridge (same image).

### New
- `VisualizationFiles(..., vtk_format='binary')` writes the ParaView files with raw binary appended data: faster to write and read, and about 3 times smaller. The default is still `'ascii'`.

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
