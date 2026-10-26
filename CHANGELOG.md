# Changelog

## v1.3.0

Clip the cells with any polygon, also a non-convex one such as a dog-bone tensile specimen.

### New
- `box_shape = "polygon"` with `box_polygon = [[x, y], ...]`: clip the cells with any simple polygon, convex or not (for example a dog-bone tensile specimen: `examples/dogbone.toml`). Cell walls that leave and re-enter a non-convex polygon are cut into their parts inside it. With `merge_operation`, the short walls are merged by a polygon-aware version of the merge: closest pairs first, fused vertices not rounded to the `merge_tol` grid, and a merge is skipped when it would move a wall out of the polygon (across a notch or neck narrower than `merge_tol`) or move a corner. The bulk volume and porosity use the area of the polygon, also with `boundaryFlag = false` (the other shapes use the convex hull of the boundary points, which leaves out the box corners when `boundaryFlag` is off). `Clipping_Box` takes the vertices as `polygon=...`, `ModelInfo` as `polygon=...`; `polygon_vertices()` checks them (at least 3 vertices, non-zero area, no crossing edges). The square, triangle and hexagon boxes give the same meshes as before.

## v1.2.0

A real entry point instead of copying a template script, and much faster mesh generation.

### New
- `ringspy run wood_cube.toml`: command line tool (also `python -m ringspy run ...`) that generates a mesh from a TOML input file. `--seed`, `--path` and `--geoName` override the values in the file.
- `ringspy.generate(config, **changes)` runs the full generation of the template scripts (placement, tessellation and clipping, extrusion, all output files) and returns a `Result` with the output folder, the cell sites and the model properties. `config` is a `ringspy.Config`, a TOML file name or a dict.
- `ringspy.Config`: dataclass of all input parameters, with the names of the template scripts. Values are checked when it is made (types, box shapes, growth rules, unknown or misspelled names), so mistakes show up before a long run. Parameters left out take the values of the `wood_cube` example; the ring widths and `z_max` follow `r_max`, `nrings` and `box_size` as in the templates.
- `seed` parameter: the same seed gives the same files. It is drawn through `numpy.random.default_rng`; without a seed the cells come from NumPy's global random state, as before. The report log has a new `seed:` line.
- `CellPlacement_Binary`, `CellPlacement_Binary_Lloyd` and `CellPlacement_Honeycomb` take an optional `rng` argument (a seed or a `numpy.random.Generator`; default: the global random state, the same numbers as before).
- Flags (`boundaryFlag`, `merge_operation`, `precrackFlag`, `stlFlag`, `inpFlag`) accept booleans. The `'on'`/`'off'` strings still work.
- New parameters: `precrack_nodes` (line segments of the precracks; the template scripts' crack is the default), `save_sites` (also save `geoName_sites.npy`/`geoName_radii.npy`, to mesh the same cells again with a `.npy` growth rule) and `sites_dir` (where to look for the `.npy` files; default: the folder of the TOML file).
- Examples in `examples/`: `wood_cube`, `wood_cube_lloyd`, `honeycomb_hex`, `honeycomb_tri` and `saved_sites`, each a TOML file and a short Python script. They replace the template scripts `tests/test_*.py`.
- `VisualizationFiles(..., vtk_format='binary')` writes the ParaView files with raw binary appended data: faster to write and read, and about 3 times smaller. The default is still `'ascii'`.

### Performance
Mesh generation is much faster, with identical output files for the same random seed (checked byte for byte against v1.1.0 for all template scripts, with and without merging, boundary and precracks). The `wood_cube` template takes 2.5 s instead of 28 s with `iter_max = 500` (3.3k cells), and 6 s instead of 106 s with `iter_max = 2000` (10k cells), and needs less memory; `benchmarks/benchmark_pipeline.py` reports the time per stage.
- STL export: the surface quads are split into triangles with `vtkTriangleFilter` before `vtkSTLWriter`, which splits them about 100 times slower itself (same STL file). Only the needed VTK modules are imported.
- `CellPlacement_Binary` checks new cells for overlaps only against the placed cells in the neighbouring cells of a grid, instead of all placed cells (same sites).
- `RebuildVoronoi_merge` finds points by their coordinates with dictionaries instead of comparing with all points, and replaces merged vertex indices in one pass.
- `VisualizationFiles`, `AbaqusFile` and the `.mesh`/`-mesh.txt` files format their data in large chunks instead of line by line, and `VisualizationFiles`, `GenerateBeamElement`, `ConnectorMeshFile` and `VertexandRidgeinfo` build their arrays without Python loops over vertices, ridges or connectors.
- The preview plot draws the Voronoi ridges as one `LineCollection` instead of one line per ridge (same image).

### Changed
- Python 3.9 is no longer supported (it is end-of-life): RingsPy requires Python >= 3.10.
- Template scripts of earlier versions still run unchanged and give the same files as `ringspy.generate` for the same random state (checked by `tests/test_legacy_template.py`).
- An unsupported `inpType` is an error when the configuration is made; the template scripts first saved the cells, then failed after the whole generation.
- New dependency on Python < 3.11: `tomli` (Python 3.11+ reads TOML with the standard library).
- Tests: the end-to-end tests run every example and check the files, that the cell vertices lie in the clipping box, and that a second run with the same seed (through the command line) gives identical files. They are marked `slow` (`pytest -m "not slow"` skips them).

### Fixes
- Wall thickness: `cellwallthickness_sparse` was applied to the heart and dense rings and `cellwallthickness_dense` to the sparse rings, opposite to the names and the README (present since v1.0). The parameters now act on the rings they name. To keep the meshes of the examples unchanged, the defaults and the examples swapped their values (sparse 0.006, dense 0.010); meshes you made with your own values need the two values swapped to stay the same. The Abaqus section sizes that used `cellwallthickness_sparse` now use `cellwallthickness_dense`, so they keep their old values with the swapped defaults.
- `merge_operation = 'on'` crashed in `GenerateBeamElement` (`IndexError`, or `ValueError` from `np.repeat`) in about 5% of the random `wood_cube` geometries. A vertex near a box corner whose neighbouring vertices are all outside the box has only cut ridges, and `RebuildVoronoi_merge` left it out of `nvertices_in`, so `npt_per_layer` was one less than the number of vertices. `RidgeMidQuarterPts` now raises a clear error if these counts disagree. Meshes that were generated before are unchanged.
- `CellPlacement_Binary_Lloyd` crashed (`ValueError: need at least one array to concatenate`) when the random ring widths made a ring narrower than one cell. Such a ring now only gets the cells on its outer edge.
- `CellPlacement_Binary_Lloyd` crashed (`QhullError`) when the heart region was narrower than about 2/3 of a cell. A heart narrower than one cell has no cells inside, so the Lloyd relaxation is now skipped for it; the innermost cells of the first ring cover the centre, as they already did for hearts between 2/3 and 1 cell wide.

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
