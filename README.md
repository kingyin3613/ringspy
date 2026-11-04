# RingsPy
[![Python 3](https://img.shields.io/static/v1?label=Python&logo=Python&color=3776AB&message=3)](https://www.python.org/)
[![PyPI version](https://badge.fury.io/py/ringspy.svg)](https://badge.fury.io/py/ringspy)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/ringspy/badges/version.svg)](https://anaconda.org/conda-forge/ringspy)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
[![Tests](https://github.com/kingyin3613/ringspy/actions/workflows/tests.yml/badge.svg)](https://github.com/kingyin3613/ringspy/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/kingyin3613/ringspy/branch/main/graph/badge.svg?token=4AZN3HGGET)](https://codecov.io/gh/kingyin3613/ringspy)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04945/status.svg)](https://doi.org/10.21105/joss.04945)

![RingsPy: cell seeds placed ring by ring, Voronoi tessellation, clipping to any polygon (a dog-bone, a star, an irregular polygon, a hexagon) and extrusion with a twist into a 3D cellular mesh, exported for ParaView, 3D printing and Abaqus](./contents/ringspy_hero.gif)

RingsPy is a Voronoi diagrams-based geometric generation tool that generates 3D meshes and models of prismatic cellular solids with radial growth rules.

## Dependencies and Installation
RingsPy needs Python 3.10 or newer and depends on mainstream Python libraries ``numpy``, ``scipy`` and ``matplotlib``. It optionally depends on library ``hexalattice``, if a regular hexagonal lattice (e.g. honeycomb) is wanted, and on ``vtk``, if the 3D model (STL) files are also wanted.
### 1. pip install
To install RingsPy with all optional features, one may use `pip`:
```bash
pip install "ringspy[all]"
``` 
`pip install ringspy` installs the core package only; add the extras `[vtk]` (STL export) or `[hex]` (`regular_hexagonal` growth rule) as needed. One may also use:
```bash
pip install "ringspy[all] @ git+https://github.com/kingyin3613/ringspy.git"
``` 
to get updates beyond the latest release. 

### 2. conda install
To install RingsPy in `conda` environments, one can also use `conda-forge` channel:
```bash
conda install -c conda-forge ringspy
``` 

### 3. Installation Check

There are some unit tests in [tests](https://github.com/kingyin3613/ringspy/tree/main/tests/). One can use ``pytest`` to check whether the installation is successful by installing the test dependencies in the cloned ``ringspy`` directory and then running ``pytest`` there (generated files go to temporary folders):

```bash
pip install -e ".[test]"
pytest
```

`pytest -m "not slow"` skips the end-to-end runs of the examples and takes well under a minute.

## Getting Started (Current version: v1.3.0)
Once all required components are installed and one is ready to begin, a path forward should be established for generating the mesh. The basic steps for running/viewing a cellular mesh are listed as the following:

    1. Edit geometry and algorithm parameters in an input file
    2. Generate mesh with the command `ringspy run` (or with `ringspy.generate` in Python)
    3. Visualize 2D view using Matplotlib or 3D model in `ParaView`
    4. (Optional) Export 3D STL model for 3D editing and/or printing
    5. (Optional) Export input file for numerical simulations with software `Abaqus`
    6. (Optional) Generate many meshes at once with `ringspy batch`, and export the cell mesh as a graph for machine learning

### 1. Geometry and Parameters
The first step to generate a cellular geometry is selecting geometry and appropriate parameters. 

### 1.1. Geometry
An input file in [TOML](https://toml.io) format, for example [`wood_cube.toml`](https://github.com/kingyin3613/ringspy/tree/main/examples/wood_cube.toml) in the [examples](https://github.com/kingyin3613/ringspy/tree/main/examples/) directory, holds the parameters for the generation of a cubic wood specimen. The other examples there are a wood specimen with more regular cells (`wood_cube_lloyd.toml`), a softwood with a gradual earlywood-latewood transition and rings that differ (`wood_gradual.toml`), honeycombs clipped by a hexagon or a triangle (`honeycomb_hex.toml`, `honeycomb_tri.toml`), a dog-bone tensile specimen clipped by a polygon (`dogbone.toml`), and a mesh of previously saved cells (`saved_sites.toml`). Copy the one closest to the wanted geometry and edit it.

*Note: The Mesh Generation Tool accepts the pre-defined boundary geometries triangle, square and hexagon, and any simple polygon given by its vertices (`box_shape = "polygon"`); importing of CAD and/or other 3D model files will be implemented in subsequent versions.*

*Note: for greatest compatibility create the geometry using all millimeters.*

### 1.2. Parameters

By opening an input file, e.g., `examples/wood_cube.toml` in any text editor, a file format similar to what is shown below will be displayed:
```toml
geoName = "wood_cube"
path = "meshes"                 # output folder: the files go to meshes/wood_cube/

# Cell placement
radial_growth_rule = "binary"   # binary, binary_lloyd, regular_hexagonal, or a saved *_sites.npy file
iter_max = 500                  # increase this number to achieve a more regular geometry
print_interval = 500            # interval for printing progress info
seed = 2023                     # random seed; remove it for a different geometry on every run

# Radial cell growth
r_min = 0                       # inner radius of generation domain
r_max = 4                       # outer radius of generation domain
nrings = 4                      # number of rings
width_heart = 0.3               # ring width of the innermost ring (left out: 0.3*(r_max-r_min)/nrings)
width_sparse = 0.7              # ring width of rings with sparse cells (left out: 0.7*(r_max-r_min)/nrings)
width_dense = 0.3               # ring width of rings with dense cells (left out: 0.3*(r_max-r_min)/nrings)
generation_center = [0, 0]      # coordinates of generation domain center

cellsize_sparse = 0.02          # a number, or one value per ring: [0.024, 0.020, 0.022, 0.018]
cellsize_dense = 0.01
cellwallthickness_sparse = 0.006
cellwallthickness_dense = 0.010
transition = "abrupt"           # abrupt, linear or sigmoid: change from sparse to dense within a ring
transition_width = 1.0          # fraction of the sparse part of a ring over which it changes (linear, sigmoid)

# Clipping box
boundaryFlag = true
box_shape = "square"            # square, triangle, hexagon or polygon (see box_polygon)
box_center = [0, 0]             # coordinates of clipping box center
box_size = 4.0                  # side length

# Longitudinal direction
nsegments = 2                   # increase this number to achieve a smoother transition in longitudinal direction if theta is large
theta_min = 0                   # unit: radian
theta_max = 0.05                # unit: radian
z_min = 0
z_max = 4.0                     # segment_length = (z_max - z_min) / nsegments (left out: box_size)
long_connector_ratio = 0.02     # longitudinal joint length = ratio * segment_length

# Material
skeleton_density = 1.5e-9       # unit: tonne/mm3

# Mesh operations
merge_operation = true          # merge cell walls shorter than merge_tol
merge_tol = 0.01

precrackFlag = false
precrack_widths = 0.1

# Output files
stlFlag = true                  # 3D model (needs vtk: pip install "ringspy[vtk]")
inpFlag = true                  # input file for simulations
inpType = "Abaqus"
save_sites = false              # also save wood_cube_sites.npy and wood_cube_radii.npy, to reuse the cells later
graphFlag = false               # also save wood_cube-graph.npz, the cell mesh as a graph (PyTorch Geometric layout)
```

A parameter that is left out takes its default value, which is its value in `wood_cube.toml`; a misspelled parameter name is an error. Flags take `true`/`false` (the `'on'`/`'off'` strings of earlier versions also work).

- `geoName` is the geometry name, `path` is the folder where the mesh files will be generated (relative to the folder where `ringspy run` is started).
- `radial_growth_rule` is the radial growth rule for cell placement. When a file name with extension`.npy` is specified, a saved cell data file will be loaded (for v1.x, choose one of these rules: `binary`, `binary_lloyd`, `regular_hexagonal`, or a file name with extension `.npy`). The `.npy` file and its `_radii.npy` partner are looked for in the folder of the input file and its parents, or in `sites_dir` if given; `save_sites = true` saves them for a generated geometry.
- `iter_max` is the max number of iteration for randomly placing non-overlapping cell particles in the 2D toroidal cell placement region. Noticing that, larger `iter_max` leads to more centroidal Voronoi cells, for more reference, see [Wiki - Centroidal Voronoi Tessellation](https://en.wikipedia.org/wiki/Centroidal_Voronoi_tessellation/).
- `print_interval` is the print interval when every n cell particles are placed in the placement region.
- `seed` makes the random cell placement repeatable: the same seed gives the same files. Without a seed, every run gives a different geometry (the cells then come from NumPy's global random state, as in earlier versions). The seed is written to the report log.
- `r_min` and `r_max` are the upper and lower bounds of radii of toroidal cell placement regions (generation rings), `nrings` is the total number of rings.
- `width_heart`, `width_sparse`, and `width_dense`, are ring widths for innermost ring, rings with sparse cells, and rings with dense cells, respectively, which all together determine the morphology of the cellular structure.
- `generation_center` is the location of the placement region.
- `cellsize_sparse`,`cellsize_dense`, `cellwallthickness_sparse`, and `cellwallthickness_dense` are parameters for the sparse and dense cells (the earlywood and latewood of a growth ring; the innermost ring, the heart, has dense cells). Each is a single number for all rings, or a list with one value per ring, innermost first (e.g. `cellwallthickness_dense = [0.009, 0.010, 0.011, 0.012]` with `nrings = 4`), so that every ring can differ; the heart takes the dense values of the first ring.
- `transition` is how the cell size and the wall thickness change from the sparse part of a ring to its dense part: `"abrupt"` (the default: a step at the boundary, as in all earlier versions), `"linear"` or `"sigmoid"` (an S-shaped curve). With the last two, the values change from the sparse to the dense values over the outer fraction `transition_width` of the sparse part (between 0 and 1, default 1: the whole sparse part), so that the two parts of a ring join without a jump, like the gradual earlywood-latewood transition of spruce or fir; the dense part, and the step at the ring boundary back to the next sparse part, are unchanged. The cells of the `binary` and `binary_lloyd` rules follow the profile when they are placed, and the wall thickness follows it at every vertex. See [examples/wood_gradual.toml](./examples/wood_gradual.toml).
- `boundaryFlag` flag can be turned on/off for generating neat boundaries consisting of grains.
- `box_shape` is the shape of cutting box: `square`, `triangle`, `hexagon`, or `polygon`.
- `box_center`, and `box_size` are for describing the cutting box (square, triangle and hexagon).
- `box_polygon` gives the vertices `[[x1, y1], [x2, y2], ...]` of the cutting box for `box_shape = "polygon"`: any simple polygon, convex or not (e.g. a dog-bone specimen, see [examples/dogbone.toml](./examples/dogbone.toml)), with its vertices in either direction around it. Without `z_max`, the specimen is as long as the polygon's largest width; pre-cracks need `precrack_nodes`. As for the triangle, the `.inp` file gets only the set of all nodes, no node sets for the sides of the box (those are written for `box_shape = "square"`, and for a hexagon whose `geoName` contains `hydrostatic_` or `uniaxial_`). Keep `boundaryFlag = true` for a non-convex polygon: without the walls along the boundary, a part narrower than about a cell (a sharp tip, say) can give loose pieces. With `merge_operation`, a merge that would move a wall out of the polygon (across a notch narrower than `merge_tol`, say) is skipped, and the corners never move, so polygon edges shorter than `merge_tol` stay as short walls: simplify a finely digitised outline first.
- `nsegments` is the number of segments consisting the prismatic cells during the z- (longitudinal) extrusion, with `segment_length = (z_max - z_min) / nsegments`.
- `theta_min` and `theta_max` determine the twisting angles (unit: radian) of the 2D mesh around the `generation_center`, during the spiral z-extrusion.
- `z_min` and `z_max` determine the boundaries of prismatic cells in z- (longitudinal) direction.
- `long_connector_ratio` is the length of longitudinal joints, with `longitudinal joint length = ratio * segment_length`.
- `skeleton_density` is the density of the skeleton (substance) material, e.g. density of wood cell walls in the wood microstructure. 
- `merge_operation` flag can be turned on/off for the merging operation, when flag is on, cell ridges that are shorter than the threshold `merge_tol` in the 2D mesh will be deleted, and their vertices will be merged respectively, the mesh will be reconstructed. This is designed for eliminating small cell ridges/walls which fall out of the resolution range of the 3D printing and for avoiding having elements with too small stable time increments in numerical simulations. 
- `precrackFlag` flag is for inserting pre-cracks, for the notched specimens. `precrack_nodes = [[x1, y1, x2, y2], ...]` gives the line segments of the cracks; when it is left out, a horizontal crack of length `0.1*box_size` through `box_center`, starting `0.12*box_size` from the left side of the box, is inserted. `precrack_widths` is the width of the cracks.
- `stlFlag` flag can be turned on/off for generating 3D STL files.
- `inpFlag` flag can be turned on/off for generating input files for numerical simulations.
- `inpType` indicates the software(s) that the mesh generation should generate input files for.
- `graphFlag` flag can be turned on/off for saving the 2D cell mesh as a graph (`geoName-graph.npz`, see section 5.2).

![growth_rule_binary](<./contents/growth_rule_binary.png>)
### 2.1. Run Mesh Generation
Open a Command Prompt or Terminal window, set the current directory to [examples](https://github.com/kingyin3613/ringspy/tree/main/examples/) or the folder of your own input file, then run the command:
```bash
ringspy run wood_cube.toml
```
(`python -m ringspy run wood_cube.toml` does the same.) The options `--seed`, `--path` and `--geoName` override the values in the file, e.g. `ringspy run wood_cube.toml --seed 7 --geoName wood_cube_7`.

Functions in `MeshGenTools` library will be called to create the initial mesh, wood cell particles following certain cell size distribution will be placed, then `Scipy.voronoi` function will be called to form the initial 2D Voronoi tessellation, additional code reforms the tesselation and generates the desired files. A successful generation will end with the line "Mesh files written to ..." in the Terminal window.

A new folder should have been created in the `.\meshes` directory with the same name as the `geoName` in `wood_cube.toml`.

The same generation can be run from Python, e.g. to loop over parameters:
```python
import ringspy

result = ringspy.generate("wood_cube.toml")                 # same as: ringspy run wood_cube.toml
for seed in range(5):
    ringspy.generate("wood_cube.toml", seed=seed, geoName=f"wood_cube_{seed}", stlFlag=False)

config = ringspy.Config(geoName="small_cube", r_max=2, box_size=1.0, seed=1)  # parameters left out: defaults
result = ringspy.generate(config)
print(result.output_dir, result.nParticles, result.porosity)
```
The `Result` also holds the clipped 2D cell mesh with its growth ring labels: `result.vertices` and `result.ridges` (the cell vertices, one beam each, and the cell walls as pairs of vertex indices), `result.vertex_ring`, `result.vertex_dense`, `result.vertex_boundary` and `result.vertex_thickness` (the growth ring of every vertex, 0 for the heart; whether it is in a dense part; whether it lies on the clipping boundary; its wall thickness), and `result.ridge_ring`, `result.ridge_dense`, `result.ridge_boundary` for the walls. Many meshes at once: section 5.1.

Each example also has a short Python script next to it, e.g. `python wood_cube.py`.

*Note: template scripts of earlier versions (e.g. `tests/test_wood_cube.py` of v1.1), which call the `MeshGenTools` functions one by one, still work. Their placement functions take an optional `rng` argument (a seed or a `numpy.random.Generator`) to make the cell placement repeatable.*

${\textbf{(New feature in version v0.4.x)}}$

The Lloyd's relaxation (a.k.a. Lloyd's algorithm) has been integrated with the mesh generation, one may use this approach by changing the following line in the input file (see [`wood_cube_lloyd.toml`](https://github.com/kingyin3613/ringspy/tree/main/examples/wood_cube_lloyd.toml)):
```toml
radial_growth_rule = "binary_lloyd"
``` 
This allows a rapid radial generation of more regular centroidal Voronoi cells from the `generation_center` (see the difference between `wood_binary` and `wood_binary_lloyd` in the preview of 2D cross-sections in section 2.2). For more detail of the algorithm, see [Wiki - Lloyd's algorithm](https://en.wikipedia.org/wiki/Lloyd's_algorithm)).

![lloyd](<./contents/lloyd.gif>)

### 2.2. Check Mesh and 3D Model Files

The following files should be generated in the `.\meshes\geoName` directory with a successful run:
- Log file
    - Log file for model generation: `wood_cube-report.log`
- Mesh files
    - Non-Uniform Rational B-Splines (NURBS) beam file: `wood_cubeIGA.txt`
    - connector data file: `wood_cube-mesh.txt`
    - Grain-ridge data file: `wood_cube-vertex.mesh`
    - Ridge data file: `wood_cube-ridge.mesh`
- Visualization files
    - Cross-sectional image file for initial 2D configuration: `wood_cube.png`
    - Paraview vtk file for initial vertex configuration: `wood_cube_vertices.vtu`
    - Paraview vtk file for initial grain solid configuration: `wood_cube_beams.vtu`
    - Paraview vtk file of initial cell ridge configuration: `wood_cube_conns.vtu`
    - Paraview vtk file of initial connector (volume) configuration: `wood_cube_conns_vol.vtu`
    - The `.vtu` files are written as text. For large models, `rpgen.VisualizationFiles(..., vtk_format='binary')` writes them as binary files instead, which are faster to write and read and about 3 times smaller.
- (Optional) 3D model files
    - Binary STL file of initial cellular solid configuration: `wood_cube.stl`
- (Optional) Abaqus input files
    - INP file of simulation input of initial cellular solid configuration in `Abaqus`: `wood_cube.inp`
- (Optional) Graph file
    - The 2D cell mesh as a graph, in the layout of PyTorch Geometric: `wood_cube-graph.npz` (`graphFlag = true`, see section 5.2)

The growth ring of every cell is written to the files: the `.vtu` files carry the cell data `Ring` (0 for the heart, then 1, 2, ...) and `Dense` (1 in the heart and in the dense part of a ring, 0 in the sparse part), the `.inp` file has an element set per ring (section 4.3), and the report log lists the radii of the rings. A beam and the wing quads around it take the ring of their vertex; a cell wall (and its transverse connectors) takes the ring of the vertex whose thickness it has, its second vertex; a longitudinal connector takes the ring of the wall of its wing, whose width it has. So the ring labels of the connectors always go with their `Connector_width`.


![Model2DPreview](<./contents/Model2DPreview.png>)
	
### 3. Visualization
A scientific visualization application `ParaView` can directly visualize the generated vtk files; It can also visualize generated 3D model STL files if the STL flag is on. `Paraview` is a free software, it can be downloaded from its official website: [https://www.paraview.org/download/](https://www.paraview.org/download/), latest version is recommeded.

### 3.1. Visualize Components of the 3D Model in ParaView
    1. Open ParaView
    2. Recommeded to temporarily turn off ray tracing
        - Uncheck "Enable Ray Tracing" (bottom left)
    3. Open File Sets in `.\meshes\geoName`
        - File > Open...
    4. Select the visualization files containing: `_vertices.vtu`, `_beams.vtu`, `_conns.vtu`
    5. Apply to visualize
        - Press Apply (left side, center)
    6. Turn on color plotting
        - Left click `_conns.vtu`, then select coloring (mid/lower left) > Connector_width
        - Or color by `Ring` (the growth ring of every beam and connector) or `Dense` (dense/sparse rings)
    7. Scale and position the image as desired
    8. Turn back on Ray Tracing
    9. Adjust Ray Tracing lighting and settings as desired
    10. Export Image
        - File > Save Screenshot
        - Enter a file name > OK
        - Leave new window as-is or increase resolution > OK

### 3.2. Visualize Volumes of the 3D Model in ParaView
    1. Open Paraview
    2. Recommeded to temporarily turn off ray tracing
        - Uncheck "Enable Ray Tracing" (bottom left)
    3. Open File Sets in `.\meshes\geoName`
        - File > Open...
    4. Select the newly created visualization files containing: `_conns_vol.vtu`
    5. Apply to visualize
        - Press Apply (left side, center)
    6. Turn on color plotting
        - Left click `_conns_vol.vtu`, then select coloring (mid/lower left) > Connector_width
    7. Scale and position the image as desired
    8. Turn back on Ray Tracing
    9. Adjust Ray Tracing lighting and settings as desired
    10. Export Image
        - File > Save Screenshot
        - Enter a file name > OK
        - Leave new window as-is or increase resolution > OK
![ModelVisualization](<./contents/ModelVisualization.png>)

### 4. (Optional) Numerical Simulation
The mesh generation tool can also prepare the input files for the numerical simulations of the cellular solid in other softwares. By now (version 1.x), the input file format, `.inp`, that is used in a finite element method (FEM) software `Abaqus` is supported, if the INP flag is on. `Abaqus` is a commerical software suite for integrated computer-aided engineering (CAE) and finite element analysis, own by `Dassault Systèmes`. One may refer to its [Wiki](https://en.wikipedia.org/wiki/Abaqus) for more about `Abaqus`, and to [Introduction](https://bertoldi.seas.harvard.edu/files/bertoldi/files/abaqusinputfilemanualv1.pdf?m=1444417191) for the introduction of Abaqus input files.

All steps for the model setup can be accomplished through manually coding the Abaqus input file in a text editor. The method used in the example procedure shown below requires access to the Abaqus GUI.

### 4.1 Create New Model
    1. Open Abaqus CAE
    2. Create a new model
        -	File > New Model Database > With Standard/Explicit Model

### 4.2 Import Meshed Part

    1. Import the meshed part
        - File > Import > Part
        - Select the newly created file which ends in `.inp` (note: you may need to change the File Filter to see this file)
        - Click OK
    2. Rename part
        - Under the model tree (left) expand `Parts`
        - Right click on the part > Rename...
        - Enter a new name less than 14 characters, no special symbols, and easily recognizable (i.e. "WoodCube")

*Note: if FEM parts will be added to the model, this RingsPy generated part must come first alphabetically. Also recommended not to include numbers in the name.*

![ModelNumericalSimulation](<./contents/ModelNumericalSimulation.png>)

### 4.3 One Material per Growth Ring
The `.inp` file defines, inside the instance, an element set for the beams, the transverse connectors and the longitudinal connectors of every growth ring: `Ring0Beams`, `Ring0TConns`, `Ring0LConns` for the heart, `Ring1Beams`, `Ring1TConns`, `Ring1LConns` for the first ring, and so on (only non-empty sets are written; `Ring5...` with `nrings = 4` holds the few elements beyond the outermost ring radius). `SparseBeams`, `SparseTConns`, `SparseLConns` and `DenseBeams`, `DenseTConns`, `DenseLConns` collect the elements of all sparse (earlywood) and all dense (latewood, and the heart) rings. In Abaqus CAE they appear as sets of the imported part, so a section with its own material can be assigned to every ring. To do the same in the input file, replace the `*Beam Section` of `AllBeams` and the `*Solid Section` lines of `AllBotConns`, `AllRegConns`, `AllTopConns` and `AllLongConns` by one section per ring set (an element must belong to one section only), for instance:
```
*Solid Section, elset=Ring1TConns, material=Earlywood
 2.5000E-03
*Solid Section, elset=Ring2TConns, material=Latewood
 2.5000E-03
```
The transverse connectors of a ring set have the wall thickness of that ring (the `width` column of `geoName-mesh.txt`), the longitudinal connectors the thickness of the wall of their wing.

### 5. (Optional) Batch Generation and Graph Export
Both are made for parameter studies and for building data sets, e.g. to train a graph neural network surrogate of the lattice model.

### 5.1 Many Meshes at Once
`ringspy batch` generates one mesh per combination of the swept parameter values and per random seed, and writes one table of their descriptors:
```bash
ringspy batch wood_cube.toml --seeds 0-9 --set cellsize_sparse=0.02,0.03 --set stlFlag=false --jobs 4
```
`--seeds` takes integers and ranges (`0-9` is 0 to 9), `--set NAME=V1,V2,...` the values of a parameter to sweep (TOML syntax: lists in brackets, e.g. `--set "cellwallthickness_dense=[0.01,0.01],[0.01,0.02]"`; a single value just changes the parameter in all runs), `--jobs` the number of runs in parallel, and `--path`/`--geoName` override the file. Run i is named `<geoName>_<i>` (`wood_cube_000`, `wood_cube_001`, ...) and gets the usual folder of files, with the first swept parameter varying slowest and the seeds fastest. The table `<path>/<geoName>-batch.csv` has one row per run with its seed, the swept parameters, and `nParticles`, `nvertex`, `nridge`, `nbeamElem`, `mass`, `bulk_volume`, `bulk_density` and `porosity`, plus the error message of a run that failed (the other runs are still generated; the command then exits with status 1). From Python:
```python
results = ringspy.batch("wood_cube.toml", seeds=range(10), sweep={"cellsize_sparse": [0.02, 0.03]},
                        jobs=4, stlFlag=False)        # a list of Results, None for a failed run
```

### 5.2 The Cell Mesh as a Graph
With `graphFlag = true` (or `--set graphFlag=true` in a batch), `geoName-graph.npz` holds the clipped 2D cell mesh as a graph in the layout of a [PyTorch Geometric](https://pytorch-geometric.readthedocs.io) `Data` object: the nodes are the cell vertices (the beams of the 3D model) and the edges the cell walls (the transverse connectors), each in both directions.
- `pos` (nvertex, 2): the vertex coordinates
- `x` (nvertex, 6): node features, named in `node_features`: growth ring, dense (1/0), boundary (1 on the clipping boundary), wall thickness at the vertex, degree (number of walls), distance to the generation center
- `edge_index` (2, 2*nridge) and `edge_attr` (2*nridge, 6): the edges and their features, named in `edge_features`: length, angle of the wall (radians, counter-clockwise from the x axis, from the first node of the edge to the second), thickness, growth ring, dense, boundary (1 for walls along the clipping boundary)
- `sites`, `radii`, `generation_center`: the cells and growth rings the mesh was made from; `nParticles`, `nbeamElem`, `mass`, `bulk_volume`, `bulk_density`, `porosity`: the model properties; `config`: all parameters of the run as a JSON string
```python
import numpy as np, torch
from torch_geometric.data import Data

g = np.load("meshes/wood_cube/wood_cube-graph.npz")
data = Data(x=torch.from_numpy(g["x"]).float(), edge_index=torch.from_numpy(g["edge_index"]),
            edge_attr=torch.from_numpy(g["edge_attr"]).float(), pos=torch.from_numpy(g["pos"]).float(),
            y=torch.tensor([float(g["porosity"])]))
```
The same arrays are returned by `ringspy.graph_arrays(result)` for any `Result`, and held in `result.graph` when the flag is on.

## Contributing
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04945/status.svg)](https://doi.org/10.21105/joss.04945)

If you make use of this code, or find it is useful for your own research , please cite the JOSS paper:
```
@article{ringspy,
    doi = {10.21105/joss.04945},
    url = {https://doi.org/10.21105/joss.04945},
    year = {2023},
    publisher = {The Open Journal},
    volume = {8},
    number = {83},
    pages = {4945},
    author = {Hao Yin and Gianluca Cusatis},
    title = {RingsPy: A Python package for Voronoi mesh generation of cellular solids with radial growth pattern},
    journal = {Journal of Open Source Software}
}
```

Contributions are always welcome!

If you wish to contribute code/algorithms to this project, or to propose a collaboration study, please send an email to haoyin2022 [at] u.northwestern.edu .
## License
![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)

Distributed under the GPL v3 license.


