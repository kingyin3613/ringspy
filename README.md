# RingsPy
[![python version](https://upload.wikimedia.org/wikipedia/commons/a/a5/Blue_Python_3.8_Shield_Badge.svg)
![PyPI version](https://badge.fury.io/py/RingsPy.svg)](https://badge.fury.io/py/RingsPy)
![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)
![Tests](https://github.com/kingyin3613/RingsPy/actions/workflows/tests.yml/badge.svg)

RingsPy is a Voronoi diagrams-based geometric generation tool that generates 3D meshes and models of prismatic cellular solids with radial growth rules.

## Dependencies and Installation

### 1. pip install
To install RingsPy, one may use:
```bash
pip install RingsPy
``` 
or use:
```bash
pip install git+https://github.com/kingyin3613/RingsPy.git
``` 
to get updates beyond the latest release. 

RingsPy depends on mainstream Python libraries ``numpy`` and ``scipy``, and optionally depends on library ``hexalattice``, if the regular hexagonal lattice (e.g. honeycomb) is wanted; also ``vtk``, if the 3D STL files are also wanted.

### 2. Installation check

There are some unit tests in [tests](https://github.com/kingyin3613/RingsPy/tree/main/tests/). You can use ``pytest`` to check whether the installation is successful by running this command:

```bash
pytest
```

## Getting Started
Once all required components are installed and one is ready to begin, a path forward should be established for generating the mesh. The basic steps for running/viewing a cellular mesh are listed as the following:

    1. Edit geometry and algorithm parameters
    2. Generate mesh using Mesh Generation Tools
    3. Visualize 2D view using Matplotlib or 3D model in Paraview
    4. (Optional) Export 3D STL model for 3D editing and/or printing 

### 1. Geometry and Parameters
The first step to generate a cellular geometry is selecting geometry and appropriate parameters. 

### 1.1. Geometry
A template file, `test_wood_cube.py` located in the [tests](https://github.com/kingyin3613/RingsPy/tree/main/tests/) directory acts as both the parameter input file, and main executable for the generation of a cubic wood specimen.

*Note: The Mesh Generation Tool by now only accepts many of pre-defined boundary geometries, importing of CAD and/or other 3D model files will be implemented in subsequent versions.*

*Note: for greatest compatibility create the geometry using all millimeters.*

### 1.2. Parameters


By opening a input file, such as `tests/test_wood_cube.py` in any text editor, a file format similar to what is shown below will be displayed:
```
geoName = 'wood_cube'
path = 'meshes'

iter_max = 100
print_interval = 500

# length unit: mm
r_min = 0   # inner radius of wood log
r_max = 4   # outer radius of wood log
nrings = 4 # number of rings
width_heart = 0.3*(r_max-r_min)/nrings # heart wood ring width
width_early = 0.7*(r_max-r_min)/nrings # early wood ring width
width_late = 0.3*(r_max-r_min)/nrings # late wood ring width
log_center = (0,0) # coordinates of log center in the global system of reference
box_center = (1.25,0) # coordinates of box center in the global system of reference
box_size = 1.0 # cube size

# if precracked
x_indent_size = box_size*0.120
y_indent_size = box_size*0.125
x_precrack_size = box_size*0.1
y_precrack_size = box_size*0.02
x_indent = x_min + x_indent_size
y_indent_min = box_center[1] - y_indent_size/2
y_indent_max = box_center[1] + y_indent_size/2
x_precrack = x_indent + x_precrack_size
y_precrack = box_center[1]

cellsize_early = 0.02
cellsize_late = 0.01
cellwallthickness_early = 0.010
cellwallthickness_late = 0.006

merge_operation = 'off'
merge_tol = 0.01

precrackFlag = 'off'
precrack_widths = 0.1

boundaryFlag = 'on'
stlFlag = 'on'
```

- `geoName` is the geometry name, `path` is the folder where the mesh files will be generated.
- `iter_max` is the max number of iteration for randomly placing a new non-overlapping wood cell particle in the 2D annual rings domain.
- `print_interval` is the print interval when every n cell particles are placed in the model domain.
- `r_min` and `r_max` are the upper and lower bounds of radii to generate 2D annual rings, `nrings` is the number of rings.
- `width_heart`, `width_early`, and `width_late`, are annual ring widths for heartwood, earlywood, and latewood, respectively, which all together determine the morphology of the wood mesostructure.
- `log_center`,`box_center`, and `box_size` are for locating the wood log and cutting box.
- `_indent`,`_precrack` parameters are related to precracked sample generation, which will be added in the next version.
- `cellsize_early`,`cellsize_late`, `cellwallthickness_early`, and `cellwallthickness_late` are parameters for the earlywood and latewood cells.
- `merge_operation` flag can be turned on/off for the merging operation, when turned on, all small wood cell ridges shorter than the threshold `merge_tol` will be merged with neighboring ridges.
- `precrackFlag` flag is for inserting a pre-crack, for the notched specimens. So far, only a single line pre-crack with the length of `precrack_widths` is supported.
- `boundaryFlag` flag can be turned on/off for generating neat boundaries consisting of grains.
- `stlFlag` flag can be turned on/off for generating 3D STL files.

![MeshGenerator](<./contents/MeshGenerator.png>)
### 2.1. Run Mesh Generation
Open a Command Prompt or Terminal window and set the current directory to [tests](https://github.com/kingyin3613/RingsPy/tree/main/tests/) or any other directory, then run the command:
```
    python test_wood_cube.py
```
Functions in `MeshGenTools` library will be called to create the initial mesh, wood cell particles following certain cell size distribution will be placed, then `Scipy.voronoi` function will be called to form the initial 2D Voronoi tessellation, additional code reforms the tesselation and generates the desired files. A successful generation will end with the line "Files generated ..." in the Terminal window.

A new folder should have been created in the `.\meshes` directory with the same name as the `geoName` in `test_wood_cube.py`.

### 2.2. Check Mesh and 3D Model Files

The following files should be generated in the `.\meshes\geoName` directory with a successful run:
- Mesh files
    - Non-Uniform Rational B-Splines (NURBS) beam file: `wood_cubeIGA.txt`
    - connector data file: `wood_cube-mesh.txt`
    - Grain-ridge data file: `wood_cube-vertex.mesh`
    - Ridge data file: `wood_cube-ridge.mesh`
- Visualization files
    - Paraview vtk file for initial vertex configuration: `wood_cube_vertices.vtu`
    - Paraview vtk file for initial grain solid configuration: `wood_cube_beams.vtu`
    - Paraview vtk file of initial cell ridge configuration: `wood_cube_conns.vtu`
    - Paraview vtk file of initial connector (volume) configuration: `wood_cube_conns_vol.vtu`
- (Optional) 3D model files
    - STL file of initial cellular solid configuration: `wood_cube.stl`

### 3. Visualization
A scientific visualization application `ParaView` can directly visualize the generated vtk files; It can also visualize generated 3D model STL files if the STL flag is on. The Paraview software can be downloaded from their official website: [https://www.paraview.org/download/](https://www.paraview.org/download/), latest version is recommeded.

### 3.1. Visualize the components of the 3D model in ParaView
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
    7. Scale and position the image as desired
    8. Turn back on Ray Tracing
    9. Adjust Ray Tracing lighting and settings as desired
    10. Export Image
        - File > Save Screenshot
        - Enter a file name > OK
        - Leave new window as-is or increase resolution > OK

### 3.2. Visualize the volumes of the 3D model in ParaView
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

## Contributing

Contributions are always welcome!

If you wish to contribute code/algorithms to this project, or to propose a collaboration study, please send an email to haoyin2022 [at] u.northwestern.edu .

## License
![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)

Distributed under the GPL v3 license. Copyright 2022 Hao Yin.


