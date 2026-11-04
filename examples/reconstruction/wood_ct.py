"""Reconstruct the lattice of a wood cube from a 3D image of its cells: same as `ringspy reconstruct wood_ct.toml`

A synthetic scan is made first (if wood_ct.npy is not there yet): the cells of the wood_cube example
are extruded along fibers that lean, twist and wobble, about a third of them ending inside the cube,
rasterised into voxels of 2 um, blurred and given noise, as a CT scan of a 0.6 mm cube would look. With a real scan, save it as a .npy/.npz array
or a .tif stack, set `volume` and `voxel_size` in the TOML file, and skip this step.
"""
from pathlib import Path

import numpy as np

import ringspy

HERE = Path(__file__).parent

if __name__ == "__main__":
    scan_file = HERE / "wood_ct.npy"
    if not scan_file.exists():
        cells = ringspy.load_config(HERE.parent / "wood_cube.toml").replace(
            r_max=2.0, nrings=3, box_size=0.6, box_center=[1.1, 0.0], z_max=0.6, iter_max=4000, seed=2023)
        scan = ringspy.voxelize(cells, voxel_size=0.002, tilt=[0.05, 0.02], twist=0.02,
                                wobble_amplitude=0.006, wobble_wavelength=0.4, noise=0.15, blur=0.8, seed=2023,
                                fiber_length=1.5, fiber_length_cov=0.4)
        np.save(scan_file, np.clip(255*scan.volume, 0, 255).astype(np.uint8))  # 8-bit grayscale, 300^3 voxels
        print("Synthetic scan written to", scan_file)
    config = ringspy.ReconstructionConfig.from_toml(HERE / "wood_ct.toml")
    # parameters can be changed here, e.g. config = config.replace(nsegments=2, stlFlag=False)
    result = ringspy.reconstruct(config)
    print("Grain angle {:.2f} degrees, {} cells, {} fibers, {} beam elements".format(
        result.grain_angle, result.nParticles, result.nfibers, result.nbeamElem))
