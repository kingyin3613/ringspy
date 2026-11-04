"""Generate the wood_cube_fibers example from Python: same as `ringspy run wood_cube_fibers.toml`"""
from pathlib import Path

import ringspy

if __name__ == "__main__":
    config = ringspy.load_config(Path(__file__).with_suffix(".toml"))
    # parameters can be changed here, e.g. config = config.replace(fiber_length=1.0, seed=1)
    ringspy.generate(config)
