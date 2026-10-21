"""Generate the honeycomb_tri example from Python: same as `ringspy run honeycomb_tri.toml`"""
from pathlib import Path

import ringspy

if __name__ == "__main__":
    config = ringspy.load_config(Path(__file__).with_suffix(".toml"))
    # parameters can be changed here, e.g. config = config.replace(seed=1, stlFlag=False)
    ringspy.generate(config)
