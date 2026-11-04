"""Generate the wood_gradual example from Python: same as `ringspy run wood_gradual.toml`"""
from pathlib import Path

import ringspy

if __name__ == "__main__":
    config = ringspy.load_config(Path(__file__).with_suffix(".toml"))
    # parameters can be changed here, e.g. config = config.replace(transition="linear", transition_width=1.0)
    ringspy.generate(config)
