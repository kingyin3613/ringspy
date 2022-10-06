"""
the magick behind ``pip install ...``
"""
from setuptools import setup
import pathlib

# The directory containing this file
HERE = pathlib.Path(__file__).parent

if __name__ == "__main__":
    setup()