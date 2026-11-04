#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys


__version__ = '1.3.0'
VERSION = __version__

from .config import Config, load_config  # noqa: E402 (MeshGenTools imports __version__)
from .pipeline import Result, generate, graph_arrays  # noqa: E402
from .batch import batch  # noqa: E402

__all__ = ['Config', 'Result', 'batch', 'generate', 'graph_arrays', 'load_config']
