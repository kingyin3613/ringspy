#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys


__version__ = '1.2.0'
VERSION = __version__

from .config import Config, load_config  # noqa: E402 (MeshGenTools imports __version__)
from .pipeline import Result, generate  # noqa: E402

__all__ = ['Config', 'Result', 'generate', 'load_config']
