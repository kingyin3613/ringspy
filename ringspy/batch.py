# -*- coding: utf-8 -*-
"""
Batch generation: many meshes from one configuration, over random seeds and parameter
sweeps, with one table of their descriptors (``ringspy batch``)

author: king_yin3613
email: haoyin2022@u.northwestern.edu
"""

import csv
import itertools
import json
import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from .config import Config

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

from .stats import SUMMARY_FIELDS

# columns of the table after the run, its name, its seed and the swept parameters: attributes of the Result,
# then the descriptors of its cells (Result.summary())
DESCRIPTORS = ['nParticles', 'nvertex', 'nridge', 'nbeamElem', 'mass', 'bulk_volume', 'bulk_density', 'porosity']


def batch(config, seeds=None, sweep=None, jobs=1, **changes):
    """Generate one mesh per combination of the sweep values and per seed, and write a table of
    their descriptors.

    config :: a Config, a TOML file name, or a dict of parameters: the base of all runs.
    seeds :: the random seeds of the runs (e.g. ``range(10)``); None: every run keeps the seed of config.
    sweep :: {parameter: [values, ...]}: all combinations of the values are generated (the first
        parameter varies slowest, the seeds fastest). A parameter with one value is simply changed.
    jobs :: number of runs in parallel (separate processes).
    changes :: parameters to override in all runs, e.g. ``stlFlag=False``.

    Run i is named ``<geoName>_<i:03d>`` and its files go to ``<path>/<geoName>_<i:03d>/``, with
    ``path`` and ``geoName`` those of config. The table ``<path>/<geoName>-batch.csv`` has one row
    per run: run, geoName, seed, the swept parameters, then nParticles, nvertex, nridge, nbeamElem,
    mass, bulk_volume, bulk_density, porosity, output_dir, and error (the message of a run that
    failed; the other runs are still generated).

    Returns the list of Results, in the order of the table (None for a failed run).
    """
    if isinstance(config, (str, os.PathLike)):
        config = Config.from_toml(config)
    elif isinstance(config, dict):
        config = Config.from_dict(config)
    elif not isinstance(config, Config):
        raise TypeError("config must be a Config, a TOML file name or a dict, not {}".format(type(config).__name__))
    if changes:
        config = config.replace(**changes)
    sweep = dict(sweep or {})
    for name, values in sweep.items():
        if isinstance(values, (str, bytes)) or not hasattr(values, '__len__') or len(values) == 0:
            raise ValueError("sweep values of {!r} must be a non-empty list, not {!r}".format(name, values))
    config.replace(**{name: values[0] for name, values in sweep.items()})  # the names are checked once, early
    if jobs < 1:
        raise ValueError("jobs must be at least 1")

    runs = []  # (run index, parameters of the run, swept values of the run)
    seed_values = [None] if seeds is None else list(seeds)
    for i, (values, seed) in enumerate(itertools.product(itertools.product(*sweep.values()), seed_values)):
        swept = dict(zip(sweep.keys(), values))
        params = dict(config.to_dict(), geoName='{}_{:03d}'.format(config.geoName, i), **swept)
        if seeds is not None:
            params['seed'] = seed
        runs.append((i, params, swept))
    table = Path(config.path, config.geoName + '-batch.csv')
    os.makedirs(table.parent, exist_ok=True)

    results = [None]*len(runs)
    rows = [None]*len(runs)
    if jobs == 1:
        outcomes = map(_run, runs)
    else:
        executor = ProcessPoolExecutor(max_workers=min(jobs, len(runs)))
        outcomes = executor.map(_run, runs)
    try:
        for i, params, swept, result, error in outcomes:
            results[i] = result
            rows[i] = _row(i, params, swept, result, error)
            print('run {:d}/{:d} ({}): {}'.format(i + 1, len(runs), params['geoName'],
                                                  'failed, ' + error.splitlines()[-1] if error else 'done'))
    finally:
        if jobs > 1:
            executor.shutdown()
    with open(table, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['run', 'geoName', 'seed'] + list(sweep) + DESCRIPTORS + SUMMARY_FIELDS
                                + ['output_dir', 'error'])
        writer.writeheader()
        writer.writerows(rows)
    print('Batch table written to {}'.format(table))
    return results


def _run(run):
    """One run of the batch, in this or another process: (i, params, swept, Result or None, error text or None)"""
    i, params, swept = run
    import matplotlib
    matplotlib.use('Agg')  # only files are written, no windows (also in worker processes)
    from .pipeline import generate
    try:
        return i, params, swept, generate(params), None
    except Exception:  # the other runs of the batch go on
        return i, params, swept, None, traceback.format_exc()


def _row(i, params, swept, result, error):
    row = {'run': i, 'geoName': params['geoName'], 'seed': params['seed']}
    row.update({name: _cell(value) for name, value in swept.items()})
    if result is not None:
        row.update({name: getattr(result, name) for name in DESCRIPTORS})
        row.update(result.summary())
        row['output_dir'] = str(result.output_dir)
    row['error'] = error.strip().splitlines()[-1] if error else ''
    return row


def _cell(value):
    """A parameter value in the table: lists as JSON, so that they read back"""
    return json.dumps(value) if isinstance(value, (list, tuple)) else value


def parse_seeds(texts):
    """Seeds of the command line: integers and ranges a-b (both included), separated by commas"""
    seeds = []
    for text in texts:
        for item in text.split(','):
            item = item.strip()
            if not item:
                continue
            lo, sep, hi = item.partition('-')
            try:
                seeds.extend(range(int(lo), int(hi) + 1) if sep else [int(item)])
            except ValueError:
                raise ValueError("--seeds takes integers and ranges like 0-9, not {!r}".format(item)) from None
    return seeds


def parse_set(texts):
    """Sweeps of the command line: 'name=value1,value2,...' with the values in TOML syntax (lists in
    brackets, strings that are not numbers or true/false may be left unquoted)"""
    sweep = {}
    for text in texts:
        name, sep, values = text.partition('=')
        name = name.strip()
        if not sep or not name:
            raise ValueError("--set takes name=value1,value2,..., not {!r}".format(text))
        try:
            parsed = tomllib.loads('v = [{}]'.format(values))['v']
        except tomllib.TOMLDecodeError:
            parsed = [v.strip() for v in values.split(',')]  # plain strings, e.g. box_shape=square,hexagon
        if not parsed:
            raise ValueError("--set {}: give at least one value".format(name))
        sweep[name] = parsed
    return sweep
