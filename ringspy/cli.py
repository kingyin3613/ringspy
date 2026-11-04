# -*- coding: utf-8 -*-
"""
Command line interface: ``ringspy run wood_cube.toml`` and ``ringspy batch wood_cube.toml --seeds 0-9``
"""

import argparse
import sys

from . import __version__


def main(argv=None):
    parser = argparse.ArgumentParser(prog='ringspy', description='RingsPy: Voronoi mesh generation of cellular solids '
                                     'with radial growth patterns.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    commands = parser.add_subparsers(dest='command', metavar='command')
    commands.required = True

    run = commands.add_parser('run', help='generate a mesh from a TOML file of input parameters',
                              description='Generate a mesh from a TOML file of input parameters. '
                              'The files go to <path>/<geoName>/ (relative to the current folder).')
    run.add_argument('config', help='TOML file, e.g. examples/wood_cube.toml')
    run.add_argument('--seed', type=int, help='random seed, overrides the one in the file')
    run.add_argument('--path', help='output folder, overrides the one in the file')
    run.add_argument('--geoName', '--geo-name', dest='geoName', help='geometry name, overrides the one in the file')

    batch = commands.add_parser('batch', help='generate many meshes over seeds and parameter sweeps, with a table '
                                'of their descriptors',
                                description='Generate one mesh per combination of the --set values and per seed, '
                                'from a TOML file of input parameters. Run i goes to <path>/<geoName>_<i>/, and '
                                '<path>/<geoName>-batch.csv lists the runs with their parameters and descriptors '
                                '(number of cells, vertices, walls and beam elements, mass, bulk volume and density, '
                                'porosity). Example: ringspy batch wood_cube.toml --seeds 0-9 '
                                '--set cellsize_sparse=0.02,0.03 --set stlFlag=false --jobs 4')
    batch.add_argument('config', help='TOML file, the base of all runs')
    batch.add_argument('--seeds', action='append', default=[], metavar='SEEDS',
                       help='random seeds of the runs: integers and ranges like 0-9, separated by commas '
                       '(may be repeated); without it, all runs keep the seed of the file')
    batch.add_argument('--set', action='append', default=[], metavar='NAME=V1,V2,...',
                       help='values of a parameter to sweep (TOML syntax, lists in brackets); one value just '
                       'changes the parameter in all runs; may be repeated, all combinations are generated')
    batch.add_argument('--jobs', '-j', type=int, default=1, help='runs in parallel (default 1)')
    batch.add_argument('--path', help='output folder, overrides the one in the file')
    batch.add_argument('--geoName', '--geo-name', dest='geoName', help='name of the batch, overrides the one in the file')

    args = parser.parse_args(argv)

    import matplotlib
    matplotlib.use('Agg')  # only files are written, no windows
    from .config import Config
    from .pipeline import generate

    if args.command == 'batch':
        from .batch import batch as run_batch, parse_seeds, parse_set
        try:
            config = Config.from_toml(args.config)
            changes = {k: getattr(args, k) for k in ['path', 'geoName'] if getattr(args, k) is not None}
            config = config.replace(**changes)
            seeds = parse_seeds(args.seeds) if args.seeds else None
            sweep = parse_set(args.set)
            if seeds is not None and not seeds:
                raise ValueError('--seeds: give at least one seed')
            results = run_batch(config, seeds=seeds, sweep=sweep, jobs=args.jobs)
        except (OSError, ValueError, TypeError) as e:
            parser.exit(2, 'ringspy: error: {}\n'.format(e))
        return 0 if all(r is not None for r in results) else 1

    try:
        config = Config.from_toml(args.config)
        changes = {k: getattr(args, k) for k in ['seed', 'path', 'geoName'] if getattr(args, k) is not None}
        config = config.replace(**changes)
    except (OSError, ValueError, TypeError) as e:
        parser.exit(2, 'ringspy: error: {}\n'.format(e))
    result = generate(config)
    print('Mesh files written to {}'.format(result.output_dir))
    return 0


if __name__ == '__main__':
    sys.exit(main())
