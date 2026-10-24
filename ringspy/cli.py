# -*- coding: utf-8 -*-
"""
Command line interface: ``ringspy run wood_cube.toml``
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

    args = parser.parse_args(argv)

    import matplotlib
    matplotlib.use('Agg')  # only files are written, no windows
    from .config import Config
    from .pipeline import generate

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
