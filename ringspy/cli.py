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

    rec = commands.add_parser('reconstruct', help='reconstruct the lattice of a specimen from a 3D image (e.g. a CT '
                              'scan), with the parameters of a TOML file',
                              description='Reconstruct the lattice (grains along the cell corners, connectors across '
                              'the cell walls) of a cellular specimen from a 3D image, with the parameters of a TOML '
                              'file. The files go to <path>/<geoName>/ (relative to the current folder).')
    rec.add_argument('config', help='TOML file, e.g. examples/reconstruction/wood_ct.toml')
    rec.add_argument('--volume', help='3D image file (.npy, .npz or .tif), overrides the one in the TOML file')
    rec.add_argument('--path', help='output folder, overrides the one in the file')
    rec.add_argument('--geoName', '--geo-name', dest='geoName', help='geometry name, overrides the one in the file')

    args = parser.parse_args(argv)

    import matplotlib
    matplotlib.use('Agg')  # only files are written, no windows

    if args.command == 'reconstruct':
        from .config import ReconstructionConfig
        from .reconstruction import reconstruct
        try:
            config = ReconstructionConfig.from_toml(args.config)
            changes = {k: getattr(args, k) for k in ['volume', 'path', 'geoName'] if getattr(args, k) is not None}
            config = config.replace(**changes)
        except (OSError, ValueError, TypeError) as e:
            parser.exit(2, 'ringspy: error: {}\n'.format(e))
        result = reconstruct(config)
        print('Mesh files written to {}'.format(result.output_dir))
        return 0

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
