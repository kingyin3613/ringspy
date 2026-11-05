# -*- coding: utf-8 -*-
"""
Command line interface: ``ringspy run wood_cube.toml``, ``ringspy batch wood_cube.toml --seeds 0-9`` and
``ringspy outline shape.svg``
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
    run.add_argument('--preview', action='store_true', help='only the 2D cell mesh, quickly: stop after the '
                     'tessellation and clipping and save <geoName>.png (the walls coloured by growth ring), without '
                     'the extrusion and the mesh files (previewFlag = true in the file does the same)')

    rec = commands.add_parser('reconstruct', help='reconstruct the lattice of a specimen from a 3D image (e.g. a CT '
                              'scan), with the parameters of a TOML file',
                              description='Reconstruct the lattice (grains along the cell corners, connectors across '
                              'the cell walls) of a cellular specimen from a 3D image, with the parameters of a TOML '
                              'file. The files go to <path>/<geoName>/ (relative to the current folder).')
    rec.add_argument('config', help='TOML file, e.g. examples/reconstruction/wood_ct.toml')
    rec.add_argument('--volume', help='3D image file (.npy, .npz or .tif), overrides the one in the TOML file')
    rec.add_argument('--path', help='output folder, overrides the one in the file')
    rec.add_argument('--geoName', '--geo-name', dest='geoName', help='geometry name, overrides the one in the file')

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
    batch.add_argument('--preview', action='store_true', help='only the 2D cell mesh of every run (as ringspy run '
                       '--preview): the table has the statistics of the cells but no mass, volume or density')

    outline = commands.add_parser('outline', help='read a clipping outline from an SVG or DXF file and print it as '
                                  'box_polygon vertices',
                                  description='Read the clipping polygon of an SVG or DXF file (as box_polygon does), '
                                  'optionally simplified, and print its size and its vertices as a TOML line.')
    outline.add_argument('file', help='SVG or DXF file')
    outline.add_argument('--tolerance', type=float, help='Douglas-Peucker simplification tolerance (mm), e.g. merge_tol')
    outline.add_argument('--scale', type=float, default=1.0, help='factor on the coordinates (after the units -> mm)')
    outline.add_argument('--center', type=float, nargs=2, metavar=('X', 'Y'), help='move the outline to this center')

    args = parser.parse_args(argv)

    if args.command == 'outline':
        from .outline import outline_toml, read_outline
        try:
            points = read_outline(args.file, tolerance=args.tolerance, scale=args.scale, center=args.center)
        except (OSError, ValueError) as e:
            parser.exit(2, 'ringspy: error: {}\n'.format(e))
        lo, hi = points.min(axis=0), points.max(axis=0)
        print('# {}: {:d} vertices, x from {:g} to {:g}, y from {:g} to {:g} (mm)'.format(
            args.file, len(points), lo[0], hi[0], lo[1], hi[1]))
        print(outline_toml(points))
        return 0

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

    if args.command == 'batch':
        from .batch import batch as run_batch, parse_seeds, parse_set
        try:
            config = Config.from_toml(args.config)
            changes = {k: getattr(args, k) for k in ['path', 'geoName'] if getattr(args, k) is not None}
            if args.preview:
                changes['previewFlag'] = True
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
        if args.preview:
            changes['previewFlag'] = True
        config = config.replace(**changes)
    except (OSError, ValueError, TypeError) as e:
        parser.exit(2, 'ringspy: error: {}\n'.format(e))
    result = generate(config)
    print('{} written to {}'.format('Preview' if config.previewFlag else 'Mesh files', result.output_dir))
    return 0


if __name__ == '__main__':
    sys.exit(main())
