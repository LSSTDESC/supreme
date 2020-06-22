#!/usr/bin/env python

import argparse
import supreme

import lsst.daf.persistence as dafPersist

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('Make observing conditions maps for a '
                                                  'combination of tract/patch/filters'))

    parser.add_argument('-c', '--configfile', action='store', type=str, required=True,
                        help='YAML config file')
    parser.add_argument('-r', '--repo', action='store', type=str, required=True,
                        help='Name of input repo')
    parser.add_argument('-o', '--outputpath', action='store', type=str, required=False,
                        default='./', help='Name of output path.')
    parser.add_argument('-t', '--tracts', action='store', type=str, required=True,
                        help='Name of tract(s)')
    parser.add_argument('-f', '--filters', action='store', type=str, required=True,
                        help='Name of filter(s), ^ delimited')
    parser.add_argument('-p', '--patches', action='store', type=str, required=False,
                        help='Name of patch(es), ^ delimited')
    parser.add_argument('-C', '--cores', action='store', type=int, required=True,
                        default=1, help='Number of cores to run on')
    parser.add_argument('-P', '--individual_patches', action='store_true', required=False,
                        help='Save individual patches instead of consolidating tracts')
    parser.add_argument('-k', '--clobber', action='store_true', required=False,
                        help='Clobber any existing files')

    args = parser.parse_args()

    butler = dafPersist.Butler(args.repo)
    config = supreme.Configuration.load_yaml(args.configfile)

    tracts = [int(t) for t in args.tracts.split('^')]
    filters = [f for f in args.filters.split('^')]
    if args.patches is not None:
        patches = [p for p in args.patches.split('^')]
    else:
        patches = None

    mapper = supreme.MultiMapper(butler, config, args.outputpath, ncores=args.cores)
    mapper(tracts, filters, patches=patches, consolidate=not args.individual_patches,
           clobber=args.clobber)
