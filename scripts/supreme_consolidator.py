#!/usr/bin/env python

import argparse
import supreme

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('Consolidate survey property maps for a '
                                                  'combination of tracts/filters/map types'))

    parser.add_argument('-c', '--configfile', action='store', type=str, required=True,
                        help='YAML config file')
    parser.add_argument('-o', '--outputpath', action='store', type=str, required=False,
                        default='./', help='Name of output path.')
    parser.add_argument('-b', '--outputbase', action='store', type=str, required=False,
                        help='Name of output base filename of consolidated maps.')
    parser.add_argument('-t', '--tracts', action='store', type=str, required=False,
                        help='Name of tract(s), ^ delimited')
    parser.add_argument('-f', '--filters', action='store', type=str, required=True,
                        help='Name of filter(s), ^ delimited')
    parser.add_argument('-n', '--nside_coverage', action='store', type=int, required=False,
                        help=('Nside of output coverage map; if not set, use same as '
                              'tract maps.'))
    parser.add_argument('-k', '--clobber', action='store_true', required=False,
                        help='Clobber any existing files')

    args = parser.parse_args()

    config = supreme.Configuration.load_yaml(args.configfile)

    if args.tracts is not None:
        tracts = [int(t) for t in args.tracts.split('^')]
    else:
        tracts = None
    filters = [f for f in args.filters.split('^')]

    consolidator = supreme.TractConsolidator(config, args.outputpath)
    _ = consolidator(filters, tracts=tracts, clobber=args.clobber,
                     nside_coverage=args.nside_coverage, outputbase=args.outputbase)
