#!/usr/bin/env python

import argparse
import supreme

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('Make mask from ds9 region file(s)'))

    parser.add_argument('-c', '--configfile', action='store', type=str, required=True,
                        help='YAML config file')
    parser.add_argument('-r', '--regionfiles', action='store', type=str, required=True,
                        help='Name of region files, ^-delimited')
    parser.add_argument('-m', '--maskfile', action='store', type=str, required=True,
                        help='Output mask file')
    parser.add_argument('-k', '--clobber', action='store_true', required=False,
                        help='Clobber any existing files')

    args = parser.parse_args()

    config = supreme.RegionConfiguration.load_yaml(args.configfile)

    regionfiles = [f for f in args.regionfiles.split('^')]

    masker = supreme.RegionMasker(config)
    masker(args.maskfile, regionfiles, clobber=args.clobber)
