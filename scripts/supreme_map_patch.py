#!/usr/bin/env python

import argparse
import supreme

import lsst.daf.persistence as dafPersist

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make observing conditions maps for a patch')

    parser.add_argument('-c', '--configfile', action='store', type=str, required=True,
                        help='YAML config file')
    parser.add_argument('-r', '--repo', action='store', type=str, required=True,
                        help='Name of input repo')
    parser.add_argument('-o', '--outputpath', action='store', type=str, required=False,
                        default='./', help='Name of output path.')
    parser.add_argument('-t', '--tract', action='store', type=int, required=True,
                        help='Name of tract')
    parser.add_argument('-p', '--patch', action='store', type=str, required=True,
                        help='Name of patch')
    parser.add_argument('-f', '--filter_name', action='store', type=str, required=True,
                        help='Name of filter')

    args = parser.parse_args()

    butler = dafPersist.Butler(args.repo)
    config = supreme.Configuration(args.configfile)

    mapper = supreme.PatchMapper(butler, config, args.outputpath)
    mapper.run(args.filter_name, args.tract, args.patch, return_values_list=False)
