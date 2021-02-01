#!/usr/bin/env python
# -*- coding: utf-8 -*-

import snakemake
import argparse
import multiprocessing
import sys

__version__ = '1.3.0'

parser = argparse.ArgumentParser(description="MOSCA's main script")
parser.add_argument("-s", "--snakefile", type=str, default="{}/Snakefile".format(sys.path[0]), help="Snakefile file")
parser.add_argument("-c", "--configfile", type=str, help="Configuration file for MOSCA (JSON or YAML)",
                    default='config.json')
parser.add_argument('--unlock', action='store_true', default=False,
                    help='If user forced termination of workflow, this might be required')
parser.add_argument('-v', '--version', action='version', version='MOSCA ' + __version__)
args = parser.parse_args()

snakemake.main("-s {} --printshellcmds --cores {} --configfile {}{}".format(
    args.snakefile, multiprocessing.cpu_count(), args.configfile, ' --unlock' if args.unlock else ''))
