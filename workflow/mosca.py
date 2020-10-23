#!/usr/bin/env python
# -*- coding: utf-8 -*-

import snakemake
import argparse
import multiprocessing

__version__ = '1.2.0'

parser = argparse.ArgumentParser(description="MOSCA's main script")
parser.add_argument("-s", "--snakefile", type=str, default="Snakefile", help="Snakefile file")
parser.add_argument("-c", "--configfile", type=str, help="Configuration file for MOSCA (JSON or YAML)",
                    default='config.json')
parser.add_argument('-v', '--version', action='version', version='MOSCA ' + __version__)
args = parser.parse_args()

snakemake.main("-s {} --cores {} --configfile {}".format(args.snakefile, multiprocessing.cpu_count(), args.configfile))
