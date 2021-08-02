#!/usr/bin/env python
# -*- coding: utf-8 -*-

import snakemake
import argparse
import sys
import json
import yaml
from time import gmtime, strftime, time

__version__ = '1.5.1'

parser = argparse.ArgumentParser(description="MOSCA's main script")
parser.add_argument("-s", "--snakefile", type=str, default=f'{sys.path[0]}/Snakefile', help="Snakefile file")
parser.add_argument("-c", "--configfile", type=str, help="Configuration file for MOSCA (JSON or YAML). "
                                                         "Obtain one at https://iquasere.github.io/MOSGUITO",
                    default='config.json')
parser.add_argument('--unlock', action='store_true', default=False,
                    help='If user forced termination of workflow, this might be required')
parser.add_argument('-v', '--version', action='version', version=f'MOSCA {__version__}')
args = parser.parse_args()


def read_config(filename):
    if filename.split('.')[-1] == 'yaml':
        with open(filename) as stream:
            try:
                return yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
    elif filename.split('.')[-1] == 'json':
        with open(filename) as f:
            return json.load(f)
    else:
        exit('Config file must end in either ".json" or ".yaml"')


start_time = time()
config = read_config(args.configfile)
if config['download_uniprot']:
    config['diamond_database'] = f'{config["resources_directory"]}/uniprot.dmnd'

snakemake.main(f"-s {args.snakefile} --printshellcmds --cores {config['threads']} --configfile {args.configfile}"
               f"{' --unlock' if args.unlock else ''}")

print(f'MOSCA analysis finished in {strftime("%Hh%Mm%Ss", gmtime(time() - start_time))}')
