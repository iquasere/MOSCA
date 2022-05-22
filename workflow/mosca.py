#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pathlib

import snakemake
import argparse
import sys
import json
import yaml
from time import gmtime, strftime, time

__version__ = '1.6.0'

parser = argparse.ArgumentParser(description="MOSCA's main script")
parser.add_argument("-s", "--snakefile", default=f'{sys.path[0]}/Snakefile', help="Snakefile file")
parser.add_argument(
    "-c", "--configfile", default='config.json',
    help="Configuration file for MOSCA (JSON or YAML). Obtain one at https://iquasere.github.io/MOSGUITO")
parser.add_argument(
    '--unlock', action='store_true', default=False,
    help='If user forced termination of workflow, this might be required')
parser.add_argument('-v', '--version', action='version', version=f'MOSCA {__version__}')
args = parser.parse_args()


def read_config(filename):
    if filename.split('.')[-1] == 'yaml':
        with open(filename) as stream:
            try:
                return yaml.safe_load(stream), 'yaml'
            except yaml.YAMLError as exc:
                print(exc)
    elif filename.split('.')[-1] == 'json':
        with open(filename) as f:
            return json.load(f), 'json'
    else:
        exit('Config file must end in either ".json" or ".yaml"')


def save_config(config_data, filename, output_format):
    with open(filename, 'w', encoding='utf-8') as f:
        if output_format == 'json':
            json.dump(config_data, f, ensure_ascii=False, indent=2)
        elif output_format == 'yaml':
            yaml.dump(config_data, f, ensure_ascii=False, indent=2)
        else:
            return NotImplementedError


def human_time(seconds):
    days = seconds // 86400
    if days > 0:
        return strftime(f"{days}d%Hh%Mm%Ss", gmtime(seconds))
    return strftime("%Hh%Mm%Ss", gmtime(seconds))


start_time = time()
config, config_format = read_config(args.configfile)
pathlib.Path(config["output"]).mkdir(parents=True, exist_ok=True)
save_config(config, f'{config["output"]}/config.json', output_format=config_format)

snakemake.main(
    f"-s {args.snakefile} --printshellcmds --cores {config['threads']} --configfile {args.configfile}"
    f"{' --unlock' if args.unlock else ''}")

print(f'MOSCA analysis finished in {human_time(time() - start_time)}')
