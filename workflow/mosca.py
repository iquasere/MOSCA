#!/usr/bin/env python

import pathlib
import snakemake
import argparse
import sys
import json
import yaml
import pandas as pd
import re

__version__ = '2.3.0'

parser = argparse.ArgumentParser(description="MOSCA's main script")
parser.add_argument("-s", "--snakefile", default=f'{sys.path[0]}/Snakefile', help="Path to Snakefile")
parser.add_argument(
    "-c", "--configfile", required=True,
    help="Configuration file for MOSCA (JSON or YAML). Obtain one at https://iquasere.github.io/MOSGUITO")
parser.add_argument("--use-singularity", action="store_true", default=False, help="Use singularity")
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
        sys.exit('ERROR: Config file must end in either ".json" or ".yaml"')


def save_config(config_data, filename, output_format):
    with open(filename, 'w', encoding='utf-8') as f:
        if output_format == 'json':
            json.dump(config_data, f, ensure_ascii=False, indent=2)
        elif output_format == 'yaml':
            yaml.dump(config_data, f, indent=2)
        else:
            return NotImplementedError


def validate_exps(exps_data):
    exps = pd.DataFrame(exps_data)
    reserved_words = [
        'if', 'else', 'repeat', 'while', 'function', 'for', 'in', 'next', 'break', 'TRUE', 'FALSE', 'NULL', 'Inf',
        'NaN', 'NA', 'NA_integer_', 'NA_real_', 'NA_complex_', 'NA_character_']
    good_pattern = re.compile(r'^(?!^\d)(?!^\.\d)([\w.]+)$')    # don't start with a number, nor a decimal (dot followed by number)
    for name in exps['Name']:
        if not name:        # if name is None, or empty string, MOSCA should be able to build one that is fine
            continue
        if name in reserved_words:
            sys.exit(f'INVALID "NAME" in "experiments": {name} is a reserved R word.')
        if not bool(good_pattern.match(name)):
            sys.exit(f'INVALID "NAME" in "experiments": {name} starts with a number or has a special character.')


def validate_config(config_data):
    if not config_data['do_assembly'] and config_data['do_binning']:
        sys.exit('ERROR: Can only do binning if assembly is performed.')
    validate_exps(config["experiments"])


user_config, config_format = read_config(args.configfile)
config = read_config(f'{sys.path[0]}/default_config.json')[0]       # default configurations
for key in config.keys():
    if key in user_config.keys():
        config[key] = user_config[key]                              # set default values
validate_config(config)
pathlib.Path(config["output"]).mkdir(parents=True, exist_ok=True)
save_config(config, f'{config["output"]}/config.json', output_format=config_format)

command = (
    f"-s {args.snakefile} --printshellcmds --cores {config['threads']} --configfile {config['output']}/config.json "
    f"--use-conda{' --use-singularity' if args.use_singularity else ''}{' --unlock' if args.unlock else ''}")

print(f"MOSCA command: snakemake {command}")
snakemake.main(command)
