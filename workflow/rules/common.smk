from snakemake.remote import FTP
from snakemake.utils import validate
import pathlib
import pandas as pd
import numpy as np
import sys

ftp = FTP.RemoteProvider()

validate(config, schema="../schemas/config.schema.yaml")

OUTPUT = config["output"]
EXPS = pd.DataFrame(config["experiments"])

def set_name(files, data_type):
    filename = files.split('/')[-1]
    if data_type == 'protein':
        return filename # which is the foldername (e.g. input/mp1 -> mp1)
    if ',' in files:
        return filename.split(',')[0].split('_R')[0]
    return filename.split('.fa')[0]

for i in range(len(EXPS)):
    if pd.isnull(EXPS.iloc[i]['Name']) or EXPS.iloc[i]['Name'] == '':
        EXPS.iloc[i, EXPS.columns.get_loc('Name')] = set_name(
            EXPS.iloc[i]['Files'], EXPS.iloc[i]['Data type'])
    if not config['do_assembly']:
        EXPS.iloc[i]['Sample'] = EXPS.iloc[i]['Name']

pathlib.Path(f"{OUTPUT}").mkdir(parents=True, exist_ok=True)
EXPS.to_csv(f"{OUTPUT}/exps.tsv", sep = '\t', index = False)

mg_exps = EXPS[EXPS["Data type"] == 'dna']
mt_exps = EXPS[EXPS["Data type"] == 'mrna']
mp_exps = EXPS[EXPS["Data type"] == 'protein']

if len(mg_exps) == 0 and len(mt_exps) != 0:
    mg_exps = mt_exps
not_mp_exps = EXPS[EXPS["Data type"] != 'protein']
not_mg_exps = EXPS[EXPS["Data type"] != 'dna']

def all_input(wildcards):
    if config['do_assembly']:
        return (
            [f"{OUTPUT}/MOSCA_Protein_Report.xlsx",
            f"{OUTPUT}/MOSCA_Entry_Report.xlsx",
            f"{OUTPUT}/technical_report.tsv",
            f"{OUTPUT}/MOSCA_General_Report.tsv",
            f"{OUTPUT}/MOSCA_results.zip",
            f"{OUTPUT}/KEGG_maps/KEGGCharter_results.tsv"] +
            (expand("{output}/Binning/{sample}/checkm.tsv", output=OUTPUT, sample=set(EXPS['Sample'])) if
                config['do_binning'] else []))
    else:
        return f"{OUTPUT}/MOSCA_Entry_Counts_Report.xlsx"

def join_reads_input(wildcards):
    df = mg_exps[mg_exps['Sample'] == wildcards.sample].reset_index()
    return [f'{OUTPUT}/Preprocess/Trimmomatic/quality_trimmed_{df.iloc[i]["Name"]}{fr}.fq'
           for i in range(len(df))
           for fr in (['_forward_paired', '_reverse_paired'] if ',' in df.iloc[i]["Files"] else [''])]

def fastq2fasta_input(wildcards):
    return expand("{output}/Preprocess/Trimmomatic/quality_trimmed_{name}{fr}.fq", output=OUTPUT,
        fr=(['_forward_paired', '_reverse_paired'] if EXPS["Files"].str.contains(',').tolist() else ''),
        name=wildcards.sample)

def gene_calling_input(wildcards):
    if config['do_assembly']:
        return expand("{output}/Assembly/{sample}/scaffolds.fasta", output=OUTPUT, sample=wildcards.sample)
    return expand(
        "{output}/Preprocess/piled_{name}.fasta", output=OUTPUT, name=wildcards.sample)

def upimapi_input(wildcards):
    if config['do_assembly']:
        return expand("{output}/Annotation/{sample}/aligned.blast", output=OUTPUT, sample=set(EXPS['Sample']))
    return expand("{output}/Annotation/{name}/aligned.blast", output=OUTPUT, name=set(EXPS['Name']))
