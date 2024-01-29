# -*- coding: utf-8 -*-
"""
MOSCA's script for quantifying reads

By João Sequeira

Dec 2022
"""

import pandas as pd
from mosca_tools import perform_alignment, normalize_counts_by_size


def quantification(exps: pd.DataFrame, output: str, sample: str, did_assembly: bool = True) -> tuple:
    """
    Perform quantification of reads with contigs as reference
    :param exps: DataFrame with the experiments
    :param output: Output directory
    :param sample: Sample name
    :param did_assembly: Whether assembly was performed or not
    """
    mg_result = mg_result_norm = pd.DataFrame(columns=['Contig'])
    mt_result = mt_result_norm = pd.DataFrame(columns=['Gene'])
    pexps = exps[(exps['Sample'] == sample)]
    for i in pexps.index:
        if pexps.loc[i]['Data type'] == 'mrna' or not did_assembly:
            reference = f"{output}/Annotation/{pexps.loc[i]['Sample']}/fgs.ffn"
        elif pexps.loc[i]['Data type'] == 'dna':
            reference = f"{output}/Assembly/{pexps.loc[i]['Sample']}/contigs.fasta"
        else:
            continue
        if ',' in pexps.loc[i]['Files']:
            reads = [f"{output}/Preprocess/Trimmomatic/quality_trimmed_{pexps.loc[i]['Name']}_{fr}_paired.fq"
                     for fr in ['forward', 'reverse']]
        else:
            reads = [f"{output}/Preprocess/Trimmomatic/quality_trimmed_{pexps.loc[i]['Name']}.fq"]
        perform_alignment(
            reference, reads, f"{output}/Quantification/{pexps.loc[i]['Name']}", threads=snakemake.threads)
        counts = pd.read_csv(
            f"{output}/Quantification/{pexps.loc[i]['Name']}.readcounts",
            sep='\t', names=['Gene' if pexps.loc[i]['Data type'] == 'mrna' else 'Contig', pexps.loc[i]['Name']])
        normalize_counts_by_size(
            f"{output}/Quantification/{pexps.loc[i]['Name']}.readcounts", reference)
        # Read the results of alignment and add them to the readcounts result at sample level
        normalized_counts = pd.read_csv(
            f"{output}/Quantification/{pexps.loc[i]['Name']}.readcounts.norm",
            sep='\t', names=['Gene' if pexps.loc[i]['Data type'] == 'mrna' else 'Contig', pexps.loc[i]['Name']])
        if pexps.loc[i]['Data type'] == 'dna':
            mg_result = pd.merge(mg_result, counts, how='outer', on='Contig')
            mg_result_norm = pd.merge(mg_result_norm, normalized_counts, how='outer', on='Contig')
        else:
            mt_result = pd.merge(mt_result, counts, how='outer', on='Gene')
            mt_result_norm = pd.merge(mt_result_norm, normalized_counts, how='outer', on='Gene')
    return mg_result, mg_result_norm, mt_result, mt_result_norm


def run():
    exps = pd.read_csv(snakemake.params.exps, sep='\t')
    did_assembly=snakemake.params.did_assembly

    for sample in set(exps['Sample']):
        mg_result, mg_result_norm, mt_result, mt_result_norm = quantification(
                exps, snakemake.params.output, sample, did_assembly=did_assembly)
        if len(mg_result) > 0:
            mg_result.to_csv(
                f"{snakemake.params.output}/Quantification/{sample}_mg.readcounts", sep='\t', index=False)
            if did_assembly:
                mg_result_norm.to_csv(
                    f"{snakemake.params.output}/Quantification/{sample}_mg_norm.tsv", sep='\t', index=False)
        if len(mt_result) > 0:
            mt_result.to_csv(
                f"{snakemake.params.output}/Quantification/{sample}_mt.readcounts", sep='\t', index=False)
            if did_assembly:
                mt_result_norm.astype(int, errors='ignore').to_csv(
                    f"{snakemake.params.output}/Quantification/{sample}_mt_norm.tsv", sep='\t', index=False)


if __name__ == '__main__':
    run()
