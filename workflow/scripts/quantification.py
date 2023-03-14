# -*- coding: utf-8 -*-
"""
MOSCA's script for quantifying reads

By João Sequeira

Dec 2022
"""

import pandas as pd
from pathlib import Path
from mosca_tools import perform_alignment, normalize_counts_by_size


def run():
    exps = pd.read_csv(snakemake.params.exps, sep='\t')

    for sample in set(exps['Sample']):
        mg_result = pd.DataFrame(columns=['Contig'])
        mt_result = pd.DataFrame(columns=['Gene'])
        pexps = exps[(exps['Sample'] == sample)]
        for i in pexps.index:
            if pexps.iloc[i]['Data type'] == 'mrna':
                reference = f"{snakemake.params.output}/Annotation/{pexps.iloc[i]['Sample']}/fgs.ffn"
            elif pexps.iloc[i]['Data type'] == 'dna':
                reference = f"{snakemake.params.output}/Assembly/{pexps.iloc[i]['Sample']}/contigs.fasta"
            else:
                continue
            if ',' in pexps.iloc[i]['Files']:
                reads = [f"{snakemake.params.output}/Preprocess/Trimmomatic/quality_trimmed_{pexps.iloc[i]['Name']}_{fr}_paired.fq"
                         for fr in ['forward', 'reverse']]
            else:
                reads = [f"{snakemake.params.output}/Preprocess/Trimmomatic/quality_trimmed_{pexps.iloc[i]['Name']}.fq"]
            perform_alignment(
                reference, reads, f"{snakemake.params.output}/Quantification/{pexps.iloc[i]['Name']}",
                threads=snakemake.threads)
            normalize_counts_by_size(
                f"{snakemake.params.output}/Quantification/{pexps.iloc[i]['Name']}.readcounts", reference)
            # Read the results of alignment and add them to the readcounts result at sample level
            counts = pd.read_csv(
                f"{snakemake.params.output}/Quantification/{pexps.iloc[i]['Name']}_normalized.readcounts", sep='\t',
                names=['Gene' if pexps.iloc[i]['Data type'] == 'mrna' else 'Contig', pexps.iloc[i]['Name']])
            if pexps.iloc[i]['Data type'] == 'dna':
                mg_result = pd.merge(mg_result, counts, how='outer', on='Contig')
            else:
                mt_result = pd.merge(mt_result, counts, how='outer', on='Gene')
        Path(f"{snakemake.params.output}/Quantification/{sample}").mkdir(parents=True, exist_ok=True)
        if len(mg_result) > 0:
            mg_result.to_csv(f"{snakemake.params.output}/Quantification/{sample}/mg.readcounts", sep='\t', index=False)
        if len(mt_result) > 0:
            mt_result.astype(int, errors='ignore').to_csv(
                f"{snakemake.params.output}/Quantification/{sample}/mt.readcounts", sep='\t', index=False)


if __name__ == '__main__':
    run()
