"""
MOSCA's Assembly package for performing Assembly with MetaSPAdes
and Megahit and Quality Control analysis of the resulting contigs

By João Sequeira

Jun 2017
"""

import os
import pathlib
import shutil
from mosca_tools import run_command, run_pipe_command, perform_alignment


class Assembler:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def run_assembler_mg(self, reads, out_dir, assembler, threads='12', memory=None):
        run_command(
            f"{assembler if assembler == 'megahit' else f'{assembler}.py'} -o {out_dir} -t {threads} "
            f"""{f'-1 {reads[0]} -2 {reads[1]}' if len(reads) == 2 else 
            f'-{"r" if assembler == "megahit" else "s"} {reads[0]}'} """
            f"-m {round(memory) if memory else ''}")

    def run_assembler_mt(self, reads, out_dir, threads, memory=None):
        reads = f'--left {reads[0]} --right {reads[1]}' if len(reads) == 2 else f'--single {reads[0]}'
        run_command(
            f"Trinity --seqType fq {reads} --CPU {threads} --max_memory {int(memory)}G --output {out_dir}/trinity")

    def percentage_of_reads(self, file):
        handler = open(file)
        lines = handler.readlines()
        return lines[-1].split('%')[0]

    def run_metaquast(self, contigs, out_dir, threads='12', max_ref_number=0):
        run_command(
            f'metaquast.py --threads {threads} --output-dir {out_dir} --max-ref-number {max_ref_number} {contigs}')

    def close_gaps(
            self, scaffolds, contigs, output_basename, read1, read2, read_length=150, insert_size=1500, std_insert=10,
            threads=14, base_qual='phred33'):
        run_command(
            f'gmcloser --target_scaf {scaffolds} --query_seq {contigs} --prefix_out {output_basename} '
            f'-r {read1} {read2} --read_len {read_length} --insert {insert_size} --sd_insert {std_insert} '
            f'--thread {threads} --base_qual {base_qual}')

    def run(self):
        # Assembly
        if snakemake.params.assembler == 'megahit':     # snakemake workflow creates output directory, and megahit don't like it
            if os.path.isdir(snakemake.params.output):
                shutil.rmtree(snakemake.params.output)

        if snakemake.params.assembler != 'trinity':
            self.run_assembler_mg(
                snakemake.params.reads, snakemake.params.output, snakemake.params.assembler, threads=snakemake.threads,
                memory=snakemake.params.max_memory)
        else:
            self.run_assembler_mt(
                snakemake.params.reads, snakemake.params.output, threads=snakemake.threads,
                memory=snakemake.params.max_memory)
        if snakemake.params.assembler == 'megahit':  # all contigs files are outputed the metaspades way
            run_pipe_command(f"awk \'{{print $1}}\' {snakemake.params.output}/final.contigs.fa",
                             # k141_714 flag=1 multi=1.0000 len=369 -> k141_714
                             output=f'{snakemake.params.output}/contigs.fasta')
            shutil.copyfile(f'{snakemake.params.output}/contigs.fasta', f'{snakemake.params.output}/scaffolds.fasta')   # TODO - put SOAPdenovo producing scaffolds from Megahit
        elif snakemake.params.assembler == 'trinity':
            run_pipe_command(f"awk \'{{print $1}}\' {snakemake.params.output}/trinity/Trinity.fasta",
                             # >TRINITY_DN19419_c0_g1_i1 len=248 path=[0:0-247] -> TRINITY_DN19419_c0_g1_i1
                             output=f'{snakemake.params.output}/contigs.fasta')
            shutil.copyfile(f'{snakemake.params.output}/contigs.fasta', f'{snakemake.params.output}/scaffolds.fasta')

        # TODO - reactivate gap closing
        '''
        self.close_gaps(
            f'{snakemake.params.output}/scaffolds.fasta', f'{snakemake.params.output}/contigs.fasta',
            f'{snakemake.params.output}/gap_close', snakemake.params.reads[0], snakemake.params.reads[1],
            read_length=snakemake.params.read_length, insert_size=snakemake.params.insert_size,
            std_insert=snakemake.params.std_insert, threads=snakemake.params.threads,
            base_qual=snakemake.params.base_qual)
        '''

        # Quality control
        pathlib.Path(f'{snakemake.params.output}/quality_control').mkdir(parents=True, exist_ok=True)

        self.run_metaquast(f'{snakemake.params.output}/contigs.fasta', f'{snakemake.params.output}/quality_control')
        perform_alignment(
            f'{snakemake.params.output}/contigs.fasta', snakemake.params.reads,
            f'{snakemake.params.output}/quality_control/alignment', threads=snakemake.threads)
        percentage_of_reads = self.percentage_of_reads(f'{snakemake.params.output}/quality_control/alignment.log')

        if os.path.isfile(f'{snakemake.params.output}/quality_control/combined_reference/report.tsv'):  # if metaquast finds references to contigs, it will output results to different folders
            shutil.copyfile(f'{snakemake.params.output}/quality_control/combined_reference/report.tsv',
                            f'{snakemake.params.output}/quality_control/report.tsv')
        with open(snakemake.params.output + '/quality_control/report.tsv', 'a') as f:
            f.write(f'Reads aligned (%)\t{percentage_of_reads}\n')


if __name__ == '__main__':
    Assembler().run()
