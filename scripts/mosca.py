# -*- coding: utf-8 -*-
"""
MOSCA main class

By João Sequeira

Sep 2017
"""

from mosca_tools import MoscaTools
from preprocess import Preprocesser
from assembly import Assembler
from annotation import Annotater
from binning import Binner
from metatranscriptomics_analyser import MetaTranscriptomicsAnalyser
from metaproteomics_analyser import MetaProteomicsAnalyser
from kegg_pathway import KEGGPathway
from time import gmtime, strftime

import argparse, pathlib, os, glob, multiprocessing

mosca_dir = os.path.dirname(os.path.realpath(__file__))

mtools = MoscaTools()

parser = argparse.ArgumentParser(description="Multi Omics Software for Community Analysis",
                                 epilog="""A tool for performing metagenomics, metatranscriptomics 
                                 and metaproteomics analysis.""")
parser.add_argument("-f","--files", type=str, nargs = '*', required=True,
                    help="Input files for analysis (mg1R1,mg1R2:mt1R1,mt1R2 mg2R1,...)")
parser.add_argument("-st","--sequencing-technology",default='paired',choices=["paired","single"],
                    action='store',type=str,help='Type of data (paired/single)-ended')
parser.add_argument("-a","--assembler",type=str,choices=["metaspades","megahit"],
                    help="Tool for assembling the reads",
                    default = "metaspades")
parser.add_argument("-db","--annotation-database",type=str,
                    help="Database for annotation (.fasta or .dmnd)", 
                     default = "Databases/annotation_databases/uniprot.fasta")
parser.add_argument("-o","--output",type=str,help="Directory for storing the results",
                    metavar = "Directory", default = "/MOSCA_analysis")
parser.add_argument("-nopp","--no-preprocessing",action = "store_true",
                    help="Don't perform preprocessing", default = False)
parser.add_argument("-noas","--no-assembly",action = "store_true",
                    help="Don't perform assembly", default = False)
parser.add_argument("-noan","--no-annotation",action = "store_true",default = False,
                    help="Don't perform annotation")
parser.add_argument("-nobin","--no-binning",action = "store_true",default = False,
                    help="Don't perform binning")
parser.add_argument("-ol","--output-level", default = 'maximum', type = str,
                    choices = ["minimum","medium","maximum"], action = 'store',
                    help=("""Level of file output from MOSCA, minimum outputs 
                          only the analysis results, medium removes intermediate 
                          files, maximum outputs all intermediate and final data"""))
parser.add_argument("-tod", "--type-of-data", default = "metatranscriptomics",
                    help = ("""If data is metagenomics integrated with metatranscriptomics
                            or metaproteomics, if not specified will be assumed to be metagenomics
                            and metatranscriptomics. This option can be ignored when dealing
                            only with metagenomics data"""), 
                            choices=["metatranscriptomics", "metaproteomics"])
parser.add_argument("-c","--conditions", type=str, nargs = '*',
                    help="""Different conditions for metatranscriptomics/metaproteomics 
                    analysis, separated by comma (,)""")
parser.add_argument("-t","--threads",type=str, metavar = "Threads", 
                    default = str(multiprocessing.cpu_count()),
                    help="Number of threads available for MOSCA. Default is number of cores available.")
parser.add_argument("-m","--memory",type=str,
                    help="Maximum memory (byte) available for MOSCA. Applied only in the assembly")
parser.add_argument("-mark","--marker-gene-set",type=str,
                    help="""Marker gene set to use for binning with MaxBin2. 
                    107 if archaea are not to be considered, 40 if data is diverse.""",
                    default = '40')
parser.add_argument("-tr","--fgs-train-file",type=str,
                    help="""File name that contains model parameters of sequencing,
                    related to sequencing technology and error rate""", default = 'illumina_5', 
                    choices=["sanger_5","sanger_10","454_10","454_30","illumina_5",
                             "illumina_10"])
parser.add_argument("-mt", "--metaquast-threshold", type = str, help="""Minimum 
                    contig length to be considered in assembly quality control""",
                    default = '500')
parser.add_argument("-bp", "--bowtie2-profile", type = str, 
                    help="""Profile for bowtie2 alignment of reads""", 
                    choices = ["very-fast", "fast", "sensitive", "very-sensitive",
                               "very-fast-local", "fast-local", "sensitive-local", 
                               "very-sensitive-local"])
parser.add_argument("-kl", "--kmer-list", type = str, 
                    help="""Kmer sizes for multi-kmer assembly""")
parser.add_argument("-qs", "--quality-score", type = str, 
                    help="""Scoring system of quality of reads""")
parser.add_argument("-bcl", "--binning-contig-legth", type = str, 
                    help="""Minimum length of contigs to be considered for binning""")
parser.add_argument("-anncols", "--annotation-columns", type = str, 
                    help="""List of UniProt columns to obtain information from""")
parser.add_argument("-anndbs", "--annotation-databases", type = str, 
                    help="""List of databases to cross-check with UniProt information""")
             
args = mtools.validate_arguments(parser)

mtools.print_arguments(args)            # TODO - refine this function
monitorization_file = args.output + '/monitorization_report.txt'

experiments = [experiment for experiment in args.files]

mtools.timed_message('MOSCA analysis has started')

print('Creating directories at ' + args.output)
directories = ([args.output + '/Preprocess/' + software for software in 
                ['FastQC', 'Trimmomatic', 'SortMeRNA']] + 
               [args.output + folder for folder in ['/Assembly','/Annotation',
                '/Metatranscriptomics' if args.type_of_data == 'metatranscriptomics' 
                else '/Metaproteomics']])

for directory in directories:
    print('Created ' + directory)
    path = pathlib.Path(directory)
    path.mkdir(parents=True, exist_ok=True)

mg_processed = list()
expression_analysed = list()

for experiment in experiments:
    pairs = experiment.split(':')
    mg = pairs[0].split(',')
    
    if len(mg) == 1 and args.sequencing_technology == 'paired':                 # if data is interleaved paired end, it will be split up to forward and reverse files
        (forward, reverse) = (args.output + '/Preprocess/' + mg[0].split('/')[-1].split('.fastq')[0]
                                + fr for fr in ['_R1.fastq','_R2.fastq'])
        mtools.timed_message('Splitting reads from {} to {} and {}'.format(
                mg[0], forward, reverse))
        mtools.divide_fq(mg[0], forward, reverse)
        mg = [forward, reverse]
            
    mg_name = mg[0].split('/')[-1].split('_R')[0]
    
    if mg_name not in mg_processed:                                             #several MT samples might correspond to the same MG sample

        '''    
        Metagenomics Preprocess        
        '''
        if not args.no_preprocessing:
            mtools.timed_message('Preprocessing metagenomic reads')
            preprocesser = Preprocesser(files = mg,
                                        paired = 'PE' if args.sequencing_technology == 'paired' else 'SE',
                                        working_dir = args.output,
                                        data = 'dna',
                                        name = mg_name,
                                        threads = args.threads)
            if hasattr(args, 'quality_score'):
                setattr(preprocesser, 'quality_score', args.quality_score)
                
            preprocesser.run()
            
            mtools.task_is_finished(task = 'Preprocessing',
                    file = monitorization_file, 
                    task_output = args.output + '/Preprocess')
            
            mg = [args.output + '/Preprocess/Trimmomatic/quality_trimmed_' + mg_name + 
                  '_' + fr + '_paired.fq' for fr in ['forward', 'reverse']]
        
        '''
        Assembly
        '''
        
        if not args.no_assembly:
            pathlib.Path(args.output + '/Assembly/' + mg_name).mkdir(parents=True, exist_ok=True)
            
            mtools.timed_message('Assembling reads')
            
            mid_name = args.output + '/Preprocess/Trimmomatic/quality_trimmed_'
            
            assembler = Assembler(out_dir = args.output,
                                 assembler = args.assembler,
                                 name = mg_name,
                                 forward = mg[0],
                                 reverse = mg[1],
                                 threads = args.threads)
            
            if (args.assembler == 'metaspades' and hasattr(args, 'quality_score')
            and getattr(args, 'quality_score') not in ['None', None]):          # Megahit doesn't accept quality score input
                setattr(assembler, 'phred_offset', args.quality_score)              # --phred-offset is the name of the parameter in MetaSPAdes
            if args.memory not in [None, 'None']:
                setattr(assembler, 'memory', args.memory)
            
            assembler.run()
            
            mtools.task_is_finished(task = 'Assembly',
                    file = monitorization_file, 
                    task_output = args.output + '/Assembly/' + mg_name)
        
        '''
        Annotation
        '''
        if not args.no_annotation:
            pathlib.Path(args.output + '/Annotation/' + mg_name).mkdir(parents=True, exist_ok=True)
            assembled = False if args.no_assembly else True
            
            mtools.timed_message('Annotating sequences')
            annotater = Annotater(file = args.output + '/Assembly/' + mg_name + '/contigs.fasta',
                                  out_dir = args.output,
                                  assembler = args.assembler,
                                  db = args.annotation_database,
                                  assembled = assembled,
                                  error_model = 'illumina_10',
                                  name = mg_name,
                                  threads = args.threads,
                                  columns = args.annotation_columns.split(','),
                                  databases = args.annotation_databases.split(','))
            
            annotater.run()
            
            mtools.task_is_finished(task = 'Annotation',
                    file = monitorization_file, 
                    task_output = args.output + '/Annotation/' + mg_name)
            
        '''
        Binning
        '''
        if not args.no_binning:
            mtools.timed_message('Binning has begun')
            
            pathlib.Path(args.output + '/Binning/' + mg_name).mkdir(parents=True, exist_ok=True)
            
            binner = Binner(output = args.output + '/Binning/' + mg_name,
                            contigs = args.output + '/Assembly/' + mg_name + '/contigs.fasta',
                            blast = args.output + '/Annotation/' + mg_name + '/aligned.blast',
                            uniprotinfo = args.output + '/Annotation/' + mg_name + '/uniprot_info.tsv',
                            threads = args.threads,
                            mg1 = args.output + '/Preprocess/Trimmomatic/quality_trimmed_' + mg_name + '_forward_paired.fq',
                            mg2 = args.output + '/Preprocess/Trimmomatic/quality_trimmed_' + mg_name + '_reverse_paired.fq',
                            markerset = args.marker_gene_set)
            binner.maxbin_workflow()

            mtools.task_is_finished(task = 'Binning',
                    file = monitorization_file, 
                    task_output = args.output + '/Binning/' + mg_name)
            
        mg_processed.append(mg_name)

    if len(pairs) > 1:
        
        if args.type_of_data == 'metatranscriptomics':
            
            ''''
            Metatranscriptomics Preprocessing
            '''
            
            mt = pairs[1].split(',')
            
            if len(mt) == 1 and args.sequencing_technology == 'paired':                           # if data is interleaved paired end, it will be split up
            
                mt_name = mt[0].split('/')[-1].split('.fastq')[0]
                
                (forward, reverse) = (args.output + '/Preprocess/' + mt[0].split('/')[-1].split('.fastq')[0]
                                        + fr for fr in ['_R1.fastq','_R2.fastq'])
                print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Splitting ' + 
                      'reads at ' + mt[0] + ' to ' + forward + ' and ' + reverse + '.')
                mtools.divide_fq(mt[0], forward, reverse)
                mt = [forward, reverse]
                
            else:
                mt_name = mt[0].split('/')[-1].split('_R')[0]
                
            mtools.timed_message('Analysis ' + mt_name + ' sample data')
                
            preprocesser = Preprocesser(files = mt,
                                        paired = 'PE',
                                        working_dir = args.output,
                                        data = 'mrna',
                                        name = mt_name,
                                        threads = args.threads)
            if hasattr(args, 'quality_score'):
                setattr(preprocesser, 'quality_score', args.quality_score)
                
            preprocesser.run()
            
            '''
            Metatranscriptomics Quantification
            '''
            mtools.timed_message('Analysing gene expression with metatranscriptomics')
            
            path = pathlib.Path(args.output + '/Metatranscriptomics/' + 
                                mg_name).mkdir(parents=True, exist_ok=True)
            
            mt = [args.output + '/Preprocess/Trimmomatic/quality_trimmed_' + mt_name + 
                  '_' + fr + '_paired.fq' for fr in ['forward', 'reverse']]
            
            mta = MetaTranscriptomicsAnalyser(out_dir = args.output + '/Metatranscriptomics',
                          contigs = args.output + '/Assembly/' + mg_name + '/contigs.fasta',
                          blast = args.output + '/Annotation/' + mg_name + '/aligned.blast',
                          reads = mt,
                          mt = mt_name,
                          threads = args.threads)
            mta.readcounts_file()
            
            mtools.task_is_finished(task = 'Metatranscriptomics analysis',
                  file = monitorization_file, 
                  task_output = args.output + '/Metatranscriptomics/' + mt_name)
            
            expression_analysed.append(mt_name)
            
        else:
            
            '''
            MetaProteomics Analysis
            '''
            mtools.timed_message('Analysing gene expression with metaproteomics')
            
            path = pathlib.Path(args.output + '/Metaproteomics/' + mt_name).mkdir(parents=True, exist_ok=True)   
            
            spectra_folder = pairs[1]
            analyser = MetaProteomicsAnalyser(faa = args.output + '/Annotation/' + mg_name + '/fgs.faa',
                          blast = args.output + '/Annotation/' + mg_name + '/aligned.blast',
                          crap_folder = '/HDDStorage/jsequeira/Thesis/Metaproteomics',
                          output = args.output + '/Metaproteomics/' + mg_name,
                          protease = 'trypsin',
                          spectra_folder = spectra_folder,
                          experiment_name = args.output,
                          sample_name = mg_name,
                          replicate_number = '1',
                          threads = args.threads)
        
            analyser.standard_run()
            
            mtools.task_is_finished(task = 'Metaproteomics analysis',
                    file = monitorization_file, 
                    task_output = args.output + '/Metaproteomics/' + mt_name)
            
'''
Join all information on one report
'''
mtools.timed_message('Integrating all information.')

annotater = Annotater(out_dir = args.output,
                      threads = args.threads)

joined = annotater.global_information()

mtools.timed_message('Integration is available at ' + args.output)

if len(experiment.split(':')) > 1:      
    if args.type_of_data == 'metatranscriptomics':
        readcount_files = glob.glob(args.output + '/Metatranscriptomics/*.readcounts')
        
        for file in readcount_files[:-1]:
            mt_name = file.split('/')[-1].split('.')[0]
            joined = mtools.define_abundance(joined, origin_of_data = 'metatranscriptomics',name = mt_name, readcounts = file)
        
        mta.merge_readcounts(readcount_files, expression_analysed, args.output + '/Metatranscriptomics/all_experiments.readcounts')

        mta.differential_analysis(args.output + '/Metatranscriptomics/all_experiments.readcounts', args.conditions[0].split(','), args.output + '/Metatranscriptomics/')
        
    else:
        # TODO - the metaproteomics workflow
        pass
    
samples = mg_processed + expression_analysed

joined[samples].to_csv(args.output + '/readcounts.table',
      sep = '\t', index = False)

joined = mtools.normalize_readcounts(args.output + '/readcounts.table', samples,
                            args.output + '/normalization_factors.txt')

joined.to_csv(args.output + '/mosca_results.tsv', sep = '\t', index = False)
joined.to_excel(args.output + '/mosca_results.xlsx', index = False)

# KEGG Pathway representations
kp = KEGGPathway(input_file = args.output + '/mosca_results.tsv',
                 output_directory = args.output + '/Metatranscriptomics/KEGG Pathway',
                 mg_samples = mg_processed,
                 mt_samples = expression_analysed)
kp.run()

mtools.timed_message('Analysis with MOSCA was concluded with success!')

mtools.write_technical_report(args.output + 'technical_report.txt')             # TODO - must be improved
print('Versions of the tools used is available at {}/technical_report.txt'.format(args.output))