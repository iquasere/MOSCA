# -*- coding: utf-8 -*-
"""
MOSCA main class

By João Sequeira

Sep 2017
"""

from preprocessing import Preprocessing
from assembling import Assembling
from annotating import Annotating
from analysing import Analysing
from mosca_tools import MoscaTools
from diamond import DIAMOND

import argparse, pathlib

parser = argparse.ArgumentParser(description="Multi Omics Software for Community Analysis",
                                 epilog="A tool for performing metagenomics and metatranscriptomics analysis.")
parser.add_argument("-mg","--metagenomic", type=str, help="Input files for metagenomic analysis", nargs = '*',
                    metavar = "Metagenomic files")
parser.add_argument("-mt","--metatranscriptomic", type=str, help="Input files for metatranscriptomic analysis", nargs = '*',
                    metavar = "Metatranscriptomic files")
parser.add_argument("-data","--type-of-data",default='paired',choices=["paired","single"],action='store',type=str,
                    help='Type of data (paired/single)-ended', metavar = "[paired/single]")
parser.add_argument("-a","--assembler",type=str,choices=["metaspades","megahit"],help="Tool for assembling the reads", nargs = 1,
                    metavar = "Assembler", default = "metaspades")
parser.add_argument("-db","--annotation-database",type=str,help="Database for annotation (.fasta or .dmnd)", nargs = 1,
                    metavar = "Database", default = "Databases/Annotation/UniProt/uniprot.dmnd")
parser.add_argument("-o","--output-dir",type=str,help="Directory for storing the results",metavar = "Directory")
parser.add_argument("-nopp","--no-preprocessing",action = "store_true",help="Don't perform preprocessing", default = False)
parser.add_argument("-noas","--no-assembly",action = "store_true",help="Don't perform assembly", default = False)
parser.add_argument("-noan","--no-annotation",action = "store_true",help="Don't perform annotation", default = False)
parser.add_argument("-node","--no-differential-expression",action = "store_true",help="Don't perform differential expression analysis",
                    default = False)
parser.add_argument("-of","--output-files",default='min',choices=["min","med","max"],action='store',type=str,
                    help=("""Level of file output from MOSCA, min outputs only the analysis results, med removes intermediate files, max outputs 
                    all intermediate and final data"""))
parser.add_argument("-mgfq", "--mg-fastq", default = '0', choices = ["0","1"], action = 'store', type = str,
                    help='If metagenomic data is either in fastq (0) or fasta (1) format')
parser.add_argument("-mtfq","--mt-fastq",default='0',choices=["0","1"],action='store',type=str,
                    help='If metatranscriptomic data is either in fastq (0) or fasta (1) format')
             
args = parser.parse_args()

nice_arguments = True
if args.metagenomic is None:
    print()
    print('Must specify which metagenomic files to use!')
    nice_arguments = False
if args.output_dir is None:
    print()
    print('Must specify which output directory to use!')
    nice_arguments = False
if nice_arguments == False:
    print()
    parser.print_help()
    exit()
    
print('Creating directories at ' + args.output_dir)
directories = ([args.output_dir + '/Preprocess/' + software for software in ['FastQC', 'Trimmomatic', 'SortMeRNA']] + 
               [args.output_dir + folder for folder in ['/Assembly/Quality_control','/Annotation','/Analysis']])

for directory in directories:
    print('Created ' + directory)
    path = pathlib.Path(directory)
    path.mkdir(parents=True, exist_ok=True)
    
mg = args.metagenomic[0].split(',')
if args.metatranscriptomic is not None:
    mt = args.metatranscriptomic[0].split(',')
else:
    mt = list()

#Preprocess
if not args.no_preprocessing:
    interleaved = False
    if args.type_of_data == 'paired' and len(mg) == 1:
        interleaved = True
        print('Supplied interleaved paired-end file being splited into forward and reverse files')
        MoscaTools.divide_fq(mg[0], args.output_dir + '/mg1.fastq', args.output_dir + '/mg2.fastq')
        mg = [directory + '/mg1.fastq', args.output_dir + '/mg2.fastq']
        
    print('Preprocessing of metagenomic reads has began')
    preprocesser = Preprocessing(files = mg,
                                 paired = 'PE',
                                 working_dir = args.output_dir,
                                 trimmomatic_dir = '~/anaconda3/jar/',
                                 data = 'dna')
    #preprocesser.run()

    if len(mt) > 0:
        if args.mt_fastq:
            print('Preprocessing of metatranscriptomic reads has began')
            preprocesser = Preprocessing(files = mt,
                                         paired = 'PE',
                                         working_dir = args.output_dir,
                                         trimmomatic_dir = '~/anaconda3/jar/',
                                         data = 'mrna')
            #preprocesser.run()

#Assembly
if not args.no_assembly:
    print('Assembly has began')
    assemble = Assembling(out_dir = args.output_dir,
                           assembler = args.assembler,
                           forward_paired = args.output_dir + '/Preprocess/SortMeRNA/forward_rejected.fastq',
                           reverse_paired = args.output_dir + '/Preprocess/SortMeRNA/reverse_rejected.fastq')
    assemble.run()

#Annotation
if not args.no_annotation:
    assembled = False if args.no_assembly else True
    print('Annotation has began')
    annotate = Annotating(out_dir = args.output_dir,
                          assembler = args.assembler,
                          db = args.annotation_database,
                          assembled = assembled,
                          error_model = 'illumina_10')
    ids = annotate.run()
    
#Data Analysis
if not args.no_differential_expression:
    print('Analysis has began')
    analyser = Analysing(out_dir = args.output_dir,
                         assembler = args.assembler,
                         ids = ids,
                         mt = mt,
                         de = True)
    analyser.full_analysis()

