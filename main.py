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

import argparse, pathlib, os, time

mtools = MoscaTools()

parser = argparse.ArgumentParser(description="Multi Omics Software for Community Analysis",
                                 epilog="""A tool for performing metagenomics, metatranscriptomics 
                                 and metaproteomics analysis.""")
parser.add_argument("-f","--files", type=str, 
                    help="Input files for analysis (mg1R1,mg1R2:mt1R1,mt1R2;)", nargs = '*')
parser.add_argument("-data","--type-of-data",default='paired',choices=["paired","single"],
                    action='store',type=str,help='Type of data (paired/single)-ended')
parser.add_argument("-a","--assembler",type=str,choices=["metaspades","megahit"],
                    help="Tool for assembling the reads", nargs = 1,metavar = "Assembler", 
                    default = "metaspades")
parser.add_argument("-db","--annotation-database",type=str,help="Database for annotation (.fasta or .dmnd)", 
                    nargs = 1,metavar = "Database", default = "Databases/Annotation/uniprot.dmnd")
parser.add_argument("-o","--output",type=str,help="Directory for storing the results",
                    metavar = "Directory")
parser.add_argument("-nopp","--no-preprocessing",action = "store_true",
                    help="Don't perform preprocessing", default = False)
parser.add_argument("-noas","--no-assembly",action = "store_true",
                    help="Don't perform assembly", default = False)
parser.add_argument("-noan","--no-annotation",action = "store_true",help="Don't perform annotation", 
                    default = False)
parser.add_argument("-nobin","--no-binning",action = "store_true",help="Don't perform binning", default = False)
parser.add_argument("-ol","--output-level",default='min',choices=["min","med","max"],
                    action='store',type=str,help=("""Level of file output from MOSCA, 
                    min outputs only the analysis results, med removes intermediate files, max outputs 
                    all intermediate and final data"""))
parser.add_argument("-mp", "--metaproteomic", action = "store_true",
                    help = ("""If data is metagenomic and metaproteomic, if not specified
                            will be assumed to be metagenomic and metatranscriptomic"""),
                    default = False)
parser.add_argument("-c","--conditions", type=str, help="Different conditions for metatranscriptomic analysis", 
                    nargs = '*')

             
args = mtools.validate_arguments(parser)

mtools.print_arguments(args)            # TODO - refine this function

experiments = [experiment for experiment in args.files]

print('Creating directories at ' + args.output)
directories = ([args.output + '/Preprocess/' + software for software in 
                ['FastQC', 'Trimmomatic', 'SortMeRNA']] + 
               [args.output + folder for folder in ['/Assembly','/Annotation',
                '/Metatranscriptomics_analysis' if not args.metaproteomic else
                '/Metaproteomics_analysis']])

for directory in directories:
    print('Created ' + directory)
    path = pathlib.Path(directory)
    path.mkdir(parents=True, exist_ok=True)

already_processed = list()
expression_names = list()

for experiment in experiments:
    pairs = experiment.split(':')
    mg = pairs[0].split(',')
    
    if len(mg) == 1 and args.type_of_data == 'paired':                           # if data is interleaved paired end, it will be split up
        (forward, reverse) = (mg[0].split('/')[-1].split('.fastq')[0] + fr for 
                             fr in ['_R1.fastq','_R2.fastq'])
        #mtools.divide_fq(mg[0], forward, reverse)
        mg = [forward, reverse]
    
    mg_name = mg[0].split('/')[-1].split('_R')[0]
    
    print(mg_name)
    
    if mg_name not in already_processed:                                        #several MT samples might correspond to the same MG sample
        '''    
        Metagenomics Preprocess        
        '''
        if not args.no_preprocessing:
            print('Preprocessing of metagenomic reads has begun')
            preprocesser = Preprocesser(files = mg,
                                         paired = 'PE',
                                         working_dir = args.output,
                                         data = 'dna',
                                         name = mg_name)
            #preprocesser.run()
        
        '''
        Assembly
        '''
        
        if not args.no_assembly:
            pathlib.Path(args.output + '/Assembly/' + mg_name).mkdir(parents=True, exist_ok=True)
            
            print('Assembly has begun')
            assembler = Assembler(out_dir = args.output,
                                 assembler = args.assembler,
                                 name = mg_name,
                                 forward = args.output + '/Preprocess/Trimmomatic/quality_trimmed_' + mg_name + '_forward_paired.fq',
                                 reverse = args.output + '/Preprocess/Trimmomatic/quality_trimmed_' + mg_name + '_reverse_paired.fq')
            #assembler.run()
        
        '''
        Annotation
        '''
        
        if not args.no_annotation:
            pathlib.Path(args.output + '/Annotation/' + mg_name).mkdir(parents=True, exist_ok=True)
            assembled = False if args.no_assembly else True
            
            print('Annotation has begun')
            annotater = Annotater(out_dir = args.output,
                                 assembler = args.assembler,
                                 db = args.annotation_database,
                                 assembled = assembled,
                                 error_model = 'illumina_10',
                                 name = mg_name)
            #annotater.run()
            
        '''
        Binning
        '''
        if not args.no_binning:
            print('Binning has begun')
            pathlib.Path(args.output + '/Binning/' + mg_name).mkdir(parents=True, exist_ok=True)
            
            binner = Binner(output = args.output + '/Binning/' + mg_name,
                            contigs = args.output + '/Assembly/' + mg_name + '/contigs.fasta',
                            blast = args.output + '/Annotation/' + mg_name + '/aligned.blast',
                            uniprotinfo = args.output + '/Annotation/' + mg_name + '/uniprot.info')
            binner.run()
            
        already_processed.append(mg)
    
    if len(pairs) > 1:
        if not args.metaproteomic:
            ''''
            Metatranscriptomics Preprocessing
            '''
            mt = pairs[1].split(',') 
            mt_name = mt[0].split('/')[-1].split('_R')[0]
            print(mt_name)
                
            print('Preprocessing of metatranscriptomic reads has begun')
            preprocesser = Preprocesser(files = mt,
                                        paired = 'PE',
                                        working_dir = args.output,
                                        trimmomatic_dir = os.path.expanduser('~/anaconda3/bin'),
                                        data = 'mrna',
                                        name = mt_name)
            preprocesser.run()
            
            '''
            Metatranscriptomics Quantification
            '''
            print('Metatranscriptomics Analysis has begun')
            
            path = pathlib.Path(args.output + '/MetaTranscriptomics/' + mt_name).mkdir(parents=True, exist_ok=True)        
            
            mta = MetaTranscriptomicsAnalyser(out_dir = args.output + '/MetaTranscriptomics',
                                              contigs = args.output + '/Assembly/' + mg_name + '/contigs.fasta',
                                              blast = args.output + '/Annotation/' + mg_name + '/aligned.blast',
                                              mg = mg_name,
                                              mt = mt_name,
                                              assembler = args.assembler)
            mta.readcounts_file()
            
            expression_names.append(mt_name)
            
        else:
            '''
            MetaProteomics Analysis
            '''
            print('MetaProteomics Analysis has begun')
            
            path = pathlib.Path(args.output + '/MetaProteomics/' + mt_name).mkdir(parents=True, exist_ok=True)   
            
            spectra_folder = pairs[1]
            analyser = MetaProteomicsAnalyser(faa = args.output + '/Annotation/' + mg_name + '/fgs.faa',
                          blast = args.output + '/Annotation/' + mg_name + '/aligned.blast',
                          crap_folder = '/HDDStorage/jsequeira/Thesis/MetaProteomics',
                          output = args.output + '/MetaProteomics/' + mg_name,
                          protease = 'trypsin',
                          spectra_folder = spectra_folder,
                          experiment_name = args.output,
                          sample_name = mg_name,
                          replicate_number = '1')
        
            analyser.standard_run()

if len(experiment.split(':')) > 1:      
    if not args.metaproteomic:
        mta.merge_readcounts(args.output + '/MetaTranscriptomics', expression_names, 
                             args.output + '/MetaTranscriptomics/all_experiments.readcounts')

        mta.differential_analysis(args.output + '/MetaTranscriptomics/all_experiments.readcounts',
                                  args.output + '/MetaTranscriptomics/', args.conditions)