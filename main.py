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
from binner import Binner
from metatranscriptomics_analyser import MetaTranscriptomicsAnalyser
from metaproteomics_analyser import MetaProteomicsAnalyser

import argparse, pathlib, os, time

mtools = MoscaTools()

parser = argparse.ArgumentParser(description="Multi Omics Software for Community Analysis",
                                 epilog="""A tool for performing metagenomics, metatranscriptomics 
                                 and metaproteomics analysis.""")
parser.add_argument("-f","--files", type=str, 
                    help="Input files for analysis (mg1R1,mg1R2:mt1R1,mt1R2;)", nargs = '*',
                    metavar = "Metagenomic files")
parser.add_argument("-data","--type-of-data",default='paired',choices=["paired","single"],
                    action='store',type=str,help='Type of data (paired/single)-ended', 
                    metavar = "[paired/single]")
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
parser.add_argument("-nobin","--no-binning",action = "store_true",
                    help="Don't perform binning", default = False)
parser.add_argument("-ol","--output-level",default='min',choices=["min","med","max"],
                    action='store',type=str,help=("""Level of file output from MOSCA, 
                                                  min outputs only the analysis results, 
                                                  med removes intermediate files, max outputs 
                                                  all intermediate and final data"""))
parser.add_argument("-mp", "--metaproteomic", action = "store_true",
                    help = ("""If data is metagenomic and metaproteomic, if not specified
                            will be assumed to be metagenomic and metatranscriptomic"""),
                    default = False)
             
args = mtools.validate_arguments(parser)

mtools.print_arguments(args)            # TODO - refine this function

experiments = [experiment for experiment in args.files]

print('Creating directories at ' + args.output)
directories = ([args.output + '/Preprocess/' + software for software in 
                ['FastQC', 'Trimmomatic', 'SortMeRNA']] + 
               [args.output + folder for folder in ['/Assembly','/Annotation']] +
               ['/Metatranscriptomics_analysis' if not parser.metaproteomic else
                '/Metaproteomics_analysis'])

for directory in directories:
    print('Created ' + directory)
    path = pathlib.Path(directory)
    path.mkdir(parents=True, exist_ok=True)

already_processed = list()

for experiment in experiments:
    pairs = experiment.split(':')
    mg = pairs[0].split(',')
    mg_name = mg[0].split('/')[-1].split('_R')[0]
    
    print(mg_name)
    if mg_name not in already_processed:                                        #several MT samples might correspond to the same MG sample
        '''    
        Metagenomics Preprocess        
        '''
        if not args.no_preprocessing:
            interleaved = False
            if args.type_of_data == 'paired' and len(mg) == 1:
                interleaved = True
                print('Supplied interleaved paired-end file being splited into forward and reverse files')
                MoscaTools.divide_fq(mg[0], args.output + '/mg1.fastq', args.output + '/mg2.fastq')
                mg = [directory + '/mg1.fastq', args.output + '/mg2.fastq']
                
            print('Preprocessing of metagenomic reads has begun')
            preprocesser = Preprocesser(files = mg,
                                         paired = 'PE',
                                         working_dir = args.output,
                                         data = 'dna',
                                         name = mg_name)
            preprocesser.run()
        
        '''
        Assembly
        '''
        pathlib.Path(args.output + '/Assembly/' + mg_name).mkdir(parents=True, exist_ok=True)
        if not args.no_assembly:
            print('Assembly has begun')
            assembler = Assembler(out_dir = args.output,
                                 assembler = args.assembler,
                                 name = mg_name,
                                 forward = args.output + '/Preprocess/Trimmomatic/quality_trimmed_' + mg_name + '_forward_paired.fq',
                                 reverse = args.output + '/Preprocess/Trimmomatic/quality_trimmed_' + mg_name + '_reverse_paired.fq')
            assembler.run()
            
        '''
        Binning
        '''
        pathlib.Path(args.output + '/Binning/' + mg_name).mkdir(parents=True, exist_ok=True)
        if not args.no_binning:
            print('Binning has begun')
            binner = Binner(output = args.output + '/Binning/' + mg_name,
                            contigs = args.output + '/Assembly/' + mg_name + '/contigs.fasta',
                            blast = args.output + '/Annotation/' + mg_name + '/aligned.blast',
                            uniprotinfo = args.output + '/Annotation/' + mg_name + 'uniprot.info')
            binner.run()
        
        '''
        Annotation
        '''
        pathlib.Path(args.output + '/Annotation/' + mg_name).mkdir(parents=True, exist_ok=True)
        if not args.no_annotation:
            assembled = False if args.no_assembly else True
            print('Annotation has begun')
            annotater = Annotater(out_dir = args.output,
                                 assembler = args.assembler,
                                 db = args.annotation_database,
                                 assembled = assembled,
                                 error_model = 'illumina_10',
                                 name = mg_name)
            annotater.run()
            
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
            print('Analysis has begun')
            
            path = pathlib.Path(args.output + '/Analysis/' + mt_name).mkdir(parents=True, exist_ok=True)        
            contigs = args.output + '/Assembly/' + mg_name + ('/contigs.fasta' if args.assembler == 'metaspades' else '/final.contigs.fa')
            
            mta = MetaTranscriptomicsAnalyser(out_dir = args.output,
                                              contigs = contigs,
                                              mg = mg_name,
                                              mt = mt_name,
                                              assembler = args.assembler)
            mta.readcounts_file()
            
        else:
            '''
            MetaProteomics Analysis
            '''
            spectra_folder = pairs[1]
            analyser = MetaProteomicsAnalyser(faa = 'MGMP/Annotation/' + mg_name + '/fgs.faa',
                          blast = 'MGMP/Annotation/' + mg_name + '/aligned.blast',
                          crap_folder = '/HDDStorage/jsequeira/Thesis/MetaProteomics',
                          output = 'MGMP/Analysis/' + mg_name,
                          protease = 'trypsin',
                          spectra_folder = 'Datasets/Proteic/' + spectra_folder,
                          experiment_name = args.output,
                          sample_name = mg_name,
                          replicate_number = '1')
        
            analyser.standard_run()