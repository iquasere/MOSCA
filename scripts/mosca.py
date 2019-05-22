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
from time import gmtime, strftime

import argparse, pathlib, os, glob, multiprocessing

mosca_dir = os.path.dirname(os.path.realpath(__file__))

mtools = MoscaTools()

parser = argparse.ArgumentParser(description="Multi Omics Software for Community Analysis",
                                 epilog="""A tool for performing metagenomics, metatranscriptomics 
                                 and metaproteomics analysis.""")
parser.add_argument("-f","--files", type=str, nargs = '*', metavar = 'Input files',
                    help="Input files for analysis (mg1R1,mg1R2:mt1R1,mt1R2 mg2R1,...)")
parser.add_argument("-data","--type-of-data",default='paired',choices=["paired","single"],
                    action='store',type=str,help='Type of data (paired/single)-ended')
parser.add_argument("-a","--assembler",type=str,choices=["metaspades","megahit"],
                    help="Tool for assembling the reads", metavar = "Assembler", 
                    nargs = 1, default = "metaspades")
parser.add_argument("-db","--annotation-database",type=str,nargs = 1,metavar = "Database",
                    help="Database for annotation (.fasta or .dmnd)", 
                     default = "Databases/annotation_databases/uniprot.dmnd")
parser.add_argument("-o","--output",type=str,help="Directory for storing the results",
                    metavar = "Directory")
parser.add_argument("-nopp","--no-preprocessing",action = "store_true",
                    help="Don't perform preprocessing", default = False)
parser.add_argument("-noas","--no-assembly",action = "store_true",
                    help="Don't perform assembly", default = False)
parser.add_argument("-noan","--no-annotation",action = "store_true",default = False,
                    help="Don't perform annotation")
parser.add_argument("-nobin","--no-binning",action = "store_true",default = False,
                    help="Don't perform binning")
parser.add_argument("-ol","--output-level",default='min',choices=["min","med","max"],
                    action='store',type=str,help=("""Level of file output from MOSCA, 
                    min outputs only the analysis results, med removes intermediate files, 
                    max outputs all intermediate and final data"""))
parser.add_argument("-mp", "--metaproteomic", action = "store_true",
                    help = ("""If data is metagenomic and metaproteomic, if not specified
                            will be assumed to be metagenomic and metatranscriptomic"""),
                    default = False)
parser.add_argument("-c","--conditions", type=str, nargs = '*',
                    help="Different conditions for metatranscriptomic analysis, separated by comma (,)")
parser.add_argument("-t","--threads",type=str,metavar = "Threads", default = str(multiprocessing.cpu_count() - 2),
                    help="Number of threads available for MOSCA")
parser.add_argument("-m","--memory",type=str,metavar = "Memory",
                    help="Maximum memory (byte) available for MOSCA. Applied only in the assembly")
             
args = mtools.validate_arguments(parser)

mtools.print_arguments(args)            # TODO - refine this function

experiments = [experiment for experiment in args.files]

print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Analysis has started')

print('Creating directories at ' + args.output)
directories = ([args.output + '/Preprocess/' + software for software in 
                ['FastQC', 'Trimmomatic', 'SortMeRNA']] + 
               [args.output + folder for folder in ['/Assembly','/Annotation',
                '/Metatranscriptomics' if not args.metaproteomic else
                '/Metaproteomics']])

for directory in directories:
    print('Created ' + directory)
    path = pathlib.Path(directory)
    path.mkdir(parents=True, exist_ok=True)

mg_processed = list()
expression_analysed = list()

for experiment in experiments:
    pairs = experiment.split(':')
    mg = pairs[0].split(',')
    
    if len(mg) == 1 and args.type_of_data == 'paired':                           # if data is interleaved paired end, it will be split up to forward and reverse files
        (forward, reverse) = (args.output + '/Preprocess/' + mg[0].split('/')[-1].split('.fastq')[0]
                                + fr for fr in ['_R1.fastq','_R2.fastq'])
        print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Splitting reads at ' 
              + mg[0] + ' to ' + forward + ' and ' + reverse)
        mtools.divide_fq(mg[0], forward, reverse)
        mg = [forward, reverse]
            
    mg_name = mg[0].split('/')[-1].split('_R')[0]
    
    if mg_name not in mg_processed:                                             #several MT samples might correspond to the same MG sample

        '''    
        Metagenomics Preprocess        
        '''
        if not args.no_preprocessing:
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Preprocessing ' + 
                  'metagenomic reads')
            preprocesser = Preprocesser(files = mg,
                                        paired = 'PE' if args.type_of_data == 'paired' else 'SE',
                                        working_dir = args.output,
                                        data = 'dna',
                                        name = mg_name,
                                        threads = args.threads)
            if hasattr(args, 'quality_score'):
                setattr(preprocesser, 'quality_score', args.quality_score)
                
            preprocesser.run()
            
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Preprocessing ' + 
          'is finished, and resulting reads are available at ' + args.output + 
          '/Preprocess/' + mg_name)
            
            mg = [args.output + '/Preprocess/Trimmomatic/quality_trimmed_' + mg_name + 
                  '_' + fr + '_paired.fq' for fr in ['forward', 'reverse']]
        
        '''
        Assembly
        '''
        
        if not args.no_assembly:
            pathlib.Path(args.output + '/Assembly/' + mg_name).mkdir(parents=True, exist_ok=True)
            
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Assembling reads')
            
            mid_name = args.output + '/Preprocess/Trimmomatic/quality_trimmed_'
            
            assembler = Assembler(out_dir = args.output,
                                 assembler = args.assembler,
                                 name = mg_name,
                                 forward = mg[0],
                                 reverse = mg[1],
                                 threads = args.threads)
            
            if args.assembler == 'metaspades' and hasattr(args, 'quality_score'):   # Megahit doesn't accept quality score input
                setattr(assembler, 'phred_offset', args.quality_score)              # --phred-offset is the name of the parameter in MetaSPAdes
            if args.memory is not None:
                setattr(assembler, 'memory', args.memory)
            
            assembler.run()
            
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Assembly is ' + 
                  'finished. Contigs are available at ' + args.output + '/Assembly/' + 
                  mg_name + '/contigs.fasta')
        
        '''
        Annotation
        '''
        if not args.no_annotation:
            pathlib.Path(args.output + '/Annotation/' + mg_name).mkdir(parents=True, exist_ok=True)
            assembled = False if args.no_assembly else True
            
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Annotating sequences')
            annotater = Annotater(file = args.output + '/Assembly/' + mg_name + '/contigs.fasta',
                                 out_dir = args.output,
                                 assembler = args.assembler,
                                 db = args.annotation_database[0],
                                 assembled = assembled,
                                 error_model = 'illumina_10',
                                 name = mg_name,
                                 threads = args.threads)
            annotater.run()
            
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Annotation is ' + 
                  'finished and results are available at ' + args.output +
                  '/Annotation/' + mg_name)
            
        '''
        Binning
        '''
        if not args.no_binning:
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Binning has begun')
            pathlib.Path(args.output + '/Binning/' + mg_name).mkdir(parents=True, exist_ok=True)
            
            binner = Binner(output = args.output + '/Binning/' + mg_name,
                            contigs = args.output + '/Assembly/' + mg_name + '/contigs.fasta',
                            blast = args.output + '/Annotation/' + mg_name + '/aligned.blast',
                            uniprotinfo = args.output + '/Annotation/' + mg_name + '/uniprot.info',
                            threads = args.threads)
            binner.run()
            
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Binning is ' + 
                  'finished and results are available at ' + args.output + 
                  '/Binning/' + mg_name)
            
        mg_processed.append(mg_name)
        
    if len(pairs) > 1:
        
        if not args.metaproteomic:
            
            ''''
            Metatranscriptomics Preprocessing
            '''
            
            mt = pairs[1].split(',')
            
            if len(mt) == 1 and args.type_of_data == 'paired':                           # if data is interleaved paired end, it will be split up
            
                mt_name = mt[0].split('/')[-1].split('.fastq')[0]
                
                (forward, reverse) = (args.output + '/Preprocess/' + mt[0].split('/')[-1].split('.fastq')[0]
                                        + fr for fr in ['_R1.fastq','_R2.fastq'])
                print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Splitting ' + 
                      'reads at ' + mt[0] + ' to ' + forward + ' and ' + reverse + '.')
                mtools.divide_fq(mt[0], forward, reverse)
                mt = [forward, reverse]
                
            else:
                mt_name = mt[0].split('/')[-1].split('_R')[0]
                
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Performing ' + mt_name + 
                  ' analysis')
                
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Preprocessing ' + 
                  'metatranscriptomics reads')
            
            preprocesser = Preprocesser(files = mt,
                                        paired = 'PE',
                                        working_dir = args.output,
                                        data = 'mrna',
                                        name = mt_name,
                                        threads = args.threads)
            if hasattr(args, 'quality_score'):
                setattr(preprocesser, 'quality_score', args.quality_score)
                
            preprocesser.run()
            
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Preprocessing ' + 
                  'of metatranscriptomics reads is finished and results are ' + 
                  'available at ' + args.output + '/Preprocess/' + mt_name)
            
            '''
            Metatranscriptomics Quantification
            '''
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Quantifying gene expression')
            
            path = pathlib.Path(args.output + '/Metatranscriptomics/' + 
                                mg_name).mkdir(parents=True, exist_ok=True)
            
            mt = [args.output + '/Preprocess/Trimmomatic/quality_trimmed_' + mt_name + 
                  '_' + fr + '_paired.fq' for fr in ['forward', 'reverse']]
            
            mta = MetaTranscriptomicsAnalyser(out_dir = args.output + '/Metatranscriptomics',
                          contigs = args.output + '/Assembly/' + mg_name + '/contigs.fasta',
                          blast = args.output + '/Annotation/' + mg_name + '/aligned.blast',
                          reads = mt,
                          mg = mg_name,
                          mt = mt_name,
                          assembler = args.assembler,
                          threads = args.threads)
            mta.readcounts_file()
            
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Gene expression ' + 
                  'quantification is finished and results are available at ' + 
                  args.output + '/Metatranscriptomics/' + mt_name)
            
            expression_analysed.append(mt_name)
            
        else:
            
            '''
            MetaProteomics Analysis
            '''
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Analysing ' + 
                  'protein expression')
            
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
                          replicate_number = '1',
                          threads = args.threads)
        
            analyser.standard_run()
            
            print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Gene expression ' + 
                  'quantification is finished and results are available at ' + 
                  args.output + '/MetaProteomics/' + mt_name)
            
'''
Join all information on one report
'''
print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Integrating all information.')

annotater = Annotater(out_dir = args.output,
                      fun = mosca_dir + '/Databases/COG/fun.txt',
                      cog = mosca_dir + '/Databases/COG/Cog',
                      cddid = mosca_dir + '/Databases/COG/cddid.tbl',
                      whog = mosca_dir + '/Databases/COG/whog',
                      cdd2cog_executable = mosca_dir + '/cdd2cog.pl',
                      threads = args.threads)

joined = annotater.global_information()

import pandas as pd
joined = pd.read_csv('debugMOSCA/joined_information.tsv',sep='\t')
print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Integration is available at '
      + args.output)

if len(experiment.split(':')) > 1:      
    if not args.metaproteomic:
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
print(samples)
joined[samples].to_csv(args.output + '/readcounts.table',
      sep = '\t', index = False)

joined = mtools.normalize_readcounts(args.output + '/readcounts.table', samples,
                            args.output + '/normalization_factors.txt')

joined.to_csv(args.output + '/mosca_resulsts.tsv', sep = '\t', index = False)
joined.to_excel(args.output + '/mosca_resulsts.xlsx', index = False)

print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Analysis with MOSCA was concluded with success!')