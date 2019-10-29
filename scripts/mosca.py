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
from metaproteomics_analyser import MetaproteomicsAnalyser
from kegg_pathway import KEGGPathway
from report import Reporter
from time import gmtime, strftime

import argparse, pathlib, os, glob, multiprocessing

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
                     default = "MOSCA/Databases/annotation_databases/uniprot.fasta")
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
parser.add_argument("-ol","--output-level", default = 'medium', type = str,
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
                    analysis, separated by comma (,) (e.g. c1,c1,c2,c2)""")
parser.add_argument("-t","--threads",type=str, metavar = "Threads", 
                    default = str(multiprocessing.cpu_count() - 2),
                    help="Number of threads available for MOSCA. Default is number of cores available.")
parser.add_argument("-m","--memory",type=str, default = 'None',                 # None is given here as default because for now, web interface needs it # TODO - drop dis 'None' necessity
                    help="Maximum memory (byte) available for MOSCA. Applied only in the assembly")
parser.add_argument("-mark","--marker-gene-set",type=str,
                    help="""Marker gene set to use for binning with MaxBin2. 
                    107 if archaea are not to be considered, 40 if data is diverse.""",
                    default = '40')
parser.add_argument("-tr","--fgs-train-file",type=str,
                    help="""File name that contains model parameters of sequencing,
                    related to sequencing technology and error rate""", default = 'illumina_10', 
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
parser.add_argument("-assstrat", "--assembly-strategy", type = str, default = 'all',
                    help="""'all' - all MG data is assembled in a single assembly\n
                    'unique' - for each MG sample, an assembly is performed\n
                    'samples' - requires the --mg-samples argument. Defines the 
                    communities corresponding to each MG sample, MG samples from 
                    the same community will be assembled together.""")
parser.add_argument("-samp", "--mg-samples", type = str,
                    help="""Groups of MG data to assemble together, separated by 
                    comma (,) (e.g. Sample1,Sample1,Sample2,Sample2). If not 
                    specified, data will be considered as coming from a single 
                    community, and assembled together in one assembly.""")

mtools = MoscaTools()
args = mtools.validate_arguments(parser)

mosca_dir = os.path.dirname(os.path.realpath(__file__))
reporter = Reporter(metatranscriptomics = True if args.type_of_data == 'metatranscriptomics' else False)

mtools.print_arguments(args)            # TODO - refine this function
monitorization_file = args.output + '/monitorization_report.txt'

experiments = [experiment.split(':') for experiment in args.files]

mtools.timed_message('MOSCA analysis has started')

print('Creating directories at ' + args.output)
directories = ([args.output + '/Preprocess/' + software for software in 
                ['FastQC', 'Trimmomatic', 'SortMeRNA']] + 
               [args.output + folder for folder in ['/Assembly','/Annotation',
                '/Metatranscriptomics' if args.type_of_data == 'metatranscriptomics' 
                else '/Metaproteomics']])

for directory in directories:
    print('Created ' + directory)
    pathlib.Path(directory).mkdir(parents=True, exist_ok=True)

mg_preprocessed = ['4478-DNA-S1613-MiSeqKapa','4478-DNA-S1616-MiSeqKapa']
mg_names = ['4478-DNA-S1613-MiSeqKapa', '4478-DNA-S1613-MiSeqKapa', '4478-DNA-S1616-MiSeqKapa']
mt_preprocessed = ['4478-R1-1-MiSeqKapa']
mt2mg = dict()

'''    
Preprocess
'''
if not args.no_preprocessing:
    mtools.timed_message('Preprocessing reads')

    for experiment in experiments:
        
        ''''
        Metagenomics preprocessing
        '''
        
        (mg, mg_name) = mtools.process_argument_file(experiment[0], 'mg', 
                                        args.output, args.sequencing_technology)
        
        if mg_name not in mg_preprocessed:                                         #several MT samples might correspond to the same MG sample
            preprocesser = Preprocesser(files = mg,
                                        paired = 'PE' if args.sequencing_technology == 'paired' else 'SE',
                                        working_dir = args.output,
                                        data = 'dna',
                                        name = mg_name,
                                        threads = args.threads)
            if hasattr(args, 'quality_score'):
                setattr(preprocesser, 'quality_score', args.quality_score)

            #preprocesser.run()
            '''
            reporter.info_from_preprocessing(args.output, mg_name)          
            mtools.remove_preprocessing_intermediates(args.output + '/Preprocess', 
                                                      args.output_level)
            '''
            mg = [args.output + '/Preprocess/Trimmomatic/quality_trimmed_' + mg_name + 
                  '_' + fr + '_paired.fq' for fr in ['forward', 'reverse']]
            
            mg_preprocessed.append(mg_name)
            
        mg_names.append(mg_name)                                                # some MT/MP might have the same MG
            
        if len(experiment) > 1:
            
            if args.type_of_data == 'metatranscriptomics':
                
                ''''
                Metatranscriptomics preprocessing
                '''
                (mt, mt_name) = mtools.process_argument_file(experiment[1], 'mt', 
                                        args.output, args.sequencing_technology)
                
                if mt_name not in mt_preprocessed:
                    if len(mt) == 1 and args.sequencing_technology == 'paired':                           # if data is interleaved paired end, it will be split up
                    
                        (forward, reverse) = (args.output + '/Preprocess/' + mt[0].split('/')[-1].split('.fastq')[0]
                                                + fr for fr in ['_R1.fastq','_R2.fastq'])
                        print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': Splitting ' + 
                              'reads at ' + mt[0] + ' to ' + forward + ' and ' + reverse + '.')
                        mtools.divide_fq(mt[0], forward, reverse)
                        mt = [forward, reverse]
                    
                    preprocesser = Preprocesser(files = mt,
                                                paired = 'PE',
                                                working_dir = args.output,
                                                data = 'mrna',
                                                name = mt_name,
                                                threads = args.threads)
                    if hasattr(args, 'quality_score'):
                        setattr(preprocesser, 'quality_score', args.quality_score)
                    
                    preprocesser.run()
                    
                    reporter.info_from_preprocessing(args.output, mt_name, 
                                                     performed_rrna_removal = True)
                    mtools.remove_preprocessing_intermediates(args.output + '/Preprocess', args.output_level)

                    mt_preprocessed.append(mt_name)
                    
                mt2mg[mt_name] = mg_name
                    
reporter.report.to_csv(args.output + '/report.tsv', sep = '\t')

mtools.task_is_finished(task = 'Preprocessing', 
                        file = monitorization_file, 
                        task_output = args.output + '/Preprocess')

if args.assembly_strategy == 'all':
    sample2name = {'Sample': [mg_name for mg_name in mg_preprocessed]}
elif args.assembly_strategy == 'unique':
    samples = {'Sample' + str(i) : mg_names[i] for i in range(len(mg_preprocessed))}
else:
    samples = args.mg_samples.split(',')
    sample2name = dict()
    for i in range(len(mg_names)):
        if samples[i] in sample2name.keys():
            sample2name[samples[i]].append(mg_names[i])
        else:
            sample2name[samples[i]] = [mg_names[i]]
        for k in sample2name.keys():
            sample2name[k] = set(sample2name[k])
            
reporter.set_samples(sample2name)

'''
Assembly
'''

if not args.no_assembly:
    mtools.timed_message('Assembling reads')
    
    for sample in sample2name.keys():
        '''
        forward_files = ['{}/Preprocess/Trimmomatic/quality_trimmed_{}_forward_paired.fq'.format(
                args.output, name) for name in sample2name[sample]]
        reverse_files = ['{}/Preprocess/Trimmomatic/quality_trimmed_{}_reverse_paired.fq'.format(
                args.output, name) for name in sample2name[sample]]
        mtools.run_command('cat ' + ' '.join(forward_files), file = '{}/Assembly/{}_forward.fastq'.format(
                args.output, sample))
        mtools.run_command('cat ' + ' '.join(reverse_files), file = '{}/Assembly/{}_reverse.fastq'.format(
                args.output, sample))
        '''
        pathlib.Path('{}/Assembly/{}'.format(args.output, sample)).mkdir(parents=True, exist_ok=True)
        assembler = Assembler(out_dir = '{}/Assembly/{}'.format(args.output, sample),
                             assembler = args.assembler,
                             forward = '{}/Assembly/{}_forward.fastq'.format(
                                     args.output, sample),
                             reverse = '{}/Assembly/{}_reverse.fastq'.format(
                                     args.output, sample),
                             threads = args.threads)
        
        if (args.assembler == 'metaspades' and hasattr(args, 'quality_score')   # Megahit doesn't accept quality score input
        and getattr(args, 'quality_score') not in ['None', None]):
            setattr(assembler, 'phred_offset', args.quality_score)              # --phred-offset is the name of the parameter in MetaSPAdes
        if args.memory != 'None':
            setattr(assembler, 'memory', args.memory)
        '''
        assembler.run()
        reporter.info_from_assembly(args.output, sample)
        '''
    #mtools.remove_assembly_intermediates(args.output, args.output_level, sample2name.keys())
    
mtools.task_is_finished(task = 'Assembly',
                        file = monitorization_file, 
                        task_output = args.output + '/Assembly')

'''
Annotation
'''
if not args.no_annotation:
    mtools.timed_message('Annotating sequences')
    
    for sample in sample2name.keys():
        pathlib.Path('{}/Annotation/{}'.format(args.output, sample)).mkdir(parents=True, exist_ok=True)
        annotater = Annotater(file = '{}/Assembly/{}/contigs.fasta'.format(
                args.output, sample),
                      out_dir = '{}/Annotation/{}'.format(args.output, sample),
                      db = args.annotation_database,
                      assembled = False if args.no_assembly else True,
                      error_model = args.fgs_train_file,
                      threads = args.threads)
    
        #annotater.run()
        annotater.cog_annotation('{}/Annotation/{}/fgs.faa'.format(args.output, sample),
                                 '{}/Annotation/{}'.format(args.output, sample),
                                 threads = args.threads)        
    
    mtools.remove_annotation_intermediates(args.output, args.output_level, 
                                           sample2name.keys())
    
mtools.task_is_finished(task = 'Annotation',
        file = monitorization_file, 
        task_output = args.output + '/Annotation')
    
'''
Binning
'''
if not args.no_binning:
    mtools.timed_message('Binning contigs')
    
    for sample in sample2name.keys():
        pathlib.Path('{}/Binning/{}'.format(args.output, sample)).mkdir(parents=True, exist_ok=True)
    
        binner = Binner(output = '{}/Binning/{}'.format(args.output, sample),
                    contigs = '{}/Assembly/{}/contigs.fasta'.format(
                            args.output, sample),
                    blast = '{}/Annotation/{}/aligned.blast'.format(args.output, sample),
                    uniprotinfo = '{}/Annotation/{}/uniprot_info.tsv'.format(args.output, sample),
                    threads = args.threads,
                    forward = '{}/Assembly/{}_forward.fastq'.format(
                            args.output, sample),
                    reverse = '{}/Assembly/{}_reverse.fastq'.format(
                            args.output, sample),
                    markerset = args.marker_gene_set)
    binner.maxbin_workflow()


    mtools.task_is_finished(task = 'Binning',
            file = monitorization_file, 
            task_output = '{}/Binning/{}'.format(args.output, sample))
    
'''
Join all information on one report
'''
mtools.timed_message('Integrating all information.')

annotater = Annotater(out_dir = args.output,
                      threads = args.threads,
                      samples = samples,
                      columns = (args.annotation_columns.split(',') if 
                      hasattr(args, 'annotation_columns') else None),
                      databases = (args.annotation_databases.split(',') if 
                      hasattr(args, 'annotation_databases') else None))
joined = annotater.global_information()

mtools.timed_message('Integration is available at ' + args.output)

expression_analysed = list()

if len(experiment[0].split(':')) > 1:                                           # this forces all analysis to be the same, but it is not the goal of MOSCA. Must allow for different types of analysis (with and without MT/MP, single/paired-ended, ...)
    if args.type_of_data == 'metatranscriptomics':
        '''
        Metatranscriptomics Quantification
        '''
        mtools.timed_message('Analysing gene expression with metatranscriptomics')
        
        for mt_name in mt_preprocessed:
            path = pathlib.Path('{}/Metatranscriptomics/{}'.format(args.output, 
                                mt_name)).mkdir(parents=True, exist_ok=True)
            
            mt = ['{}/Preprocess/Trimmomatic/quality_trimmed_{}_{}_paired.fq'.format(
                    args.output, mt_name, fr) for fr in ['forward', 'reverse']]
            mg_name = mt2mg[mt_name]
            mta = MetaTranscriptomicsAnalyser(out_dir = args.output + '/Metatranscriptomics',
                          contigs = '{}/Assembly/{}/contigs.fasta'.format(args.output, mg_name),
                          blast = '{}/Annotation/{}/aligned.blast'.format(args.output, mg_name),
                          reads = mt,
                          mt = mt_name,
                          threads = args.threads)
            mta.readcounts_file()
            
            joined = mtools.define_abundance(joined, origin_of_data = 'metatranscriptomics',
                                             name = mt_name, readcounts = '{}/Metatranscriptomics/{}.readcounts'.format(
                                                     args.output, mt_name))
            
        readcount_files = glob.glob(args.output + '/Metatranscriptomics/*.readcounts')
        mta.generate_expression_matrix(readcount_files, expression_analysed, 
            args.output + '/Metatranscriptomics/all_experiments.readcounts')

        mta.differential_analysis(args.output + '/Metatranscriptomics/all_experiments.readcounts', 
                        args.conditions[0].split(','), args.output + '/Metatranscriptomics/')
            
        mtools.task_is_finished(task = 'Metatranscriptomics analysis',
              file = monitorization_file, 
              task_output = args.output + '/Metatranscriptomics')
            
        expression_analysed.append(mt_name)
    
    else:
        
        '''
        MetaProteomics Analysis
        '''
        mtools.timed_message('Analysing gene expression with metaproteomics')
        
        path = pathlib.Path(args.output + '/Metaproteomics/' + mt_name).mkdir(parents=True, exist_ok=True)   
        
        spectra_folder = experiment.split(':')[1]
        analyser = MetaproteomicsAnalyser(faa = args.output + '/Annotation/' + mg_name + '/fgs.faa',
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

joined[mg_preprocessed].to_csv(args.output + '/mg_preprocessed_readcounts.table',
      sep = '\t', index = False)
joined = mtools.normalize_readcounts(args.output + '/mg_preprocessed_readcounts.table', 
            mg_preprocessed, args.output + '/mg_preprocessed_normalization_factors.txt')
joined[expression_analysed].to_csv(args.output + '/expression_analysed_readcounts.table',
      sep = '\t', index = False)
joined = mtools.normalize_readcounts(args.output + '/expression_analysed_readcounts.table', 
            expression_analysed, args.output + '/expression_analysed_normalization_factors.txt')

joined.to_csv(args.output + '/mosca_results.tsv', sep = '\t', index = False)
joined.to_excel(args.output + '/mosca_results.xlsx', index = False)

# KEGG Pathway representations
kp = KEGGPathway(input_file = args.output + '/mosca_results.tsv',
                 output_directory = args.output + '/Metatranscriptomics/KEGG Pathway',
                 mg_samples = mg_preprocessed,
                 mt_samples = expression_analysed)
kp.run()

mtools.timed_message('Analysis with MOSCA was concluded with success!')

mtools.write_technical_report(args.output + 'technical_report.txt')             # TODO - must be improved
print('Versions of the tools used is available at {}/technical_report.txt'.format(args.output))