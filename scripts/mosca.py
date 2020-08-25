#!/usr/bin/env python
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
from report import Reporter
from time import gmtime, strftime
import argparse, pathlib, os, multiprocessing, pandas as pd, sys, numpy as np, glob
from tqdm import tqdm

__version__ = '1.1.1'

parser = argparse.ArgumentParser(description="Meta-Omics Software for Community Analysis",
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
parser.add_argument("-rd","--resources-directory",type=str,
                    help="Directory for storing temporary files and databases",
                    default = os.path.expanduser('~/resources_directory'))
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
parser.add_argument("-anncols", "--annotation-columns", type = str, default = '',
                    help="""List of UniProt columns to obtain information from""")
parser.add_argument("-anndbs", "--annotation-databases", type = str, default = '',
                    help="""List of databases to cross-check with UniProt information""")
parser.add_argument("-assstrat", "--assembly-strategy", type = str, default = 'all',
                    help="""'all' - all MG data is assembled in a single assembly\n
                    'unique' - for each MG sample, an assembly is performed\n
                    'samples' - requires the --mg-samples argument. Defines the 
                    communities corresponding to each MG sample, MG samples from 
                    the same community will be assembled together.""",
                    choices=['all', 'unique', 'samples'])
parser.add_argument("-samp", "--mg-samples", type = str,
                    help="""Groups of MG data to assemble together, separated by 
                    comma (,) (e.g. Sample1,Sample1,Sample2,Sample2). If not 
                    specified, data will be considered as coming from a single 
                    community, and assembled together in one assembly.""")
parser.add_argument("-mpw","--metaproteomics-workflow", default = 'compomics', type = str,
                    choices = ["compomics","maxquant"], action = 'store',
                    help=("""MOSCA possesses two workflows for metaproteomics,
                          one is an adaptation of software developed at Compomics
                          group and the other an integration of MaxQuant"""))
parser.add_argument('-v', '--version', action='version', version='MOSCA ' + __version__)

mtools = MoscaTools()
args = mtools.validate_arguments(parser)

reporter = Reporter(lists_dir = sys.path[0] + '/../Databases/reporter')
reporter.initialize_report()

reporter.write_technical_report(args.output + '/technical_report.txt')
print('Versions of the tools used is available at {}/technical_report.txt'.format(args.output))

mtools.print_arguments(args)                                                    # TODO - refine this function
monitorization_file = args.output + '/monitorization_report.txt'

experiments = [experiment.split(':') for experiment in args.files]
mtools.timed_message('MOSCA analysis has started')

print('Creating directories at ' + args.output)
directories = ([args.output + '/Preprocess/' + software for software in 
                ['FastQC', 'Trimmomatic', 'SortMeRNA']] + 
               [args.output + folder for folder in ['/Assembly','/Annotation',
                '/Metatranscriptomics' if args.type_of_data == 'metatranscriptomics' 
                else '/Metaproteomics', args.resources_directory]])

for directory in directories:
    print('Created ' + directory)
    pathlib.Path(directory).mkdir(parents=True, exist_ok=True)

mg_preprocessed = list()
mg_names = list()
mt_preprocessed = list()
if args.type_of_data == 'metatranscriptomics':
    mt2mg = dict()
else:
    mp2mg = dict()

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
            
            preprocesser.run()
            reporter.info_from_preprocessing(args.output, mg_name, mg[0])          
            mtools.remove_preprocessing_intermediates(args.output + '/Preprocess', args.output_level)
            
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
                reporter.info_from_preprocessing(args.output, mt_name, mt[0],
                                                 performed_rrna_removal = True)
                mtools.remove_preprocessing_intermediates(args.output + '/Preprocess', args.output_level)
                
                mt_preprocessed.append(mt_name)
                    
                mt2mg[mt_name] = mg_name

mtools.task_is_finished(task = 'Preprocessing', 
                        file = monitorization_file, 
                        task_output = args.output + '/Preprocess')

if args.assembly_strategy == 'all':
    sample2name = {'Sample': [mg_name for mg_name in mg_preprocessed]}
elif args.assembly_strategy == 'unique':
    samples = {'Sample' + str(i) : mg_names[i] for i in range(len(mg_preprocessed))}
elif args.assembly_strategy == 'samples':
    samples = args.mg_samples.split(',')
    sample2name = dict()
    for i in range(len(mg_names)):
        if samples[i] in sample2name.keys():
            sample2name[samples[i]].append(mg_names[i])
        else:
            sample2name[samples[i]] = [mg_names[i]]
        for k in sample2name.keys():
            sample2name[k] = set(sample2name[k])
name2sample = {vx : k for k, v in sample2name.items() for vx in v}

mg2mt = dict()
for mt_name, mg_name in mt2mg.items():
    if mg_name in mg2mt.keys():
        mg2mt[mg_name].append(mt_name)
    else:
        mg2mt[mg_name] = [mt_name]
reporter.set_samples(name2sample, mg2mt)

'''
Assembly
'''
if not args.no_assembly:
    mtools.timed_message('Assembling reads')
    
    for sample in sample2name.keys():
        
        forward_files = ['{}/Preprocess/Trimmomatic/quality_trimmed_{}_forward_paired.fq'.format(
                args.output, name) for name in sample2name[sample]]
        reverse_files = ['{}/Preprocess/Trimmomatic/quality_trimmed_{}_reverse_paired.fq'.format(
                args.output, name) for name in sample2name[sample]]
        
        mtools.run_command('cat ' + ' '.join(forward_files), file = '{}/Assembly/{}_forward.fastq'.format(
                args.output, sample))
        mtools.run_command('cat ' + ' '.join(reverse_files), file = '{}/Assembly/{}_reverse.fastq'.format(
                args.output, sample))
        
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
        
        assembler.run()
        reporter.info_from_assembly(args.output, sample)
        
    mtools.remove_assembly_intermediates(args.output, args.output_level, sample2name.keys())
    
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
        
        annotater.run()
        
        # Functional annotation with reCOGnizer
        mtools.run_command('recognizer.py -f {0}/Annotation/{1}/fgs.faa -o {0}/Annotation/{1} -t {2} --tsv -rd {3}'.format(
                args.output, sample, str(args.threads), args.resources_directory))
        
        reporter.info_from_annotation(args.output, sample)
    mtools.remove_annotation_intermediates(args.output, args.output_level, 
                                               sample2name.keys())
    
mtools.task_is_finished(task = 'Annotation',
        file = monitorization_file, 
        task_output = args.output + '/Annotation')

# Retrieval of information from UniProt IDs with UPIMAPI
mtools.timed_message('Retrieving information from UniProt IDs.')
blast_files = glob.glob(args.output + '/Annotation/*/aligned.blast')
ids = list()
for file in blast_files: ids += mtools.parse_blast(file)['sseqid'].tolist()
ids = list(set(ids)); ids.remove('*')
with open('{0}/Annotation/ids.txt'.format(args.output), 'w') as f: f.write('\n'.join(ids))
ann_cols = args.annotation_columns.split(','); ann_dbs = args.annotation_databases.split(',')
mtools.run_command('upimapi.py -i {0}/Annotation/ids.txt -o {0}/Annotation/uniprotinfo --full-id{1}{2}'.format(
        args.output, ' -anncols {}'.format(ann_cols) if ann_cols != [''] else '',     # if columns are set, they will be inputed
        ' -anndbs {}'.format(ann_dbs) if ann_dbs != [''] else ''))               # if databases are set, they will be inputed

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
        reporter.info_from_binning(args.output, sample)

mtools.task_is_finished(task = 'Binning',
        file = monitorization_file, 
        task_output = '{}/Binning/{}'.format(args.output, sample))

'''
Metatranscriptomics and metaproteomics quantification
'''
expression_analysed = list()

if len(experiment[0]) > 1:                                                     # this forces all analysis to be the same, but it is not the goal of MOSCA. Must allow for different types of analysis (with and without MT/MP, single/paired-ended, ...)
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
                          contigs = '{}/Assembly/{}/contigs.fasta'.format(args.output, name2sample[mg_name]),
                          blast = '{}/Annotation/{}/aligned.blast'.format(args.output, name2sample[mg_name]),
                          reads = mt,
                          mt = mt_name,
                          threads = args.threads)
            
            mta.readcounts_file()
            reporter.info_from_alignment(args.output, mt_name)
            
            expression_analysed.append(mt_name)
            
        readcount_files = ['{}/Metatranscriptomics/{}.readcounts'.format(args.output, mt)
                            for mt in expression_analysed]

        mtools.task_is_finished(task = 'Metatranscriptomics analysis',
              file = monitorization_file, 
              task_output = args.output + '/Metatranscriptomics')
            
    else:
        '''
        MetaProteomics Analysis
        '''
        mtools.timed_message('Analysing gene expression with metaproteomics')
        
        spectra_folder = experiment[1]
        (mg, mg_name) = mtools.process_argument_file(experiment[0], 'mg', 
                                    args.output, args.sequencing_technology) 
        mp_name = spectra_folder.split('/')[-1]
        path = pathlib.Path(args.output + '/Metaproteomics/' + mp_name).mkdir(parents=True, exist_ok=True)   
        
        analyser = MetaproteomicsAnalyser(faa = args.output + '/Annotation/' + mg_name + '/fgs.faa',
                      blast = args.output + '/Annotation/' + mg_name + '/aligned.blast',
                      crap_folder = '/HDDStorage/jsequeira/Thesis/Metaproteomics',
                      output = args.output + '/Metaproteomics/' + mp_name,
                      protease = 'trypsin',
                      spectra_folder = spectra_folder,
                      experiment_name = args.output,
                      sample_name = mg_name,
                      replicate_number = '1',
                      threads = args.threads,
                      workflow = args.metaproteomics_workflow)
    
        analyser.run()
        
        mtools.task_is_finished(task = 'Metaproteomics analysis',
                file = monitorization_file, 
                task_output = args.output + '/Metaproteomics/' + mt_name)
        
'''
Join all information in Protein and Entry Reports
'''
taxonomy_columns = ['Taxonomic lineage (' + level + ')' for level in 
        ['SUPERKINGDOM', 'PHYLUM', 'CLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES']]
functional_columns = ['COG general functional category', 'COG functional category',
       'COG protein description', 'cog']

for sample in sample2name.keys():
    mtools.timed_message('Joining data for sample: {}'.format(sample))
    
    # Join BLAST and reCOGnizer outputs
    data = pd.merge(mtools.parse_blast('{}/Annotation/{}/aligned.blast'.format(args.output, sample)),
                    pd.read_csv('{}/Annotation/{}/protein2cog.tsv'.format(args.output, sample), sep = '\t'), on = 'qseqid',
                    how = 'left')
    data['sseqid'] = [ide.split('|')[1] if ide != '*' else ide for ide in data['sseqid']]
    data.columns = [column.replace('_x', ' (DIAMOND)').replace('_y', ' (reCOGnizer)') for column in data.columns]  # after merging, BLAST columns will be repeated, and this makes explicit their origin
    data.columns = ['EC number (reCOGnizer)' if column == 'EC number' else column for column in data.columns]
    # Add UniProt information
    uniprotinfo = pd.read_csv('{}/Annotation/uniprotinfo.tsv'.format(args.output), sep = '\t')
    data = pd.merge(data, uniprotinfo, left_on = 'sseqid', right_on = 'Entry', how = 'left')
    data['Contig'] = [qseqid.split('_')[1] for qseqid in data['qseqid']]
    
    # MG quantification for each MG name of each Sample
    for mg_name in sample2name[sample]:
        
        mtools.perform_alignment('{}/Assembly/{}/contigs.fasta'.format(args.output, sample),
                ['{}/Preprocess/Trimmomatic/quality_trimmed_{}_{}_paired.fq'.format(args.output, mg_name, fr)
                for fr in ['forward', 'reverse']], '{}/Annotation/{}/{}'.format(args.output, sample, mg_name),
                threads = args.threads)
        
        mtools.normalize_mg_readcounts_by_size(
            '{}/Annotation/{}/{}.readcounts'.format(args.output, sample, mg_name), 
            '{}/Assembly/{}/contigs.fasta'.format(args.output, sample))
        data = mtools.add_abundance(data, 
                '{}/Annotation/{}/{}_normalized.readcounts'.format(args.output, sample, mg_name),
                mg_name, origin_of_data = 'metagenomics', 
                readcounts_has_tail = False)                                    # readcounts tail is removed in the normalization function
    
        for mt_name in expression_analysed:
            data = mtools.add_abundance(data, 
                '{}/Metatranscriptomics/{}.readcounts'.format(args.output, mt_name),
                mt_name, origin_of_data = 'metatranscriptomics')
    
    mtools.multi_sheet_excel('{}/MOSCA_Protein_Report.xlsx'.format(args.output),
                             data, sheet_name = sample)
    
    print('Finding consensus COG for each Entry of Sample: {}'.format(sample))
    tqdm.pandas()
    cogs_df = data.groupby('Entry')['cog'].progress_apply(
        lambda x:x.value_counts().index[0] if len(x.value_counts().index) > 0 
        else np.nan)
    cogs_df.reset_index(inplace = True)
    cogs_categories = data[functional_columns].drop_duplicates()
    
    # Aggregate information for each Entry, keep UniProt information, sum MG and MT or MP quantification
    data = data.groupby('Entry')[list(mg2mt.keys()) + list(mt2mg.keys())].sum().reset_index()
    data = pd.merge(data, uniprotinfo, on = 'Entry', how = 'left')              # couldn't groupby appropriately by uniprotinfo columns
    data = pd.merge(data, cogs_df, on = 'Entry', how = 'left')
    data = pd.merge(data, cogs_categories, on = 'cog', how = 'left')
    data = data[uniprotinfo.columns.tolist() + functional_columns + 
                list(mg2mt.keys()) + list(mt2mg.keys())]
    
    # MG normalization by sample and protein abundance
    data[mg_preprocessed].to_csv(args.output + '/mg_preprocessed_readcounts.table',
          sep = '\t', index = False)
    data = pd.concat([data, mtools.normalize_readcounts(
            args.output + '/mg_preprocessed_readcounts.table', mg_preprocessed, 
            args.output + '/mg_preprocessed_normalization_factors.txt')[[
            col + '_normalized' for col in mg_preprocessed]]], axis = 1)

    # MT normalization by sample and protein expression - normalization is repeated here because it's not stored from DESeq2 analysis
    data[expression_analysed].to_csv(args.output + '/expression_analysed_readcounts.table',
              sep = '\t', index = False)
    data = pd.concat([data, mtools.normalize_readcounts(
            args.output + '/expression_analysed_readcounts.table', expression_analysed, 
            args.output + '/expression_analysed_normalization_factors.txt')[[
            col + '_normalized' for col in expression_analysed]]], axis = 1)
    
    # For each sample, write an Entry Report
    mtools.multi_sheet_excel('{}/MOSCA_Entry_Report.xlsx'.format(args.output),
                                 data, sheet_name = sample)
    
    for mg_name in sample2name[sample]:
        # Draw the taxonomy krona plot
        data.groupby(taxonomy_columns)[mg_name].sum().reset_index()[
            [mg_name] + taxonomy_columns].to_csv('{}_{}_tax.tsv'.format(
                args.output, mg_name), sep = '\t', index = False, header = False)
        mtools.run_command('ktImportText {0}.tsv -o {0}.html'.format(
            '{}_{}_tax'.format(args.output, mg_name)))
                                                
        # Draw the functional krona plot
        data.groupby(functional_columns)[mg_name].sum().reset_index()[
            [mg_name] + functional_columns].to_csv('{}_{}_fun.tsv'.format(
                args.output, mg_name), sep = '\t', index = False, header = False)
        mtools.run_command('ktImportText {0}.tsv -o {0}.html'.format(
            '{}_{}_fun'.format(args.output, mg_name)))
    
    '''
    KEGG Pathway representations
    '''
    pathlib.Path(args.output + '/KEGG_metabolic_maps').mkdir(parents=True, exist_ok=True)
    
    mg_cols = [col + '_normalized' for col in mg2mt.keys()]
    mt_cols = [col + '_normalized' for col in mt2mg.keys()]
    data = pd.read_csv('{}/MOSCA_{}_Entry_Report.tsv'.format(args.output, sample), sep ='\t')
    
    mtools.run_command('kegg_charter.py -f {0}/MOSCA_{1}_Entry_Report.tsv -o {0}/KEGG_metabolic_maps -mgc {2} -mtc {3} -not 10 -utc -tl GENUS'.format(
        args.output, sample, ','.join(mg_cols), ','.join(mt_cols)))
    
    # TODO - solve ec numbers - UniProt vs reCOGnizer vs KEGGCharter
        
'''
Gene Expression Differential Analysis
'''
mtools.timed_message('Performing differential expression analysis.')

mta.generate_expression_matrix(readcount_files, expression_analysed, 
            args.output + '/Metatranscriptomics/all_experiments.readcounts')
mta.differential_analysis(args.output + '/Metatranscriptomics/all_experiments.readcounts', 
                args.conditions[0].split(','), args.output + '/Metatranscriptomics/')

reporter.info_from_differential_expression(args.output, 'Sample')               # Differential expression analysis will also have to be by sample

reporter.report.to_excel(args.output + '/MOSCA_General_Report.xlsx')

mtools.timed_message('MOSCA results written to {0}/MOSCA_Protein_Report.tsv and {0}/MOSCA_Entry_Report.xlsx'.format(args.output))
mtools.timed_message('Analysis with MOSCA was concluded with success!')
