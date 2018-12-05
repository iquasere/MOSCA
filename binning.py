# -*- coding: utf-8 -*-
"""
MOSCA's Binning package for clustering of 
contigs into Operational Taxonomic Units

By João Sequeira

Nov 2018
"""

from sklearn.cluster import DBSCAN
from mosca_tools import MoscaTools
from diamond import DIAMOND
from progressbar import ProgressBar
import os
import pandas as pd

mtools = MoscaTools()

class Binning:
    
    def __init__ (self, **kwargs):
         self.__dict__ = kwargs
    
    '''
    *** The following functions concern binning with MyCC ***
    '''
    
    '''
    Input: 
        cluster_fasta: name of a FASTA file containing the contigs of a cluster
        blast: name of the annotation BLAST file with UniProt IDs
        uniprotinfo: name of the TSV file containing UniProt information about 
        the IDs identified
        number: a string which will identify the specific cluster. If not specified,
        it will be determined from the cluster_fasta file name
    Output:
        a pandas.DataFrame object, multi-indexed as (cluster, taxa_level), with
        columns [major_taxa, relative abundance]
    '''
    def cluster_validation(self, cluster_fasta, blast, uniprotinfo, number = None):
        if number is None:
            number = cluster_fasta.split('Cluster.')[-1].split('.fasta')[0]
        cluster_fasta = mtools.parse_fasta(cluster_fasta)
        blast = DIAMOND(out = blast).parse_result()
        uniprotinfo = pd.read_csv(uniprotinfo, sep = '\t')
        blast.index = [ide.split('|')[-1] for ide in blast.sseqid]
        blast['Node'] = ['_'.join(ide.split('_')[0:6]) for ide in blast.qseqid]
        blast = pd.merge(blast, uniprotinfo, left_index=True, right_on = ['Query'], how = 'outer')
        tax_columns = ['Taxonomic lineage (SUPERKINGDOM)','Taxonomic lineage (PHYLUM)',
                       'Taxonomic lineage (CLASS)','Taxonomic lineage (ORDER)',
                       'Taxonomic lineage (FAMILY)','Taxonomic lineage (GENUS)',
                       'Taxonomic lineage (SPECIES)']
        altered_tax_col = [str(i+1) + '.' + tax_columns[i] for i in range(len(tax_columns))]
        taxa = blast[blast['Node'].isin(cluster_fasta.keys())]
        taxa['count'] = [float(ide.split('_')[5]) for ide in taxa.qseqid]
        taxa = taxa[tax_columns + ['count']].groupby(tax_columns)['count'].sum().reset_index()
        major_taxon_abundance = dict()
        for i in range(len(tax_columns)):
            df = taxa[[tax_columns[i], 'count']]
            df = df.groupby(tax_columns[i])['count'].sum().reset_index()
            major_taxon = df[df['count'] == df['count'].max()][tax_columns[i]].tolist()            #if there is more than one taxon with maximum abundance, major_taxon will be list of them all
            major_taxon_abundance[('Cluster ' + number, altered_tax_col[i])] = [', '.join(major_taxon),             #produces a dictionary that for each contig, for each taxon level, has the value of [taxa, relative abundance of that taxa]
                                  df['count'].max() / df['count'].sum()]
        return major_taxon_abundance
    
    def parse_cluster_validation(self, file):
        return pd.read_excel(file, index_col=[0,1])
    
    def get_cluster_validation(self, validation, tax_level):
        return validation.xs(tax_level, level = 1)[1].mean()
    
    # TODO - this needs tweaking
    def bin_validation(self, cluster_validation):
        cluster_result = pd.DataFrame.from_dict(cluster_result).transpose()
        for file in files:
            cluster_result = annotater.cluster_validation(file, 
                                            'MOSCAfinal/Annotation/' + sample + '/aligned.blast',
                                            'MOSCAfinal/Annotation/' + sample + '/uniprot.info',
                                            file.split('Cluster.')[-1].split('.fasta')[0])
            cluster_result = pd.DataFrame.from_dict(cluster_result).transpose()
            sample_result = pd.concat([sample_result, cluster_result])
        sample_result.to_excel('Binning/' + sample + '/binning_validation.xlsx')
    
    
    '''
    *** The following functions concern binning with VizBin ***
    '''
    
    '''
    Input: 
        cluster_fasta: name of a FASTA file containing the contigs to bin
        output: name of folder to output files
        cutoff: minimum contig length
        kmer: k-mer length
        threads: number of threads to use
        vizbin_executable: name of vizbin executable
    Output:
        output file will be generated with coordinates of VizBin binning
        new FASTA file with contigs considered on the binning will be generated
        at the output directory, and named over_cutoff_contigs.fasta
    '''
    def run_vizbin_binning(self, contigs, output, cutoff = 1000, kmer = 5, threads = 2,
                           vizbin_executable = '~/VizBin/VizBin-dist.jar'):
        print('Performing binning with VizBin')
        bashCommand = ('java -jar ' + os.path.expanduser(vizbin_executable) + ' -i ' 
                       + contigs + ' -o ' + output + '/points.txt -c ' + str(cutoff) 
                       + ' -k ' + str(kmer) + ' -t ' + str(threads))
        mtools.run_command(bashCommand)
        print('Binning performed')
        # From down here, a new FASTA file will be generated with the contigs considered in the binning
        contigs = mtools.parse_fasta(contigs)
        over_cutoff_contigs = open(output + '/over_cutoff_contigs.fasta', 'w')
        for key, value in contigs.items():
            if len(value) >= cutoff:
                over_cutoff_contigs.write('>' + key + '\n' + value + '\n')

    '''
    Input: 
        coords: name of a coords file containing the coordinates of VizBin binning
        output: name of points file to generate
        eps: the threshold of distance under which the two points will be considered
        in the same cluster
        min_samples: the minPts heuristic, with default to 4 following the estimation
        of minPts = 2 · dim, where dim will always be 2 in this system of coordinates
    Output:
        returns a pandas.DataFrame with columns [contig, cluster]
    '''
    def cluster_coords(self, contigs, coords, eps = 0.3, min_samples = 4):
        contigs = [k for k in mtools.parse_fasta(contigs).keys()]
        contigs = pd.DataFrame(contigs)
        coords = pd.read_csv(coords, header = None)
        coords = coords.as_matrix()
        db = DBSCAN(eps, min_samples).fit(coords)
        clusters = pd.DataFrame(db.labels_.tolist())
        result = pd.concat([contigs, clusters], axis = 1)
        result.columns = ['contig','cluster']
        return result
    
    '''
    Input: 
        contigs_clusters: the pd.DataFrame object from Binning.cluster_coords
        blast: name of the annotation BLAST file with UniProt IDs
        uniprotinfo: name of the TSV file containing UniProt information about 
        the IDs identified
        by: if the averages are to be calculated by contig of by cluster (default)
        I haven't tested if this is the same thing, could be! No one really reads
        documentation nowadays... I'm watching the test of the Binning.calculate_epsilon
        function and it is fun to watch! Number of contigs going up and down...
        Well, betcha this was the most interesting part of documentation you have
        read today. Nice one, ha? ;)
    Output:
        returns a pd.DataFrame object with the most abundant taxa and its
        relative abundance for each cluster
    '''
    def estimate_mistake(self, contigs_clusters, blast, uniprotinfo, by = 'cluster'):
        blast = DIAMOND(out = blast).parse_result()
        blast['contig'] = ['_'.join(ide.split('_')[:6]) for ide in blast['qseqid']]
        blast['Query'] = [ide.split('|')[-1] for ide in blast.sseqid]
        contigs_clusters = pd.merge(contigs_clusters, blast, on='contig')
        uniprotinfo = pd.read_csv(uniprotinfo, sep = '\t')
        uniprotinfo = uniprotinfo.drop_duplicates()
        contigs_clusters = pd.merge(contigs_clusters, uniprotinfo, on='Query')
        contigs_clusters['count'] = [float(ide.split('_')[5]) for ide in contigs_clusters.qseqid]
        tax_columns = ['Taxonomic lineage (SUPERKINGDOM)','Taxonomic lineage (PHYLUM)',
                       'Taxonomic lineage (CLASS)','Taxonomic lineage (ORDER)',
                       'Taxonomic lineage (FAMILY)','Taxonomic lineage (GENUS)',
                       'Taxonomic lineage (SPECIES)']
        altered_tax_col = [str(i+1) + '.' + tax_columns[i] for i in range(len(tax_columns))]
        major_taxon_abundances = dict()
        print('Retrieving relative abundances of major taxa detected in the contigs of each cluster.')
        if by == 'cluster':
            pbar = ProgressBar()
            for cluster in pbar(set(contigs_clusters.cluster)):
                for i in range(len(tax_columns)):
                    df = contigs_clusters[contigs_clusters.cluster == cluster][[tax_columns[i], 'count']]
                    df = df.groupby(tax_columns[i])['count'].sum().reset_index()
                    major_taxon = df[df['count'] == df['count'].max()][tax_columns[i]].tolist()            #if there is more than one taxon with maximum abundance, major_taxon will be list of them all
                    major_taxon_abundances[(cluster, altered_tax_col[i])] = [', '.join(major_taxon),             #produces a dictionary that for each contig, for each taxon level, has the value of [taxa, relative abundance of that taxa]
                                          df['count'].max() / df['count'].sum()]
            # TODO - check if this is worth it, to do by contig
        elif by == 'contig':
            for cluster in set(contigs_clusters.cluster):
                print('Cluster ' + str(cluster) + '.')
                pbar = ProgressBar()
                for contig in pbar(set(contigs_clusters[contigs_clusters.cluster == cluster].contig)):
                    contigdf = contigs_clusters[contigs_clusters.contig == contig]
                    for i in range(len(tax_columns)):
                        df = contigdf[[tax_columns[i], 'count']]
                        df = df.groupby(tax_columns[i])['count'].sum().reset_index()
                        major_taxon = df[df['count'] == df['count'].max()][tax_columns[i]].tolist()            #if there is more than one taxon with maximum abundance, major_taxon will be list of them all
                        major_taxon_abundances[(cluster, altered_tax_col[i])] = ['; '.join(major_taxon),             #produces a dictionary that for each contig, for each taxon level, has the value of [taxa, relative abundance of that taxa]
                                              df['count'].max() / df['count'].sum()]
            else:
                print("'by' must be either 'cluster' (default) or 'contig'.")
            
        return pd.DataFrame(major_taxon_abundances).transpose()

    '''
    Input: 
        major_taxon_abundances: the pd.DataFrame originated from Binning.estimate_mistake
        eps: the eps value used
    Output:
        TWO objects: INT of the number of clusters obtained and a pd.DataFrame
        object with the mean, standard deviation and fragmentation for each taxa 
        level for each cluster
    '''
    def calculate_clustering_metrics(self, major_taxon_abundances, eps):
        n_clusters = round(len(major_taxon_abundances) / 7)
        major_taxa_abundance_metrics = dict()
        for tax_level in set(major_taxon_abundances.index.get_level_values(1)):
            mean_major_taxa_abundance = major_taxon_abundances.xs(tax_level, level = 1)[1].mean()
            std_major_taxa_abundance = major_taxon_abundances.xs(tax_level, level = 1)[1].std()
            fragmentation = 1 - (len(set(major_taxon_abundances.xs(tax_level, level = 1)[0])) /
                            len(major_taxon_abundances.xs(tax_level, level = 1)[0]))            # fragmentation will be 0 if the number of different taxa equals the number of contigs
            major_taxa_abundance_metrics[(eps, tax_level)] = [mean_major_taxa_abundance, 
                                                std_major_taxa_abundance, fragmentation]
        metrics_df = pd.DataFrame.from_dict(major_taxa_abundance_metrics).transpose()
        metrics_df.columns = ['Mean', 'Std dev', 'Fragmentation']
        return n_clusters, metrics_df
    
    '''
    Input: 
        contigs: name of the FASTA file with contigs generated in Binning.run_vizbin_binning
        coords: name of the TXT points file generated in Binning.run_vizbin_binning
        blast: name of the annotation BLAST file with UniProt IDs
        uniprotinfo: name of the TSV file containing UniProt information about 
        the IDs identified
        output: folder where to create result files
        start, end, step: the eps values to search for
        safety_multiplier: python's range doesn't allow floats, so this is used
        as a workaround. Only have to change it if using very low eps values!
    Output:
        This function will iterate over each eps from 'start' to 'end' 'step' by 'step'
        and for each it will bin the contigs based on the coordinates given, using the
        DBSCAN algorithm. Then it will calculate the success of the binning based on 
        the average percentage of taxa in each cluster, thus trying to minimize contamination.
        It will output an Excel file with several metrics, in several sheets:
            -the major taxa and its relative abundance for each cluster for the
            eps that was considered the best (sheet "Abundances for eps *best eps*)
            -mean, standard deviation and fragmentation for each taxa level for each 
            cluster. Fragmentation = 1 - number of different taxa / number of cluster,
            so if for each cluster there is one different taxa, fragmentation = 0.
            This value has no expression for superkingdom level, but from then, 
            it shows how much were the different taxa differentiated - if we get 
            10 clusters on a sample with 10 different phylum, it could be possible 
            to obtain fragmentation = 0.
            -number of clusters obtained and % of contigs used for each binning.
            % of contigs used = 100 * number of contigs not on cluster -1 / 
            number of contigs
        It will also output a TSV file with the contigs and corresponding clusters
    '''
    def calculate_epsilon(self, contigs, coords, blast, uniprotinfo, output, 
                          start = 0.01, end = 0.31, step = 0.01, safety_multiplier = 10000):
        best_success = 0
        [start, end, step] = map(lambda x: int(x*safety_multiplier), [start, end, step])
        all_major_taxa_abundance_metrics = pd.DataFrame()
        ns_clusters = list()
        for eps in range(start, end, step):
            print('Calculating for eps = ' + str(eps / safety_multiplier))
            contigs_clusters = self.cluster_coords(contigs, coords, eps = eps / safety_multiplier)
            clusters_used = round(sum(contigs_clusters.cluster != -1) / len(contigs_clusters) * 100, 2)
            major_taxon_abundances = self.estimate_mistake(contigs_clusters, blast, uniprotinfo)
            n_clusters, major_taxa_abundance_metrics = self.calculate_clustering_metrics(major_taxon_abundances, str(eps / safety_multiplier))
            print(str(n_clusters) + ' clusters were obtained.')
            ns_clusters.append([eps / safety_multiplier, n_clusters, clusters_used])
            all_major_taxa_abundance_metrics = pd.concat([all_major_taxa_abundance_metrics, major_taxa_abundance_metrics])
            success = float(major_taxa_abundance_metrics.xs('6.Taxonomic lineage (GENUS)', level = 1)['Mean'] *
                            clusters_used * major_taxa_abundance_metrics['Fragmentation'])
            if clusters_used > 50 and success > best_success:
                best_success = success
                best_eps = eps
                best_major_taxon_abundances = major_taxon_abundances
                best_clusters = contigs_clusters['contig','cluster']
        ns_clusters = pd.DataFrame(ns_clusters)
        ns_clusters.columns = ['Epsilon', 'Number of clusters', '% of contigs used']
        print('Best eps was ' + str(best_eps / safety_multiplier))
        best_clusters_output = output + '/best_clusters.tsv'
        best_clusters.to_csv(best_clusters_output, sep = '\t')
        print('Best clusters are outputed to ' + best_clusters_output)
        writer = pd.ExcelWriter(output + '/binning_results.xlsx', engine='xlsxwriter')
        best_major_taxon_abundances.index.names = ['Eps','Taxa level']
        best_major_taxon_abundances.columns = ['Taxa', 'Relative abundance']
        best_major_taxon_abundances.to_excel(writer, sheet_name = 'Abundances for eps ' + str(best_eps / safety_multiplier))
        all_major_taxa_abundance_metrics.index.names = ['Eps','Taxa level']
        all_major_taxa_abundance_metrics.to_excel(writer, sheet_name = 'Validation metrics')
        ns_clusters.to_excel(writer, sheet_name = 'Cluster metrics', index = False)
        writer.save()
        print('Taxonomic results of best binning are available, along with binning validation metrics, at ' + output)
        
    def run(self):
        self.run_vizbin_binning(self.contigs, self.output)
        self.calculate_epsilon(self.contigs, self.output + '/points.txt', 
                               self.blast, self.uniprotinfo, self.output)
        
if __name__ == '__main__':
    
    binner = Binning(contigs = 'Binning/all_info_contigs.fasta', 
                     output = 'Binning/VizBin/all_info',
                     blast = 'Binning/all_info_aligned.blast',
                     uniprotinfo = 'MOSCAfinal/Annotation/uniprot.info')
    
    binner.run()

    '''
    clusters = binner.cluster_coords('Binning/VizBin/' + sample + '/over_cutoff_contigs.fasta',
                          'Binning/VizBin/' + sample + '/points.txt')
    
    clusters.to_csv('Binning/VizBin/' + sample + '/clusters.tsv', sep = '\t')
    '''
