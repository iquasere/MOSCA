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
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm


mtools = MoscaTools()

class Binner:
    
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
    def run_vizbin_binning(self, contigs, output, cutoff = 1000, kmer = 5, threads = 2):
        print('Performing binning with VizBin')
        bashCommand = ('java -jar VizBin/VizBin-dist.jar -i ' 
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
        contigs: output of mtools.parse_fasta on the contigs file
        coords: name of a coords file containing the coordinates of VizBin binning
        eps: the threshold of distance under which the two points will be considered
        in the same cluster
        min_samples: the minPts heuristic, with default to 4 following the estimation
        of minPts = 2 · dim, where dim will always be 2 in this system of coordinates
    Output:
        returns a pandas.DataFrame with columns [contig, cluster]
    '''
    def cluster_coords(self, contigs, coords, eps = 0.3, min_samples = 4):
        contigs = pd.DataFrame(contigs)
        coords = pd.read_csv(coords, header = None)
        coords = coords.values
        db = DBSCAN(eps, min_samples).fit(coords)
        clusters = pd.DataFrame(db.labels_.tolist())
        result = pd.concat([contigs, clusters], axis = 1)
        result.columns = ['contig','cluster']
        result = result[result.cluster.notnull()]
        result.cluster = result.cluster.astype(int)
        print('Clustering finished successfully!')
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
        blast['Entry'] = [ide.split('|')[1] if ide != '*' else ide for ide in blast.sseqid]
        contigs_clusters = pd.merge(contigs_clusters, blast, on='contig')
        uniprotinfo = pd.read_csv(uniprotinfo, sep = '\t')
        uniprotinfo = uniprotinfo.drop_duplicates()
        contigs_clusters = pd.merge(contigs_clusters, uniprotinfo, on='Entry')
        contigs_clusters['count'] = [float(ide.split('_')[5]) for ide in contigs_clusters.qseqid]
        tax_columns = ['Taxonomic lineage (SUPERKINGDOM)','Taxonomic lineage (PHYLUM)',
                       'Taxonomic lineage (CLASS)','Taxonomic lineage (ORDER)',
                       'Taxonomic lineage (FAMILY)','Taxonomic lineage (GENUS)',
                       'Taxonomic lineage (SPECIES)']
        altered_tax_col = [str(i+1) + '.' + tax_columns[i] for i in range(len(tax_columns))]                #tax columns are altered for sorting in final output
        major_taxon_abundances = dict()
        print('Retrieving relative abundances of major taxa detected in the contigs of each cluster.')
        if by == 'cluster':
            pbar = ProgressBar()
            for cluster in pbar(set(contigs_clusters.cluster)):
                for i in range(len(tax_columns)):
                    df = contigs_clusters[contigs_clusters.cluster == cluster][[tax_columns[i], 'count']]
                    df = df.groupby(tax_columns[i])['count'].sum().reset_index()
                    major_taxon = df[df['count'] == df['count'].max()][tax_columns[i]].tolist()             #if there is more than one taxon with maximum abundance, major_taxon will be list of them all
                    major_taxon_abundances[(cluster, altered_tax_col[i])] = [', '.join(major_taxon),        #produces a dictionary that for each contig, for each taxon level, has the value of [taxa, relative abundance of that taxa]
                                          df['count'].max() / df['count'].sum() / len(major_taxon)]         #in the case of several most abundant taxon, it will return the relative abundance of each taxon
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
            mta = major_taxon_abundances                                            # TODO - try to restrict to where the taxa is not null [major_taxon_abundances[tax_level].notnull()]
            mean_major_taxa_abundance = mta.xs(tax_level, level = 1)[1].mean()
            std_major_taxa_abundance = mta.xs(tax_level, level = 1)[1].std()
            fragmentation = 1 - (len(set(mta.xs(tax_level, level = 1)[0])) /
                            len(mta.xs(tax_level, level = 1)[0]))            # fragmentation will be 0 if the number of different taxa equals the number of contigs
            score = mean_major_taxa_abundance / fragmentation
            major_taxa_abundance_metrics[(eps, tax_level)] = [mean_major_taxa_abundance, 
                                         std_major_taxa_abundance, fragmentation, score]
        metrics_df = pd.DataFrame.from_dict(major_taxa_abundance_metrics).transpose()
        metrics_df.columns = ['Mean', 'Std dev', 'Fragmentation', 'Score']
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
        This function will iterate over each eps from 'start' to 'end', 'step' by 'step'
        and for each it will bin the contigs based on the coordinates given, using the
        DBSCAN algorithm. Then it will calculate the success of the binning based on 
        the average percentage of taxa in each cluster, the fragmentation of the clusters
        and the number of contigs clustered, for a minimum of 50% contigs clustered.
        It will output an Excel file with several metrics, in several sheets:
            -the major taxa and its relative abundance for each cluster for the
            eps that was considered the best (sheet "Abundances for eps *best eps*)
            -mean, standard deviation and fragmentation for each taxa level for each 
            cluster. Fragmentation = 1 - number of different taxa / number of cluster,
            so if for each different taxa there is one cluster, fragmentation = 0.
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
                          start = 0.01, end = 0.41, step = 0.01, safety_multiplier = 10000):
        [start, end, step] = map(lambda x: int(x*safety_multiplier), [start, end, step])
        all_major_taxa_abundance_metrics = pd.DataFrame()
        ns_clusters = list()
        contigs = [k for k in mtools.parse_fasta(contigs).keys()]
        best_eps = None
        best_score = -float('inf')
        for eps in range(start, end, step):
            print('Calculating for eps = ' + str(eps / safety_multiplier))
            contigs_clusters = self.cluster_coords(contigs, coords, eps = eps / safety_multiplier)
            contigs_used = round(sum(contigs_clusters.cluster != -1) / len(contigs_clusters) * 100, 2)
            major_taxon_abundances = self.estimate_mistake(contigs_clusters, blast, uniprotinfo)
            n_clusters, major_taxa_abundance_metrics = self.calculate_clustering_metrics(
                    major_taxon_abundances, str(eps / safety_multiplier))
            ns_clusters.append([eps / safety_multiplier, n_clusters, contigs_used])
            all_major_taxa_abundance_metrics = pd.concat([all_major_taxa_abundance_metrics, 
                                                          major_taxa_abundance_metrics])
            if contigs_used > 50 and float(major_taxa_abundance_metrics.xs(
                    '6.Taxonomic lineage (GENUS)', level = 1)['Score']) > best_score:
                best_score = float(major_taxa_abundance_metrics.xs(
                        '6.Taxonomic lineage (GENUS)', level = 1)['Score'])
                best_eps = eps
                best_major_taxon_abundances = major_taxon_abundances
                best_clusters = contigs_clusters
        if best_eps is not None:
            ns_clusters = pd.DataFrame(ns_clusters)
            ns_clusters.columns = ['Epsilon', 'Number of clusters', '% of contigs used']
            print('Best eps was ' + str(best_eps / safety_multiplier))
            best_clusters_output = output + '/best_clusters.tsv'
            best_clusters.to_csv(best_clusters_output, sep = '\t', index = False)
            print('Best clusters are outputed to ' + best_clusters_output)
            writer = pd.ExcelWriter(output + '/binning_results.xlsx', engine='xlsxwriter')
            best_major_taxon_abundances.index.names = ['Eps','Taxa level']
            best_major_taxon_abundances = best_major_taxon_abundances.sort_index(
                    level = 0, sort_remaining = True)
            best_major_taxon_abundances.columns = ['Taxa', 'Relative abundance']
            best_major_taxon_abundances.to_excel(writer, sheet_name = 'Abundances for eps ' 
                                                 + str(best_eps / safety_multiplier))
            all_major_taxa_abundance_metrics.index.names = ['Eps','Taxa level']
            all_major_taxa_abundance_metrics = all_major_taxa_abundance_metrics.sort_index(
                    level = 0, sort_remaining = True)
            all_major_taxa_abundance_metrics.to_excel(writer, sheet_name = 'Validation metrics')
            ns_clusters.to_excel(writer, sheet_name = 'Cluster metrics', index = False)
            writer.save()
            print('Taxonomic results of best binning are available, along with binning validation metrics, at ' + output)
        else:
            print('No clustering used at least 50% of contigs. Binning could not be performed.')
    
    '''
    Input: 
        clusters: name of the TSV file outputed by Binner.calculate_epsilon,
        associating a cluster to each contig
        blast: name of BLAST annotation file
        uniprotinfo: name of TSV file with UniProt information
        output: 
        taxa_level: the taxa_level to be described
        deepness: how many different taxa to report for each cluster
    Output:
        A TSV file named output will describe the main taxon (number of taxon
        is deepness value) and corresponding percentage of abundance in each cluster.
        Cluster    most abundance genus    % abundance of most abundance genus  2nd most ...
    '''
    def describe_taxa_level(self, clusters, blast, uniprotinfo, output, 
                            taxa_level = 'genus', deepness = 5):
        clusters = pd.read_csv(clusters, sep = '\t')
        column = 'Taxonomic lineage (' + taxa_level.upper() + ')'
        blast = mtools.parse_blast(blast)
        print('Organizing blast information.')
        pbar = ProgressBar()
        blast_part = pd.DataFrame([['_'.join(blast.iloc[i]['qseqid'].split('_')[:6]), 
                      float(blast.iloc[i]['qseqid'].split('_')[5]),
                      blast.iloc[i]['sseqid'].split('|')[1] if blast.iloc[i]['sseqid'] != '*'
                      else blast.iloc[i]['sseqid']] for i in pbar(range(len(blast)))])
        blast_part.columns = ['contig','coverage','Entry']
        joined = pd.merge(clusters, blast_part, on='contig')
        uniprotinfo = pd.read_csv(uniprotinfo, sep = '\t').drop_duplicates()
        joined = pd.merge(joined, uniprotinfo, on = 'Entry')
        handler = open(output, 'w')
        pbar = ProgressBar()
        for cluster in pbar(list(set(joined['cluster']))):
            joined_partial = joined[joined['cluster'] == cluster]
            joined_partial = joined_partial.groupby(['cluster', column])['coverage'].sum().reset_index()
            joined_partial = joined_partial.sort_values('coverage', ascending = False)
            joined_partial['coverage'] /= joined_partial['coverage'].sum()
            handler.write(str(cluster) + '\t' + '\t'.join([joined_partial.iloc[i][column] + '\t' + 
            str(joined_partial.iloc[i]['coverage']) for i in range(len(joined_partial))]) + '\n')
        handler.close()
        
    '''
    Input: 
        data: EXCEL file outputed by Binning.calculate_epsilon
        taxa_level: The taxa level to be described - superkingdom, phylum, class, 
        order, family, genus or species
        best_clusters: TSV file with best clusters outputed by Binning.calculate_epsilon
        points: TXT file outputed by VizBin with contigs coordinates after binning
        output: name of plot file to output, should terminate in .png
    Output:
        
    '''
    def plot_clusters(self, data, taxa_level, best_clusters, points, output, 
                      label = True, subtitle_size = 20):
        points = pd.read_csv(points, header = None)
        points.columns = ['lat','lon']
        best_clusters = pd.read_csv(best_clusters, sep = '\t')
        points = pd.concat([points, best_clusters], axis = 1)
        numeration = {'superkingdom':'1', 'phylum':'2', 'class':'3', 'order':'4',
                      'family':'5', 'genus':'6', 'species':'7'}
        column = numeration[taxa_level] + '.Taxonomic lineage (' + taxa_level.upper() + ')'
        data = pd.read_excel(data, index_col = [0,1])
        partial = data.xs(column, level = 1)
        points = pd.merge(points, partial, left_on = 'cluster', right_index = True)
        points = points[['lat','lon','Taxa']]
        points = points.fillna(value = 'Not identified')
        
        taxa = list(set(points['Taxa']))
        colors = iter(cm.rainbow(np.linspace(0, 1, len(taxa))))
        plt.gcf().clear()
        for i in range(len(taxa)):
            partial_points = points[points['Taxa'] == taxa[i]]
            plt.scatter(partial_points['lat'], partial_points['lon'], 0.1,
                        color = next(colors), label = taxa[i], marker = 'o')
        if label: 
            label = plt.legend(loc='best')
            for i in range(len(taxa)):
                label.legendHandles[i]._sizes = [subtitle_size]
        plt.savefig(output, bbox_inches='tight')
        
    def run(self):
        self.run_vizbin_binning(self.contigs, self.output)
        self.calculate_epsilon(self.contigs, self.output + '/points.txt', 
                               self.blast, self.uniprotinfo, self.output)
        
    '''
    *** The following functions concern binning with MaxBin2 ***
    '''
    
    '''
    Input:
        contigs: FASTA file with contigs
        output: basename of output
        threads: number of threads to use by Maxbin
        mg1: name of forward reads file used in the assembly
        mg2: name of reverse reads file used in the assembly
        abundance: name of abundance file (format is contig\tabundance)
        marketset: either '107' marker genes present in >95% of bacteria, or
        '40' marker gene sets that are universal among bacteria and archaea. 
        '40' may be better suited for environment dominated by archaea; 
        however it tend to split genomes into more bins.
    Output:
        bins named basename + .n.fasta
        abundance of each contig in each sample (mg1 and mg2) named basename + .abund1/2
        abundace of each bin for both samples named basename + .abundance
        log of workflow named basename + .log
        markergenes used to compose the bins named basename + .marker
        contigs not included in any bin named basename + .noclass
        
    '''
            
    def run_maxbin(self, contigs, output, threads = 8, mg1 = None, mg2 = None,
                   abundance = None, markerset = '107'):
        bashCommand = 'run_MaxBin.pl -contig ' + contigs + ' -out ' + output
        parameter_dictionary = {'mg1':'-reads','mg2':'-reads2','abundance':'-abundance'}
        for parameter in [mg1, mg2, abundance]:
            if parameter is not None:
                bashCommand += ' ' + parameter_dictionary[parameter] + ' ' + parameter
        bashCommand += ' -thread ' + threads + ' -markerset ' + markerset
        mtools.run_command(bashCommand)
        
        
        
if __name__ == '__main__':
    
    binner = Binner(contigs = 'Binning/all_info_contigs.fasta', 
                    output = 'SimulatedMGMT/Binning/genus_taxa_description.tsv',
                    blast = 'SimulatedMGMT/Annotation/aligned.blast',
                    uniprotinfo = 'SimulatedMGMT/Annotation/uniprot.info')
    
    taxa_list = ['superkingdom','phylum','class','order','family','genus','species']
    
    pbar = ProgressBar()
    
    for taxon in pbar(taxa_list):
        binner.plot_clusters('SimulatedMGMT/Binning/grinder-reads/binning_results.xlsx',
                             taxon, 'SimulatedMGMT/Binning/grinder-reads/best_clusters.tsv',
                             'SimulatedMGMT/Binning/grinder-reads/points.txt',
                             'SimulatedMGMT/Binning/grinder-reads/' + taxon + '.png',
                             label = True)
    
    
    
    '''
    clusters = binner.cluster_coords('Binning/VizBin/' + sample + '/over_cutoff_contigs.fasta',
                          'Binning/VizBin/' + sample + '/points.txt')
    
    clusters.to_csv('Binning/VizBin/' + sample + '/clusters.tsv', sep = '\t')
    '''

#MaxBin2 must be used with run_MaxBin.pl