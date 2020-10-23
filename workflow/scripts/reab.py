# -*- coding: utf-8 -*-
"""
MOSCA's ReAssembly with Binning package for iterative improvement of assembly

By João Sequeira

Feb 2020
"""

# must add 
'''
conda install -c bioconda samtools
conda install -c bioconda bedtools
'''

from MOSCA.scripts.mosca_tools import MoscaTools
import glob, pathlib, shutil, os, pandas as pd, numpy as np

mtools = MoscaTools()

class reAB:
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs

    def kbase_reab(self): #assembly, bins, output_directory, reads1, reads2 = None, threads = '15'):
        '''
        # first, reads had to be fixed. There were many @1000000/1 because the files were merged when first simulated.
        awk 'BEGIN{i=1}{split($0,a,"/");print "@"i++"/"a[2];getline;print $0;getline;print $0;getline;print $0}' SimulatedMGMT/Assembly/Sample_forward.fastq > SimulatedMGMT/Assembly/Sample_forward_fixed.fastq
        awk 'BEGIN{i=1}{split($0,a,"/");print "@"i++"/"a[2];getline;print $0;getline;print $0;getline;print $0}' SimulatedMGMT/Assembly/Sample_reverse.fastq > SimulatedMGMT/Assembly/Sample_reverse_fixed.fastq
        mv SimulatedMGMT/Assembly/Sample_forward_fixed.fastq SimulatedMGMT/Assembly/Sample_forward.fastq
        mv SimulatedMGMT/Assembly/Sample_reverse_fixed.fastq SimulatedMGMT/Assembly/Sample_reverse.fastq
        
        # alignment had to be repeated because of the mistake in reads names. Not assembly, because reads were all the same.
        bowtie2 -x SimulatedMGMT/Assembly/Sample/contigs_index -1 SimulatedMGMT/Assembly/Sample_forward.fastq -2 SimulatedMGMT/Assembly/Sample_reverse.fastq -S test_simulated_reab/alignment.sam -p 15 1> test_simulated_reab/alignment.report 2> test_simulated_reab/alignment.log
        
        # SAM is converted to sorted and indexed BAM for splicing by bin
        samtools view -S -@ 15 -b test_simulated_reab/alignment.sam > test_simulated_reab/alignment.bam
        samtools sort test_simulated_reab/alignment.bam -o test_simulated_reab/alignment_sorted.bam -@ 15
        samtools index test_simulated_reab/alignment_sorted.bam
        
        # Reads used in the assembly are joined
        awk '{printf substr($0,2);getline;printf"\t"$0;getline;getline;print "\t"$0}' MOSCAfinal/Assembly/Sample_forward.fastq > test_reab/realMGMT/reads1.txt
        awk '{printf substr($0,2);getline;printf"\t"$0;getline;getline;print "\t"$0}' MOSCAfinal/Assembly/Sample_reverse.fastq > test_reab/realMGMT/reads2.txt
        join test_reab/realMGMT/reads1.txt test_reab/realMGMT/reads2.txt | sort > test_reab/realMGMT/joined_reads.txt
        '''
        threads = '15'
        bin_files = glob.glob('MOSCAfinal/Binning/Sample/*.fasta')
        
        beans_eaten = [name.split('/')[2] for name in glob.glob('test_reab/realMGMT/*/')]
        beans = list()
        
        for file in bin_files:
            bean = file.split('/Sample.')[-1].split('.fasta')[0]                # has to be reworked for when working with several different communities
            beans.append(bean)
            
            if bean not in beans_eaten:
                # for each bin, its alignment, through the contigs composing it.
                mtools.run_pipe_command("grep '>' {} | samtools view -b -@ 15 test_reab/realMGMT/alignment.bam $(awk '{{print substr($0,2)}}')".format(file), file = 'test_reab/realMGMT/bin.bam')
                mtools.run_command("samtools view -@ 15 -h -o test_reab/realMGMT/bin.sam test_reab/realMGMT/bin.bam")
                
                
                # filter the file with all the reads for the reads of the bin
                mtools.run_pipe_command("grep -v '@' test_reab/realMGMT/bin.sam | awk '{{print $0}}' | sort | uniq | join test_reab/realMGMT/joined_reads.txt -", file = 'test_reab/realMGMT/reads_bin.txt')
                mtools.run_pipe_command("awk '{{print \"@\"$1\" \"$2\"\\n\"$3\"\\n+\\n\"$4 > \"test_reab/realMGMT/bin{0}_forward.fastq\";print \"@\"$1\" \"$5\"\\n\"$6\"\\n+\\n\"$7 > \"test_reab/realMGMT/bin{0}_reverse.fastq\"}}' test_reab/realMGMT/reads_bin.txt".format(bean))
                
                try:
                    # re-assembly the reads by bin
                    mtools.run_command('metaspades.py -1 test_reab/realMGMT/bin{0}_forward.fastq -2 test_reab/realMGMT/bin{0}_reverse.fastq -o test_reab/realMGMT/bin{0} -t {1}'.format(bean, threads))
                    #shutil.copy('test_reab/realMGMT/bin{0}/contigs.fasta'.format(bean), 'test_reab/realMGMT/second_bins/bin{0}.fasta'.format(bean))
                    
                    binning = False
                    try:
                        pathlib.Path('test_reab/realMGMT/bin{}/maxbin'.format(bean)).mkdir(parents=True, exist_ok=True)
                        mtools.run_command('run_MaxBin.pl -contig test_reab/realMGMT/bin{0}/contigs.fasta -out test_reab/realMGMT/bin{0}/maxbin/bin.{0} -thread {1} -markerset 40 -reads test_reab/realMGMT/bin{0}_forward.fastq -reads2 test_reab/realMGMT/bin{0}_reverse.fastq'.format(bean, threads))
                        binning = True
                    except:
                        print('Binning has failed! Likely because of lack of marker genes in contigs (avg < 1).')
                        
                    if binning:
                        pathlib.Path('test_reab/realMGMT/bin{}/checkm'.format(bean)).mkdir(parents=True, exist_ok=True)
                        mtools.run_command('checkm lineage_wf -x fasta -r --ali --nt -t {0} --pplacer_threads {0} test_reab/realMGMT/bin{1}/maxbin test_reab/realMGMT/bin{1}/checkm --tab_table --file test_reab/realMGMT/bin{1}/checkm/output_table.tab'.format(threads, bean))
                except:
                    print('Assembly failed! Dunno why! :O')
        #mtools.run_command('checkm lineage_wf -x fasta -r --ali --nt -t {0} --pplacer_threads {0} test_simulated_reab/second_bins test_simulated_reab/second_bins/checkm --tab_table --file test_simulated_reab/second_bins/chekm/output_table.tab'.format(threads))
        '''
        soup = pd.DataFrame()
        for bean in beans:
            print(bean)
            partial = pd.read_csv('test_simulated_reab/bin{}/checkm/output_table.tab'.format(bean), sep = '\t')
            partial.index = [bean]*len(partial)
            soup = pd.concat([soup, partial])
        soup.to_excel('test_simulated_reab/second_checkm.xlsx')
        '''
    def confirm_reads_origin(self):
        for file in glob.glob('SimulatedMGMT/Binning/Sample/*.fasta'):
            bean = file.split('/Sample.')[-1].split('.fasta')[0]
            mtools.run_pipe_command("cat SimulatedMGMT/Binning/Sample/Sample.{0}.fasta | grep '>' | samtools view -b -@ 15 test_simulated_reab/alignment_sorted.bam $(awk '{{print substr($0,2)}}') | samtools view -@ 15 -h -S - | grep -v '@' | awk '{{split($1, a, \"/\"); print a[1]}}' | sort -f | uniq | join - test_simulated_reab/id2taxid.tsv | awk '{{print $2}}' | sort | uniq -c > test_simulated_reab/bin{0}.counts".format(bean))
        result = pd.DataFrame(columns = ['ASM1344v1','ASM20441v1','gtlEnvA5udCFS',
            'DSM1535','ASM1496v1','ASM1340v1','Mb_MS2','ASM96103v1','ASM1504v1',
            'ASM1548v1','ASM73042v1','ASM148326v1','ASM131834v1','ASM2420v1',
            'ASM1056v1','ASM1472v1','Syntrophomonas_Zehnderi_OL-4','spAn4DRAFT_v1',
            'BacMEV2011_1.0','ASM83310v2','ASM38963v1','ASM162387v1','ASM16579v1'])
        counts_files = glob.glob('test_simulated_reab/bin*.counts')
        
        for file in counts_files:
            partial = pd.read_csv(file,sep='\s',header=None,index_col=1).transpose()
            partial.index = [file.split('test_simulated_reab/')[1].split('.counts')[0]]
            result = pd.concat([result, partial])
        result.fillna(value=0,inplace=True)
        result = result[['ASM1344v1','ASM20441v1','gtlEnvA5udCFS','DSM1535','ASM1496v1','ASM1340v1','Mb_MS2','ASM96103v1','ASM1504v1','ASM1548v1','ASM73042v1','ASM148326v1','ASM131834v1','ASM2420v1','ASM1056v1','ASM1472v1','Syntrophomonas_Zehnderi_OL-4','spAn4DRAFT_v1','BacMEV2011_1.0','ASM83310v2','ASM38963v1','ASM162387v1','ASM16579v1']]
        result.sort_index(inplace=True)
        
    def iterative_binning(self, minimum = 0.4, maximum = 1, step = 0.05):
        
        directory = 'test_reab/simulated/2nd_first_bins'
        forward = 'MOSCAfinal/Assembly/Sample_forward.fastq'
        reverse = 'MOSCAfinal/Assembly/Sample_reverse.fastq'
        contigs = 'MOSCAfinal/Assembly/Sample/contigs.fasta'
        h = open('{}/checkm_assessment.txt'.format(directory), 'w')
        h.write('PT\tHigh\tMedium\tLow\n')
        h.close()
        
        pts_failed = list()
        
        for pt in np.arange(minimum, maximum + step, step):
            try:
                pathlib.Path('{}/{}/checkm'.format(directory, str(pt))).mkdir(parents=True, exist_ok=True)
                mtools.run_command('run_MaxBin.pl -contig {4} -out {1}/{0}/{0} -thread 15 -markerset 40 -reads {2} -reads2 {3} -prob_threshold {0}'.format(str(pt), directory, forward, reverse, contigs), print_message = False, verbose = False)
                mtools.run_command('checkm lineage_wf -x fasta -r --ali --nt -t 15 --pplacer 15 {1}/{0} {1}/{0}/checkm --tab_table --file {1}/{0}_output_table.tab'.format(str(pt), directory), print_message = False, verbose = False)
                
                data = pd.read_csv('{}/{}_output_table.tab'.format(directory, str(pt)),sep='\t')
                high = ((data['Completeness'].astype(float) >= 90) & (data['Contamination'].astype(float) < 5)).sum()
                medium = ((data['Completeness'].astype(float) >= 50) & (data['Completeness'].astype(float) < 90) & (data['Contamination'].astype(float) < 10)).sum()
                low = ((data['Completeness'].astype(float) < 50) & (data['Contamination'].astype(float) < 10)).sum()
                
                with open('{}/checkm_assessment.txt'.format(directory), 'a') as f:
                    f.write('{}\t{}\t{}\t{}\n'.format(str(pt), str(high), str(medium), str(low)))
            except:
                pts_failed.append(pt)
        print('PTs failed: ' + ','.join(pts_failed))

    def albertsen_reab(self):
        #conda install -c bioconda cytoscape
        pass

'''  
# SimulatedMGMT

# from "work" directory

# first, reads had to be fixed. There were many @1000000/1 because the files were merged when first simulated.
awk 'BEGIN{i=1}{split($0,a,"/");print "@"i++"/"a[2];getline;print $0;getline;print $0;getline;print $0}' SimulatedMGMT/Assembly/Sample_forward.fastq > SimulatedMGMT/Assembly/Sample_forward_fixed.fastq
awk 'BEGIN{i=1}{split($0,a,"/");print "@"i++"/"a[2];getline;print $0;getline;print $0;getline;print $0}' SimulatedMGMT/Assembly/Sample_reverse.fastq > SimulatedMGMT/Assembly/Sample_reverse_fixed.fastq
mv SimulatedMGMT/Assembly/Sample_forward_fixed.fastq SimulatedMGMT/Assembly/Sample_forward.fastq
mv SimulatedMGMT/Assembly/Sample_reverse_fixed.fastq SimulatedMGMT/Assembly/Sample_reverse.fastq

# alignment had to be repeated. Not assembly, because reads were all the same.
bowtie2 -x SimulatedMGMT/Assembly/Sample/contigs_index -1 SimulatedMGMT/Assembly/Sample_forward.fastq -2 SimulatedMGMT/Assembly/Sample_reverse.fastq -S test_simulated_reab/alignment.sam -p 15 1> test_simulated_reab/alignment.report 2> test_simulated_reab/alignment.log

# SAM is converted to sorted and indexed BAM for splicing by bin
samtools view -S -@ 15 -b test_simulated_reab/alignment.sam > test_simulated_reab/alignment.bam
samtools sort test_reab/alignment.bam -o test_reab/alignment_sorted.bam -@ 15
samtools index test_reab/alignment_sorted.bam

# get the numbers of the reads in specific bin -> I can do this way because I had only one pair of MG files, otherwise would have to combine file/sample information in the reads headers
awk '{split($0,a,"/"); printf substr(a[1],2,length);getline;printf"\t"$0;getline;getline;print "\t"$0}' SimulatedMGMT/Assembly/Sample_forward.fastq > test_reab/reads1.txt
awk '{split($0,a,"/"); printf substr(a[1],2,length);getline;printf"\t"$0;getline;getline;print "\t"$0}' SimulatedMGMT/Assembly/Sample_reverse.fastq > test_reab/reads2.txt
join -i test_reab/reads1.txt test_simulated_reab/reads2.txt > test_reab/joined_reads.txt
sort -f test_reab/joined_reads.txt > test_reab/joined_reads_sorted.txt

# for each bin, its alignment, through the contigs composing it.
cat SimulatedMGMT/Binning/Sample/Sample.001.fasta | grep '>' | samtools view -b -@ 15 test_reab/alignment_sorted.bam $(awk '{print substr($0,2)}') > test_reab/bin1.bam
samtools view -@ 15 -h -o test_reab/bin1.sam test_reab/bin1.bam

# get the reads' IDs of the bin
grep -v '@' test_reab/bin1.sam | awk '{split($1, a, "/"); print a[1]}' > test_reab/bin1_ids.txt
sort -f test_reab/bin1_ids.txt > test_reab/bin1_ids_sorted.txt
uniq test_reab/bin1_ids_sorted.txt > test_reab/bin1_ids_sorted_uniq.txt # not needed? Got 4305321 unique ids / 8534786 ids = 50.4% -> almost always mate pair alignment

# filter the file with all the reads for the reads of the bin
join -i test_reab/bin1_ids_sorted_uniq.txt test_reab/joined_reads_sorted.txt > test_reab/reads_bin1.txt
awk '{print "@"$1"/1\n"$2"\n+\n"$3 > "test_reab/bin1_forward.fastq";print "@"$1"/2\n"$4"\n+\n"$5 > "test_reab/bin1_reverse.fastq"}' test_reab/reads_bin1.txt

# associate read id to species id
awk '{split($1,a,"/");split($2,b,":");print substr(a[1],2,length)"\t"b[3];getline;getline;getline}' SimulatedMGMT/Assembly/Sample_forward.fastq | sort > test_reab/id2taxid.tsv

cat SimulatedMGMT/Binning/Sample/Sample.002.fasta | grep '>' | samtools view -b -@ 15 test_reab/alignment_sorted.bam $(awk '{print substr($0,2)}') | samtools view -@ 15 -h -S - | grep -v '@' | awk '{split($1, a, "/"); print a[1]}' | sort -f | uniq | join - test_simulated_reab/id2taxid.tsv | awk '{print $2}' | sort | uniq -c > test2.counts

checkm lineage_wf -x fasta -r --ali --nt -t 15 --pplacer 15 test_reab/second_bins/ test_reab/second_bins/checkm --tab_table --file test_reab/second_bins/checkm/output_table.tab
'''
'''

# MOSCAfinal
# SAM is converted to sorted and indexed BAM for splicing by bin
samtools view -S -@ 15 -b MOSCAfinal/Assembly/Sample/quality_control/alignment.sam | samtools sort -@ 15 > test_reab/realMGMT/alignment.bam
samtools index test_reab/realMGMT/alignment.bam

# get the numbers of the reads in specific bin -> I can do this way because I had only one pair of MG files, otherwise would have to combine file/sample information in the reads headers -> maybe not, if real data, and not simulated
awk '{split($0,a,"/"); printf substr(a[1],2,length);getline;printf"\t"$0;getline;getline;print "\t"$0}' SimulatedMGMT/Assembly/Sample_forward.fastq > test_reab/realMGMT/reads1.txt
awk '{split($0,a,"/"); printf substr(a[1],2,length);getline;printf"\t"$0;getline;getline;print "\t"$0}' SimulatedMGMT/Assembly/Sample_reverse.fastq > test_reab/realMGMT/reads2.txt
join -i test_reab/realMGMT/reads1.txt test_simulated_reab/reads2.txt > test_reab/realMGMT/joined_reads.txt
sort -f test_reab/realMGMT/joined_reads.txt > test_reab/realMGMT/joined_reads_sorted.txt

# for each bin, its alignment, through the contigs composing it.
cat SimulatedMGMT/Binning/Sample/Sample.001.fasta | grep '>' | samtools view -b -@ 15 test_reab/realMGMT/alignment.bam $(awk '{print substr($0,2)}') > test_reab/realMGMT/bin1.bam
samtools view -@ 15 -h -o test_reab/realMGMT/bin1.sam test_reab/realMGMT/bin1.bam

# get the reads' IDs of the bin
grep -v '@' test_reab/realMGMT/bin1.sam | awk '{split($1, a, "/"); print a[1]}' > test_reab/realMGMT/bin1_ids.txt
sort -f test_reab/realMGMT/bin1_ids.txt > test_reab/realMGMT/bin1_ids_sorted.txt
uniq test_reab/realMGMT/bin1_ids_sorted.txt > test_reab/realMGMT/bin1_ids_sorted_uniq.txt # not needed? Got 4305321 unique ids / 8534786 ids = 50.4% -> almost always mate pair alignment

# filter the file with all the reads for the reads of the bin
join -i test_reab/realMGMT/bin1_ids_sorted_uniq.txt test_reab/realMGMT/joined_reads_sorted.txt > test_reab/realMGMT/reads_bin1.txt
awk '{print "@"$1"/1\n"$2"\n+\n"$3 > "test_reab/realMGMT/bin1_forward.fastq";print "@"$1"/2\n"$4"\n+\n"$5 > "test_reab/realMGMT/bin1_reverse.fastq"}' test_reab/realMGMT/reads_bin1.txt

# associate read id to species id
awk '{split($1,a,"/");split($2,b,":");print substr(a[1],2,length)"\t"b[3];getline;getline;getline}' SimulatedMGMT/Assembly/Sample_forward.fastq | sort > test_reab/realMGMT/id2taxid.tsv

cat SimulatedMGMT/Binning/Sample/Sample.002.fasta | grep '>' | samtools view -b -@ 15 test_reab/realMGMT/alignment.bam $(awk '{print substr($0,2)}') | samtools view -@ 15 -h -S - | grep -v '@' | awk '{split($1, a, "/"); print a[1]}' | sort -f | uniq| join - test_simulated_reab/id2taxid.tsv | awk '{print $2}' | sort | uniq -c > test2.counts

checkm lineage_wf -x fasta -r --ali --nt -t 15 --pplacer 15 test_reab/realMGMT/second_bins/ test_reab/realMGMT/second_bins/checkm --tab_table --file test_reab/realMGMT/second_bins/checkm/output_table.tab
'''


if __name__ == '__main__':
    
    reab = reAB()
    
    reab.iterative_binning()
    
    #reab.confirm_reads_origin()

