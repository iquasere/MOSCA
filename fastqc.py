# -*- coding: utf-8 -*-
'''
Trimmomatic python API

By João Sequeira

7th March 2017
'''

class FastQC:
    def __init__ (self, **kwargs):
        self.attributes = kwargs
		
    def bash_command(self):
        result = 'fastqc' 
        if 'outdir' in self.attributes.keys():                          #directory of results
            result += ' --outdir ' + self.attributes['outdir']
        if 'casava' in self.attributes.keys():                          #if files come from casava output, value is list of files
            result += ' --casava '
            for file in self.attributes['casava']:
                result += file + ' '
            result = result.rstrip(' ')
        if self.attributes['extract'] == True:                          #zipped output file will be uncompressed
            result += ' --extract'
        else:                                                           #or not
            result += ' --noextract'            
        if 'java' in self.attributes.keys():                            #full path to java binary, PATH if nokey
            result += ' --java ' + self.attributes['java']
        if 'nogroup' == True:                                           #disable grouping of bases for reads >50bp
            result += ' --nogroup'
        if 'format' in self.attributes.keys():                          #valid options are bam,sam,bam_mapped,sam_mapped and fastq
            result += ' --format' + self.attributes['format']           
        if 'threads' in self.attributes.keys():                         #number of files processed simultaneously (each thread is allocated 250mb of memory to)
            result += ' --threads' + str(self.attributes['threads'])
        if 'contaminants' in self.attributes.keys():                    #file in format name[tab]sequence of contaminants
            result += ' --contaminants' + self.attributes['contaminants']
        if 'kmers' in self.attributes.keys():                           #length of kmer to look in kmer content (between 2 and 10, default is 5)
            result += ' --kmers' + self.attributes['kmers']
        for file in self.attributes['files']:
            result += ' ' + file
        print(result)
        return result

    def run(self):
        import subprocess
        bashCommand = self.bash_command()
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error