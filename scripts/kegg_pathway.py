#!/usr/bin/env python

from mosca_tools import MoscaTools
from Bio.KEGG.REST import kegg_get, kegg_link, kegg_list
from Bio.KEGG.KGML import KGML_parser, KGML_pathway
from Bio.Graphics.KGML_vis import KGMLCanvas
from matplotlib import cm
from matplotlib.colors import to_hex
from io import StringIO
import matplotlib.pyplot as plt
import re, pandas as pd, numpy as np, PIL, os
from progressbar import ProgressBar

mtools = MoscaTools()

__author__ = "Tiago Oliveira"
__credits__ = ["Tiago Oliveira", "Joao Sequeira"]
__version__ = "1.0"
__maintainer__ = "Joao Sequeira"
__email__ = "maildosequeira@gmail.com"
__status__ = "Production"

################################################################################
class KEGGPathway:
    '''
    This class concerns with the storage of information manually retrieved from 
    KEGG documentation for reference, and the conversion of certain IDs to 
    different IDs from the KEGG DB
    '''
    
    def __init__(self, **kwargs):
        '''
        Initialize object
        :param **kwargs - none is necessary for this class
        '''
        
        self.__dict__ = kwargs
        
        self.maps = dict()
        for line in open('MOSCA/Databases/kegg_pathway/metabolic_maps.txt'
                         ).read().split('\n'):
            if not line.startswith('#'):
                parts = line.split('\t')
                self.maps[parts[0]] = parts[1]
    
        self.default_maps = open('MOSCA/Databases/kegg_pathway/default_maps.txt'
                                 ).read().split('\n')
        
    def keggid2ko(self, kegg_ids, step = 150):
        '''
        Converts KEGG_ID genes to Ortholog KO ID from KEGG
        :param KEGG_ID: (list) - KEGG ID genes
        :param step: (int) - will convert "step" KEGG IDs at a time
        :return: (list) - (list,list) - KEGG ID genes converted and ko IDs
        '''
        print('Converting {:d} KEGG IDs to KOs through the KEGG API.'.format(len(kegg_ids)))
        result = list()
        pbar = ProgressBar()
        for i in pbar(range(0, len(kegg_ids) - step, step)):
            try:
                result += kegg_link("ko", kegg_ids[i:i+step]).read().split("\n")[:-1]
            except:
                print('KEGG ID to KO broke at index ' + str(i))
                result = [[part[0] + ';', part[1].strip('ko:')] for part in
                   [relation.split('\t') for relation in result]]
                return pd.DataFrame(result, columns = ['Cross-reference (KEGG)', 'KO (KEGG Pathway)'])
        result += kegg_link("ko", kegg_ids[len(kegg_ids) - step:]).read().split("\n")[:-1]
        result = [[part[0] + ';', part[1].strip('ko:')] for part in
                   [relation.split('\t') for relation in result]]
        return pd.DataFrame(result, columns = ['Cross-reference (KEGG)', 'KO (KEGG Pathway)'])

    def ko2ec(self, kos, step = 150):
        '''
        Converts KOs to EC numbers
        :param kos: list of kegg orthologs
        :return: dic associating ortholog kegg id with list
        of assotiated EC numbers
        '''
        print('Retrieving EC numbers from {} KOs.'.format(len(kos)))
        result = list()
        pbar = ProgressBar()
        for i in pbar(range(0, len(kos), step)):
            try:
                result += kegg_link("enzyme", kos[i:i+step]).read().split("\n")[:-1]
            except:
                print('KO to EC number broke at index ' + str(i))
                result = [relation.split('\t') for relation in result]
                return list(map(list, zip(*result)))
        result += kegg_link("enzyme", kos[len(kos) - step:]).read().split("\n")[:-1]
        result = [[part[0].strip('ko:'),part[1].upper()] for part in
                   [relation.split('\t') for relation in result]]
        return pd.DataFrame(result, columns = ['KO (KEGG Pathway)', 'EC number (KEGG Pathway)'])
    
    def most_abundant_taxa(self, data, samples, number_of_taxa = 10, 
                           level_of_taxa = 'GENUS'):
        '''
        Calculates top genus from samples
        :param samples: list of samples to consider for quantification of genus abundance
        :param n: number of top genus to return
        :return: list of top genus
        '''
        data = data.groupby("Taxonomic lineage (" + level_of_taxa + ")")[samples].sum()
        data["sums"] = data.sum(axis=1)
        data = data.sort_values(by=["sums"], ascending=False)
        if number_of_taxa > len(data.index.tolist()):
            number_of_taxa = len(data.index.tolist())
        return data.index.tolist()[:number_of_taxa]
    
    def kegg_maps_available(self):
        '''
        Creates a dic with all specific kegg pathway maps and their description
        :return: pandas.DataFrame with Map ID as index and maps names as
        sole column
        '''
        maps = pd.read_csv(StringIO(kegg_list("pathway").read()), sep='\t',
                           names = ['ID', 'Name'])
        
        maps['ID'] = maps['ID'].apply(lambda x: x.split(':')[1])
        maps.set_index('ID', inplace = True)
        
        return maps
    
    def taxa_colors(self, colors = None, ncolor = 1):
        '''
        Creates list of hex colors to be used, using matplotlib or using custom colors
        :param colors: list of hex colors
        :param ncolor: int indicating the ammount of hex colors should be created
        :return: returns list with hex color codes
        '''
        if not colors:                                                          # if no colors are given creates a list of hex colors with ncolor from matplotlib discrete colormaps
            color_scheme = (cm.get_cmap('Pastel2', 8) if ncolor <= 8
                            else cm.get_cmap("Set3", 12) if ncolor <= 12
                            else cm.get_cmap("rainbow", ncolor))                # if ncolor > 12 a continuous colormap is used instead
            return [to_hex(color_scheme(i)) for i in range(ncolor)]
        else:                                                                   # validates hex values and returns the original list
            isvalidhex = True
            for hexvalue in colors:
                if not re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', hexvalue):
                    isvalidhex = False
            if isvalidhex:
                return colors
            else:
                raise Exception("Colors aren't valid hex codes")
                
    def create_potential_legend(self, colors, labels, filename, resize_factor = 10):
        '''
        Draws the color to taxa labels of genomic potential representations
        :param colors: list - list of colors of the different taxa
        :param labels: list - list of taxa corresponding to the colors
        :param filename: string - filename to output
        :param size: int - how big you want your legend?
        '''
        f = lambda m, c: plt.plot([], [], marker = m, color = c, ls = "none")[0]
        handles = [f("s", color) for color in colors]
        legend = plt.legend(handles, labels, loc=3, framealpha=1, frameon = True)
        fig = legend.figure
        fig.canvas.draw()
        # The bbox manipulation removes the axis
        bbox = legend.get_window_extent()
        bbox = bbox.from_extents(*(bbox.extents + np.array([-2,-2,2,2])))
        bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
        
        fig.savefig(filename, dpi="figure", bbox_inches=bbox)

    def add_blank_space(self, image_pil, width, height, image_mode = 'RGBA'):
        '''
        Resizes an image with white background, keeping image size ratio
        :param image_pil: PIL.Image - image to be resized
        :param width: int - width of final image
        :param height: int - heigth of final image
        :param image_mode: str - image mode of image (RGBA, RGB, ...)
        '''
        ratio_w = width / image_pil.width
        ratio_h = height / image_pil.height
        if ratio_w < ratio_h:
            # It must be fixed by width
            resize_width = width
            resize_height = round(ratio_w * image_pil.height)
        else:
            # Fixed by height
            resize_width = round(ratio_h * image_pil.width)
            resize_height = height
        image_resize = image_pil.resize((resize_width, resize_height), 
                                        PIL.Image.ANTIALIAS)
        background = PIL.Image.new('RGBA', (width, height), (255, 255, 255, 255))
        offset = (round((width - resize_width) / 2), 
                  round((height - resize_height) / 2))
        background.paste(image_resize, offset)
        return background.convert(image_mode)
    
    def resize_image(self, image_pil, ratio = None, width = None, height = None):
        '''
        Resizes an image with alteration to image size ratio
        :param ratio: int - ratio of resize - how bigger or smaller will the output be?
        :param image_pil: PIL.Image - image to be resized
        :param width: int - width of final image
        :param height: int - heigth of final image
        '''
        if ratio:
            return image_pil.resize((image_pil.width * ratio, 
                                     image_pil.height * ratio), PIL.Image.ANTIALIAS)
        elif width and height:
            return image_pil.resize((width, height), PIL.Image.ANTIALIAS)
        else:
            return None
        
    def pdf2png(self, pdf_filename):
        '''
        Converts a pdf file to a png file, RGB format - name changes, .pdf to .png
        :param pdf_filename: str - filename of PDF file
        '''
        mtools.run_command('pdftoppm {} {} -png'.format(pdf_filename, 
                           pdf_filename.split('.pdf')[0]), print_message = False)
        os.rename(pdf_filename.replace('.pdf', '-1.png'), 
                  pdf_filename.replace('.pdf', '.png'))
    
    def add_legend(self, kegg_map_file, legend_file, output):
        '''
        Merges the two files - KEGG metabolic map and respective legend - into
        one file file
        :param kegg_map_file: str - filename of PDF kegg metabolic map
        :param legend_file: str - filename of PNG legend
        '''
        #self.pdf2png(kegg_map_file)
        imgs = [PIL.Image.open(file) for file in 
                [kegg_map_file.replace('.pdf', '.png'), legend_file]]
        imgs[0] = imgs[0].convert('RGBA')                                       # KEGG Maps are converted to RGB by pdftoppm, dunno if converting to RGBA adds any transparency
        imgs[1] = self.resize_image(imgs[1], ratio = 5)
        imgs[1] = self.add_blank_space(imgs[1], imgs[1].width, imgs[0].height)
        imgs_comb = np.hstack([np.asarray(i) for i in imgs])
        
        # save that beautiful picture
        imgs_comb = PIL.Image.fromarray(imgs_comb)
        imgs_comb.save(output)

    def genomic_potential_taxa(self, data, samples, genera = None, 
                               number_of_taxa = 10, level_of_taxa = 'GENUS', 
                               metabolic_map = None, output_basename = None, 
                               maxshared = 10):
        '''
        Represents the genomic potential of the dataset for a certain taxa level,
        by coloring each taxon with a unique color
        :param data: pandas.DataFrame with data already processed by KEGGPathway
        :param samples: list of str column names of the dataset correspoding to
        expression values
        :param genera: list of genus to represent
        :param number_of_taxa: int representing the number of diferent taxa to 
        be represented in the maps, in case the taxa are not specified (will always
        be used in the common MOSCa pipeline)
        :param level_of_taxa: str - taxonomic level to represent - SPECIES,
        SUPERKINGDOM, ...
        :param output_basename: str - basename for map outputs
        :param maxshared: int - maximum number of different taxa to represent
        in a single map box
        '''
        if genera is None:
            genera = self.most_abundant_taxa(data, samples, number_of_taxa)
            
        colors = self.taxa_colors(ncolor = len(genera))
        dic_colors = {genera[i] : colors[i] for i in range(len(genera))}
        '''
        pathway = KeggMap(data, metabolic_map)
        taxa_in_box = dict()
        for genus in genera:
            df = data[data["Taxonomic lineage (" + level_of_taxa + ")"] == genus][samples + ['KO (KEGG Pathway)']]
            df = df[df.any(axis=1)]
            for ortholog in df['KO (KEGG Pathway)']:
                if ortholog in pathway.ko_boxes.keys():
                    for box in pathway.ko_boxes[ortholog]:
                        if box in taxa_in_box.keys():
                            if genus not in taxa_in_box[box]:
                                taxa_in_box[box].append(genus)
                        else:
                            taxa_in_box[box] = [genus]

        for box in taxa_in_box.keys():
            nrboxes = len(taxa_in_box[box])
            if nrboxes > maxshared:
                nrboxes = maxshared
            paired = True if nrboxes % 2 == 0 else False
            for i in range(nrboxes):
                newrecord = pathway.create_box_heatmap(pathway.pathway.orthologs[box], nrboxes, 
                                                    i * 2 - (nrboxes - 1) if paired     # if nrboxes = 8, i * 2 - (nrboxes - 1) = -7,-5,-3,-1,1,3,5,7
                                                    else i - int(nrboxes / 2),          # if nrboxes = 9, i - int(nrboxes / 2) = -4,-3,-2,-1,0,1,2,3,4
                                                    paired = paired)
                newrecord.bgcolor = dic_colors[taxa_in_box[box][i]]
                pathway.pathway.orthologs[box].graphics.append(newrecord)
            pathway.create_tile_box(pathway.pathway.orthologs[box])
        
        df = data[samples + ['KO (KEGG Pathway)']]
        df = df[df.any(axis=1)]
        grey_boxes = list()
        for ortholog in df['KO (KEGG Pathway)']:
            if ortholog in pathway.ko_boxes.keys():
                for box in pathway.ko_boxes[ortholog]:
                    if box not in taxa_in_box.keys() and box not in grey_boxes:
                        grey_boxes.append(box)
        pathway.grey_boxes(grey_boxes)
        
        name_pdf = '{}_{}.pdf'.format(output_basename, metabolic_map)              #maps[metabolic_map].replace('/','|').replace(' ','_'), 
        pathway.pathway_pdf(name_pdf)
        '''
        grey_boxes = [1]
        name_pdf = '{}_{}.pdf'.format(output_basename, metabolic_map)
        if len(grey_boxes) > 0:
            colors.append("#7c7272")
            genera.append("Present in samples")
        
        self.create_potential_legend(colors, genera, 
                                     name_pdf.replace('.pdf','_legend.png'))
            
        self.add_legend(name_pdf, name_pdf.replace('.pdf','_legend.png'), 
                        name_pdf.replace(metabolic_map, self.maps[metabolic_map].replace('/','_')))
        
        for file in [name_pdf, name_pdf.replace('.pdf','_legend.png'), 
                     name_pdf.replace('.pdf', '.png')]:
            os.remove(file)
            
    def differential_expression_sample(self, data, samples, output_basename = None,
                                       log = False, metabolic_map = None):
        '''
        Represents in small heatmaps the expression levels of each sample on the
        dataset present in the given pathway map. The values can be transford to
        a log10 scale
        :param data: pandas.DataFrame with data already processed by KEGGPathway
        :param samples: list - column names of the dataset corresponding to
        expression values
        :param output_folder: string - name of folder to store pdfs
        :param log: bol - convert the expression values to logarithmic scale?
        '''
        pathway = KeggMap(data, metabolic_map)
        df = data.groupby('KO (KEGG Pathway)')[samples + ['KO (KEGG Pathway)']].sum()
        df = df[df.any(axis=1)]
        df['Boxes'] = [pathway.ko_boxes[ko] if ko in pathway.ko_boxes.keys()
                        else np.nan for ko in df.index]
        df = df[df['Boxes'].notnull()]
        df = mtools.expand_by_list_column(df, column = 'Boxes')                 #                sample1
                                                                                # 6              151800.0
        df = df.groupby('Boxes')[samples].sum()

        pathway.pathway_boxes_diferential(df, log)
        
        name_pdf = '{}_{}{}.pdf'.format(output_basename, metabolic_map,               #maps[metabolic_map].replace('/','|').replace(' ','_'), 
                    '_log' if log else '')
        pathway.pathway_pdf(name_pdf)
            
        self.differential_colorbar(df, name_pdf.replace(".pdf",'_legend.png'))
        
        self.add_legend(name_pdf, name_pdf.replace('.pdf','_legend.png'), 
                        name_pdf.replace(metabolic_map + '.pdf', 
                            self.maps[metabolic_map].replace('/','_') + '.png'))
        
        for file in [name_pdf, name_pdf.replace('.pdf','_legend.png'), 
                     name_pdf.replace('.pdf', '.png')]:
            os.remove(file)


    def differential_colorbar(self, dataframe, filename):
        FIGSIZE = (2,3)
        mpb = plt.pcolormesh(dataframe,cmap='coolwarm')
        fig,ax = plt.subplots(figsize=FIGSIZE)
        plt.colorbar(mpb,ax=ax)
        ax.remove()
        plt.savefig(filename,bbox_inches='tight')
        
    def solve_ec_numbers(self, data):
        '''
        The KEGG API manipulation requires alterations of the data that end up
        with many repeated lines for where there are many EC numbers. 
        Also, UniProt's mapping retrieves a list of EC numbers that may (and will,
        most likely) differ from the list from the KEGG API. In this method, all
        all this information is merged in one column, "EC numbers", and the 
        information of both methods is compared and merged.
        :param data: pd.DataFrame - with 'EC number' and 'EC number (KEGG Pathway)'
        columns
        :returns data with new column 'EC numbers' with integrated information
        from UniProt mapping and KEGG API
        '''
        notnan = data[data['EC number (KEGG Pathway)'].notnull()]
        notnan.drop_duplicates(inplace = True)
        notnan = notnan.groupby('Entry')['EC number (KEGG Pathway)'].apply(
                lambda ecs:'; '.join([ec.split(':')[1] for ec in ecs]))
        del data['EC number (KEGG Pathway)']
        data.drop_duplicates(inplace = True)
        data = pd.merge(data, notnan, on = 'Entry', how = 'outer')
        joined_ecs = list()
        for i in range(len(data)):
            if data.iloc[i]['EC number'] == data.iloc[i]['EC number (KEGG Pathway)']:   # EC number is the same for both methods
                joined_ecs.append(data.iloc[i]['EC number'])
            else:
                if type(data.iloc[i]['EC number']) == float:                            # No EC number from UniProt's mapping...
                    if type(data.iloc[i]['EC number (KEGG Pathway)']) == float:         # ... and none from KEGG API either
                        joined_ecs.append(np.nan)
                    else:                                                               # ... and KEGG API's got it
                        joined_ecs.append(data.iloc[i]['EC number (KEGG Pathway)'])
                else:                                                                   # There is EC number from UniProt ...
                    if type(data.iloc[i]['EC number (KEGG Pathway)']) == float:         # ... and there is none from KEGG API
                        joined_ecs.append(data.iloc[i]['EC number'])
                    else:                                                               # ... and there a different one from KEGG API. uh-oh
                        uniprot_ecs = data.iloc[i]['EC number'].split('; ')
                        kegg_ecs = data.iloc[i]['EC number (KEGG Pathway)'].split('; ')
                        result = str()
                        if len([ec for ec in uniprot_ecs if ec in kegg_ecs]) > 0:       # The EC numbers present in both methods are added as is
                            result += '; '.join(ec for ec in uniprot_ecs if ec in kegg_ecs) + '; '
                        if len(set(uniprot_ecs) - set(kegg_ecs)) > 0:                   # Checks if there are any UniProt EC numbers not in KEGG API and adds those
                            result += ' (uniprot_map); '.join([ec for ec in uniprot_ecs if ec not in kegg_ecs]) + ' (uniprot_map); '
                        if len(set(kegg_ecs) - set(uniprot_ecs)) > 0:                   # Checks if there are any KEGG API EC numbers not in UniProt and adds those
                            result += ' (kegg_api); '.join([ec for ec in kegg_ecs if ec not in uniprot_ecs]) + ' (kegg_api); '
                        joined_ecs.append(result[:-2])
        data['EC numbers'] = joined_ecs
        return data
    
    def get_organisms(self, file):
        '''
        Data obtained from www.genome.jp/kegg/catalog/org_list.html was organized
        to retrieve the identifiers of the organisms available at KEGG
        :param file: string - filename containing the information
        '''
        return pd.read_csv(file, sep = '\t', index_col = 0, header = None)
        
    
    def solve_kegg_ids(self, data):
        '''
        UniProt ID mapping returns columns of 'Cross-reference (KEGG)' with messed
        up values like mfc:BRM9_0145;mfi:DSM1535_1468;. This makes no sense when
        there are taxonomic columns that clearly define the taxon of the protein
        and so this method cleans the values that make no sense
        Seems for now this must be put in wait, as UniProt doesn't report on the
        strain of the taxons
        :param data: pd.DataFrame - fresh from MOSCA's analysis, with messed up
        'Cross-reference (KEGG)' column
        '''
        for i in range(len(data)):
            if type(data.iloc[i]['Cross-reference (KEGG)']) != float:
                kegg_ids = data.iloc[i]['Cross-reference (KEGG)'].split(':')
                if len(kegg_ids) > 1:
                    correct_taxon = data.iloc[i]['Taxonomic lineage (SPECIES)']
                    print(correct_taxon)
        return data

    def run(self, input_file, output_directory, mg_samples = None, mt_samples = None, 
            metabolic_maps = None):
        '''
        Represents in small heatmaps the expression levels of each sample on the
        dataset present in the given pathway map. The values can be transford to
        a log10 scale
        :param input_file: CSV data file outputed by MOSCA analysis
        :param output_directory: str - name of folder for outputs
        :param mg_samples: list [str] - name of columns containing MG quantification
        for genomic potential representations
        :param mt_samples: list [str] - name of columns containing MT quantification
        for differential expression representations
        '''
        
        if metabolic_maps is None:
            metabolic_maps = self.default_maps
        
        data = pd.read_csv(input_file, sep = '\t', low_memory = False)
        '''
        kegg_ids = data[data['Cross-reference (KEGG)'].notnull()]['Cross-reference (KEGG)']
        kegg_ids = [ide.split(';')[0] for ide in kegg_ids]                      # should be fixed by 'solve_kegg_ids', but not yet possible because of UniProt
        kos = self.keggid2ko(kegg_ids)
        data = pd.merge(data, kos, on = 'Cross-reference (KEGG)', how = 'outer')

        kos = data[data['KO (KEGG Pathway)'].notnull()]['KO (KEGG Pathway)'].tolist()
        ecs = self.ko2ec(kos)
        data = pd.merge(data, ecs, on = 'KO (KEGG Pathway)', how = 'outer')
        
        data = self.solve_ec_numbers(data)
        data.to_csv(input_file, sep = '\t', index = False)
        '''
        print('Creating KEGG Pathway representations for ' + str(len(metabolic_maps)) + 
              ' metabolic pathways.')
        genomic = list(); differential = list()
        mt_samples = None
        for metabolic_map in metabolic_maps:
            print('Creating representation for pathway: ' + self.maps[metabolic_map])
            if mg_samples:
                try:
                    self.genomic_potential_taxa(data, mg_samples, 
                            output_basename = output_directory + '/potential',
                            metabolic_map = metabolic_map)
                except:
                    genomic.append(metabolic_map)
            if mt_samples:
                try:
                    self.differential_expression_sample(data, mt_samples, 
                        output_basename = output_directory + '/differential',
                        metabolic_map = metabolic_map)
                except:
                    differential.append(metabolic_map)
            plt.cla()
        genomic_failed = output_directory + '/genomic_failed_maps.txt'
        differential_failed = output_directory + '/differential_failed_maps.txt'
        print(('Failed {} maps for genomic potential representation. You can consult' + 
              'which ones at {}').format(len(genomic), genomic_failed))
        open(genomic_failed, 'w').write('\n'.join(genomic))
        print(('Failed {} maps for differential expression representation. You can consult' + 
              'which ones at {}').format(len(differential), differential_failed))
        open(differential_failed, 'w').write('\n'.join(genomic))
        
class KeggMap():
    '''
    This class retrieves and manipulates KEGG metabolic maps from KEGG Pathway
    '''

    def __init__(self, data, pathway_ID, **kwargs):
        '''
        Initialize object
        :param data: pd.DataFrame - data from MOSCA analysis
        :param pathway_ID: (str) - KEGG Pathway ID
        '''
        self.__dict__ = kwargs
        self.pathway_ID = pathway_ID[-5:] if len(pathway_ID) > 5 else pathway_ID
        self.set_pathway(data, pathway_ID)

    ############################################################################
    ####                              Helper                                ####
    ############################################################################
    
    def set_bgcolor(self, pathway_element, color):
        '''
        Sets graphic element background color
        :param pathway_element: kegg pathway xml element object
        :param color: color to be used in rgb
        '''
        pathway_element.graphics[0].bgcolor = color

    def set_fgcolor(self, pathway_element, color):
        '''
        Sets graphic element contour color
        :param pathway_element: kegg pathway xml element object
        :param color:  color to be used in rgb
        '''
        pathway_element.graphics[0].fgcolor = color
        
    def set_colors(self, colors = [], ncolor = 1):
        '''
        Creates list of hex colors to be used, using matplotlib or using custom colors
        :param colors: list of hex colors
        :param ncolor: int indicating the amount of hex colors that should be created
        :return: returns list with hex color codes
        '''
        if len(colors) == 0:
            # if no colors are given creates a list of hex colors with ncolor
            # ammount from matplotlib discrete colormaps. If ncolor > 12
            # a continuous colormap is used instead.
            if ncolor <= 8:
                pastel2 = cm.get_cmap('Pastel2', 8)
                return [to_hex(pastel2(i)) for i in range(ncolor)]
            elif ncolor <= 12:
                set3 = cm.get_cmap("Set3", 12)
                return [to_hex(set3(i)) for i in range(ncolor)]
            else:
                rainbow = cm.get_cmap("rainbow",ncolor)
                return [to_hex(rainbow(i)) for i in range(ncolor)]
        else:
            # validates hex values and returns the original list
            isvalidhex = True
            i = 0
            while isvalidhex:
                match = re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', str(colors[i]))
                if not match:
                    isvalidhex = False
                i += 1
            if isvalidhex:
                return colors
            else:
                raise Exception("Colors aren't valid hex codes")

    def conv_value_rgb(self, value, colormap, norm):
        '''
        Normalizes values in a vector and transforms it into corresponding
        hex color using given colormap
        :param value: numpy vector from dataframe apply function with expression
        values
        :param colormap: matplotlib colormap object
        :param norm: matplotlib Normalize or LogNorm object
        :return: returns hex colors in vector format
        '''
        return value.apply(norm).apply(colormap)

    def conv_rgb_hex(self, rgb):
        '''
        converts rgb into hex color code in vector
        :param rgb: rgb value in vector resulting from pandas dataframe apply
        :return: vector with converted hex color code
        '''
        return rgb.apply(to_hex)

    def load_kegg_map(self, pathway_ID, organism_ID = ""):
        '''
        Downloads pathway kgml from KEGG and reads it
        :param pathway_ID: (str) - suffix Pathway ID, int part
        :param organism_ID: (str) - preffix Pathway ID, str part
        :return: (object) - pathway kgml parsed
        '''
        if not organism_ID:
            pathway = kegg_get("ko" + pathway_ID, "kgml")
            return KGML_parser.read(pathway)
        else:
            try:
                pathway = kegg_get(organism_ID + pathway_ID, "kgml")
                return KGML_parser.read(pathway)
            except:
                print("Invalid KEGG map ID")
                
    def reset_pathway(self):
        '''
        Resets pathway state
        '''
        self.set_pathway(self.pathway_ID)

    def organismo_genes(self, kegg_ID):
        '''
        returns list of organism ids from a given list of kegg gene id
        :param kegg_ID:list with kegg gene id
        :return:list with kegg org ids
        '''
        organisms = []
        for gene_ID in kegg_ID:
            organisms.append(gene_ID.split(":")[0])
        return list(set(organisms))

    def ortholog_dic(self):
        '''
        Dictionary associating all ortholog names in pathway element
        with its first ortholog
        :return: dic associaiting orthologs present in pathway elements
        '''
        orthologs_dic = {}
        for ortholog_rec in self.pathway.orthologs:
            ortholog_name = ortholog_rec.name.split(" ")[0]
            for ortholog in ortholog_rec.name.split(" "):
                if ortholog not in orthologs_dic.keys():
                    orthologs_dic[ortholog] = ortholog_name
        return orthologs_dic
    
    ############################################################################
    ####                            Sets                                    ####
    ############################################################################

    def set_pathway(self, data, pathway_ID):
        '''
        Set pathway with Kegg Pathway ID
        :param pathway_ID: (str) Kegg Pathway ID
        '''
        self.pathway = self.load_kegg_map(self.pathway_ID)                       # get the KGML
        ko = []
        self.ko_boxes = {}
        for i in range(len(self.pathway.orthologs)):
            self.set_bgcolor(self.pathway.orthologs[i], "#ffffff")              # set all boxes to white
            self.set_fgcolor(self.pathway.orthologs[i], "#ffffff")              # ditto
            orthologs_in_box = [ide[3:] for ide in self.pathway.orthologs[i].name.split()]  # 'ko:K16157 ko:K16158 ko:K16159' -> ['K16157', 'K16158', 'K16159']
            for ortholog in orthologs_in_box:
                if ortholog not in self.ko_boxes.keys():
                    self.ko_boxes[ortholog] = list()
                self.ko_boxes[ortholog].append(i)                               # {'K16157':[0,13,432], 'K16158':[4,13,545]}
            ko.append(self.pathway.orthologs[i].graphics[0].name.rstrip("."))   # 'K16157...' -> 'K16157'
        
        # Set text in boxes to EC numbers
        data = data[data['EC number (KEGG Pathway)'].notnull()][[
                'KO (KEGG Pathway)', 'EC number (KEGG Pathway)']]
        ko_to_ec = {data.iloc[i]['KO (KEGG Pathway)']:data.iloc[i]['EC number (KEGG Pathway)']
                    for i in range(len(data))}     # {'K16157':'ec:1.14.13.25'}
        for ortholog_rec in self.pathway.orthologs:
            ko = ortholog_rec.graphics[0].name.strip(".")
            if ko in ko_to_ec.keys():
                ortholog_rec.graphics[0].name = ko_to_ec[ko]

    ############################################################################
    ####                    Graphical Manipulation                          ####
    ############################################################################

    def create_tile_box(self, record):
        '''
        Create box graphical element in pathway to draw the box countour and
        give the correct name
        :param record: graphical element to be created
        '''
        newrecord = KGML_pathway.Graphics(record)
        newrecord.name = record.graphics[0].name
        newrecord.type = "rectangle"
        newrecord.width = record.graphics[0].width
        newrecord.height = record.graphics[0].height
        newrecord.y = record.graphics[0].y
        newrecord.x = record.graphics[0].x
        newrecord.bgcolor = "#FFFFFF00"
        newrecord.fgcolor = "#000000"
        record.graphics.append(newrecord)
        record.graphics[0].bgcolor = "#FFFFFF00"
        record.graphics[0].fgcolor = "#FFFFFF00"
        record.graphics[0].name = ""

    def create_box_heatmap(self, rec_old, nrboxes, i, paired = True):
        '''
        Helper function for creating heatmap, draws one expression value in its
        correct position on the bigger parent box
        :param rec_old: graphical element object to be used as reference
        :param nrboxes: int nr of boxes to be drawed
        :param i: int internal number of movements of the box given by the for loop
        :return: graphical element object
        '''
        movement_steps = rec_old.graphics[0].width / (nrboxes * (2 if paired else 1))
        newrecord = KGML_pathway.Graphics(rec_old)
        newrecord.name = ""
        newrecord.type = "rectangle"
        adjustment_factor = 1.3 if nrboxes > 2 else 1.1 if nrboxes > 1 else 1   # sub-boxes width, adjusted by a factor that experimentally fixed well in the representations
        newrecord.width = movement_steps * adjustment_factor * (2 if paired else 1)
        newrecord.height = rec_old.graphics[0].height
        newrecord.y = rec_old.graphics[0].y
        newrecord.x = (i * movement_steps) + rec_old.graphics[0].x
        newrecord.fgcolor = "#FFFFFF00"
        return newrecord

    ############################################################################
    ####                          Operations                                ####
    ############################################################################

    def pathway_pdf(self, filename, imagemap = True, orthologs = True, 
                    compounds = True, maps = True, reactions = True):
        '''
        Prints current pathway to PDF file
        :param filename: (str) - PDF filename
        :param imagemap: (bol) - Print imagemap
        :param orthologs: (bol) - Print orthologs
        :param compounds: (bol) - Print compounds
        :param maps: (bol) - Print maps
        :param reactions: (bol) - Print reactions ???
        :return: creates PDF file with current pathway
        '''
        #TODO Verificar o parametro reactions
        KGMLCanvas(self.pathway, 
                   import_imagemap = imagemap,
                   label_orthologs = orthologs, 
                   label_compounds = compounds,
                   label_maps = maps, 
                   label_reaction_entries = reactions).draw(filename)

    def pathway_box_list(self, taxa_in_box, taxa, maxshared = 10, colors = None):
        '''
        Represents items in the pathway map
        :param taxa_in_box: dict - {box : list of taxa in box}
        :param taxa: list - of taxa to be represented and given a specific color
        :param maxshared: int - maximum number of taxa sharing one box
        :param color: list of costum colors to be used to color the elements
        '''
        if colors is None:
            colors = KEGGPathway.taxa_colors(ncolor = len(taxa))
        dic_colors = {taxa[i] : colors[i] for i in range(len(taxa))}

        for box in taxa_in_box.keys():
            nrboxes = len(taxa_in_box[box])
            if nrboxes > maxshared:
                nrboxes = maxshared
                
            paired = True if nrboxes % 2 == 0 else False
            for i in range(nrboxes):
                newrecord = self.create_box_heatmap(self.pathway.orthologs[box], nrboxes, 
                                                    i * 2 - (nrboxes - 1) if paired     # if nrboxes = 8, i * 2 - (nrboxes - 1) = -7,-5,-3,-1,1,3,5,7
                                                    else i - int(nrboxes / 2),          # if nrboxes = 9, i - int(nrboxes / 2) = -4,-3,-2,-1,0,1,2,3,4
                                                    paired = paired)
                newrecord.bgcolor = dic_colors[taxa_in_box[box][i]]
                self.pathway.orthologs[box].graphics.append(newrecord)
            self.create_tile_box(self.pathway.orthologs[box])

    def pathway_boxes_diferential(self, dataframe, log = False, colormap = "coolwarm"):
        '''
        Represents expression values present in a dataframe in the
        pathway map
        :param dataframe: pandas DataFrame with each column representing a sample
        and index corresponding to int list index of the ortholog element in the
        pathway
        :param log: bol providing the option for a log normalization of data
        :param colormap: str representing a costum matplotlib colormap to be used
        '''

        if log:
            norm = cm.colors.LogNorm(vmin=dataframe.min().min(), vmax=dataframe.max().max())
        else:
            norm = cm.colors.Normalize(vmin=dataframe.min().min(), vmax=dataframe.max().max())

        colormap = cm.get_cmap(colormap)
        dataframe = dataframe.apply(self.conv_value_rgb, args=(colormap, norm))
        dataframe = dataframe.apply(self.conv_rgb_hex)
        
        dataframe = dataframe[dataframe.columns.tolist()]
        
        nrboxes = len(dataframe.columns.tolist())                               # number of samples
        
        for box in dataframe.index.tolist():
            colors = dataframe.loc[box].tolist()
            paired = True if nrboxes % 2 == 0 else False
            for i in range(nrboxes):
                newrecord = self.create_box_heatmap(self.pathway.orthologs[box], nrboxes, 
                                                    i * 2 - (nrboxes - 1) if paired     # if nrboxes = 8, i * 2 - (nrboxes - 1) = -7,-5,-3,-1,1,3,5,7
                                                    else i - int(nrboxes / 2),          # if nrboxes = 9, i - int(nrboxes / 2) = -4,-3,-2,-1,0,1,2,3,4
                                                    paired = paired)
                newrecord.bgcolor = colors[i]
                self.pathway.orthologs[box].graphics.append(newrecord)
            self.create_tile_box(self.pathway.orthologs[box])

    def grey_boxes(self, box_list):
        for i in box_list:
            self.set_bgcolor(self.pathway.orthologs[i], "#7c7272")
            self.set_fgcolor(self.pathway.orthologs[i], "#7c7272")

if __name__ == '__main__':

    kp = KEGGPathway()
    '''
    kp.run(input_file = 'MGMP/all_info_normalized.tsv',
           output_directory = 'MGMP/KEGGPathway',
           mg_samples = ['EST6_S1_L001', 'OL6_S3_L001', 'OLDES6_S4_L001', 'PAL6_S2_L001'])
    '''
    kp.run(input_file = 'MOSCAfinal/all_info_normalized.tsv',
           output_directory = 'MOSCAfinal/KEGGPathway',
           mg_samples = ['4478-DNA-S1613-MiSeqKapa', '4478-DNA-S1616-MiSeqKapa', 
                         '4478-DNA-S1618-MiSeqKapa'],
           mt_samples = ['4478-R1-1-MiSeqKapa_normalized','4478-R2-1-MiSeqKapa_normalized', 
                         '4478-R3-1-MiSeqKapa_normalized','4478-R4-1-MiSeqKapa_normalized'])