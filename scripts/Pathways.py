#!/usr/bin/env python

from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.KEGG.KGML import KGML_pathway
from matplotlib import cm
from matplotlib.colors import to_hex
import matplotlib.pyplot as plt
import re
import pandas as pd
import numpy as np
import progressbar

__author__ = "Tiago Oliveira"
__credits__ = ["Joao Sequeira"]
__version__ = "1.0"
__maintainer__ = "Tiago Oliveira"
__email__ = "tiagomanuel28@gmail.com"
__status__ = "Production"
################################################################################
class Pathway():

    def __init__(self, pathway_ID):
        '''
        Initialize object
        :param pathway_ID: (str) - KEGG Pathway ID
        '''
        self.set_pathway(pathway_ID)
        self.KEGG_ID = []
        self.ortholog_ID = []

    def __str__(self):
        '''
        Print pathway
        :return: (print) - global pathway info
        '''
        print(self.pathway)

    ############################################################################
    ####                            Hiden                                   ####
    ############################################################################
    def _pathway_ID_handle(self, pathway_ID):
        '''
        Parses the pathway ID, extracting the numeric parth
        :param pathway_ID: (str) - full pathway ID
        :return: (str) pathway ID, int part
        '''
        if len(pathway_ID) >5:
            pathway_ID = pathway_ID[-5:]
        return str(pathway_ID)

    def _set_bgcolor(self, pathway_element, color):
        '''
        Sets graphic element background color
        :param pathway_element: kegg pathway xml element object
        :param color: color to be used in rgb
        '''
        pathway_element.graphics[0].bgcolor = color

    def _set_fgcolor(self, pathway_element, color):
        '''
        Sets graphic element contour color
        :param pathway_element: kegg pathway xml element object
        :param color:  color to be used in rgb
        '''
        pathway_element.graphics[0].bgcolor = color

    def _set_colors(self, colors = [], ncolor = 1):
        '''
        Creates list of hex colors to be used,
        using matplotlib or using costum colors
        :param colors: list of hex colors
        :param ncolor: int indicating the ammount of hex colors should be created
        :return: returns list with hex color codes
        '''
        ret_color = []
        if len(colors) == 0:
            # if no colors are given creates a list of hex colors with ncolor
            # ammount from matplotlib discrete colormaps. If ncolor > 12
            # a continuous colormap is used instead.
            if ncolor <= 8:
                pastel2 = cm.get_cmap('Pastel2', 8)
                for i in range(ncolor):
                    ret_color.append(to_hex(pastel2(i)))
                return ret_color
            elif ncolor <= 12:
                set3 = cm.get_cmap("Set3", 12)
                for i in range(ncolor):
                    ret_color.append(to_hex(set3(i)))
                return ret_color
            else:
                rainbow = cm.get_cmap("rainbow",ncolor)
                for i in range(ncolor):
                    ret_color.append(to_hex(rainbow(i)))
                return ret_color
        else:
            # validates hex values and returns the original list
            isvalidhex = True
            for hexvalue in colors:
                match = re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', str(hexvalue))
                if not match:
                    isvalidhex = False
            if isvalidhex:
                return colors
            else:
                raise Exception("Colors aren't valid hex codes")

    def _conv_value_rgb(self, value, colormap, norm):
        '''
        Normalizes values in a vector and transforms it into corresponding
        hex color using given colormap
        :param value: numpy vector from dataframe apply function with expression
        values
        :param colormap: matplotlib colormap object
        :param norm: matplotlib Normalize or LogNorm object
        :return: returns hex colors in vector format
        '''
        value = value.apply(norm)
        return value.apply(colormap)

    def _conv_rgb_hex(self, rgb):
        '''
        converts rgb into hex color code in vector
        :param rgb: rgb value in vector resulting from pandas dataframe apply
        :return: vector with converted hex color code
        '''
        return rgb.apply(to_hex)

    ############################################################################
    ####                              Helper                                ####
    ############################################################################

    def load_pathway(self, pathway_ID, organism_ID = ""):
        '''
        Downloads pathway kgml from KEGG and readis it
        :param pathway_ID: (str) - suffix Pathway ID, int part
        :param organism_ID: (str) - preffix Pathway ID, str part
        :return: (object) - pathway kgml parsed
        '''
        if not organism_ID:
            print(pathway_ID)
            pathway = kegg_get(str("ko"+pathway_ID), "kgml")
            return KGML_parser.read(pathway)
        else:
            try:
                pathway = kegg_get(str(organism_ID + pathway_ID), "kgml")
                return KGML_parser.read(pathway)
            except:
                print("Invalid IDs")

    def reset_pathway(self):
        '''
        Resets pathway state
        '''
        self.set_pathway(self.pathway_ID)
        self.KEGG_ID = []
        self.ortholog_ID = []

    def covnert_KEGGID_koID(self, KEGG_ID):
        '''
        Converts KEGG_ID genes to Ortholog KO ID from KEGG
        :param KEGG_ID: (list) - KEGG ID genes
        :return: (tupple) - (list,list) - KEGG ID genes converted and ko IDs
        '''
        #TODO Adicionar excecoes para gerir KEGG IDs incorretos

        coverted_KEGG_ID = []
        orthologs_ID = []
        links = kegg_link("ko", KEGG_ID).read().split("\n")
        for link in links:
            if len(link) > 1:
                coverted_KEGG_ID.append(link.split("\t")[0])
                orthologs_ID.append(link.split("\t")[1])
        return (coverted_KEGG_ID, orthologs_ID)

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

    def convert_koID_EC(self,ko_ID):
        '''
        Converts kegg ortholog id to EC numers
        :param ko_ID: list of kegg ortholog ids
        :return: dic associating ortholog kegg id with list
        of assotiated EC numbers
        '''
        ko_dic = {}
        conversion = kegg_link("enzyme", ko_ID).read().split("\n")
        for result in conversion:
            if len(result) > 1:
                ko_dic[result.split("\t")[0].strip("ko:")] = result.split("\t")[1].upper()
        return ko_dic

    ############################################################################
    ####                            Sets                                    ####
    ############################################################################

    def set_pathway(self,pathway_ID):
        '''
        Set pathway with Kegg Pathway ID
        :param pathway_ID: (str) Kegg Pathway ID
        '''
        #TODO depricate self.boxes_ko
        self.pathway_ID = self._pathway_ID_handle(pathway_ID)
        self.pathway = self.load_pathway(self.pathway_ID)
        ko = []

        self.boxes_ko = {}
        self.ko_boxes = {}
        for i in range(0, len(self.pathway.orthologs)):
            self._set_bgcolor(self.pathway.orthologs[i], "#9e9e9e")
            self._set_fgcolor(self.pathway.orthologs[i], "#9e9e9e")
            orthologs_in_box = self.pathway.orthologs[i].name.split(" ")
            for n in range(0, len(orthologs_in_box)):
                orthologs_in_box[n] = orthologs_in_box[n][3:]
                if orthologs_in_box[n] in self.ko_boxes.keys():
                    self.ko_boxes[orthologs_in_box[n]].append(i)
                else:
                    self.ko_boxes[orthologs_in_box[n]] = [i]
            self.boxes_ko[i] = orthologs_in_box
            ko.append(self.pathway.orthologs[i].graphics[0].name.strip("."))

        ko_to_ec = self.convert_koID_EC(ko)
        for ortholog_rec in self.pathway.orthologs:
            ko = ortholog_rec.graphics[0].name.strip(".")
            if ko in ko_to_ec.keys():
                ortholog_rec.graphics[0].name = ko_to_ec[ko]

    def set_kegg_ID(self, kegg_ID):
        '''
        Sets keeg ids
        :param kegg_ID: list of kegg ids
        '''
        self.KEGG_ID = kegg_ID

    def set_ortholog_ID(self, ortholog_ID):
        '''
        Set ortholog ids
        :param ortholog_ID: list of ortholog ids
        '''
        self.ortholog_ID = ortholog_ID

    ############################################################################
    ####                            Gets                                    ####
    ############################################################################

    def get_pathway_ID(self):
        '''
        Get pathway ID
        :return: pathway ID
        '''
        return self.pathway_ID

    def get_pathway(self):
        '''
        returns pathway
        :return: pathway object
        '''
        return self.pathway

    def get_orthologs_ID(self):
        '''
        returns ortholog id list
        :return: list ortholog id
        '''
        return self.ortholog_ID

    def get_ko_boxes(self):
        '''
        returns dictionary assotiating ortholog kegg id to list of
        ortholog element index
        :return: dic with ortholog kegg id as keys and values as list of
        ortholog element index number
        '''
        return self.ko_boxes

    def get_boxes_ko(self):
        '''
        returns dictionary assotiang ortholog element index number to list of
        ortholog kegg ids
        :return: dic with ortholog element id int as keys and values as list of
        ortholog kegg ids
        '''
        return self.boxes_ko

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

    def create_pair_boxes(self, rec_old, nrboxes, i):
        '''
        Creates small box when the number of boxes to draw is an even number
        :param rec_old: graphical element object to be used as reference
        :param nrboxes: int nr of boxes to be drawed
        :param i: int internal number of movements of the box given by the for loop
        :return: graphical element object
        '''
        movement_steps = rec_old.graphics[0].width / (4 * (nrboxes / 2))
        width = rec_old.graphics[0].width / nrboxes
        rec_new = KGML_pathway.Graphics(rec_old)
        rec_new.name = ""
        rec_new.type = "rectangle"
        rec_new.width = width
        rec_new.height = rec_old.graphics[0].height
        rec_new.y = rec_old.graphics[0].y
        rec_new.x = (i * movement_steps) + rec_old.graphics[0].x
        rec_new.fgcolor = "#FFFFFF00"
        return rec_new

    def create_odd_boxes(self, rec_old, nrboxes, i):
        '''
        Creates small box when the number of boxes to draw is an odd number
        :param rec_old:graphical element object to be used as reference
        :param nrboxes:int nr of boxes to be drawed
        :param i:int internal number of movements of the box given by the for loop
        :return: graphical element object
        '''
        width = rec_old.graphics[0].width / nrboxes
        movement_steps = width
        n = 0
        newrecord = KGML_pathway.Graphics(rec_old)
        newrecord.name = ""
        newrecord.type = "rectangle"
        newrecord.width = width
        newrecord.height = rec_old.graphics[0].height
        newrecord.y = rec_old.graphics[0].y
        newrecord.x = (i * movement_steps) + rec_old.graphics[0].x
        newrecord.fgcolor = "#FFFFFF00"
        return newrecord

    ############################################################################
    ####                          Operations                                ####
    ############################################################################

    def pathway_pdf(self, filename = "", imagemap = True, orthologs = True, compounds = True, maps = True, reactions = True):
        '''
        Prints current pathway to PDF file
        :param filename: (str) - PDF filename
        :param imagemap: (bol) - Print imagemap
        :param orthologs: (bol) - Print orthologs
        :param compounds: (bol) - Print compounds
        :param maps: (bol) - Print maps
        :param reactions: (bol) - Print reactions ???
        :return: PDF file with current pathway
        '''
        #TODO Verificar o parametro reactions
        pathway_pdf = KGMLCanvas(self.pathway)
        pathway_pdf.import_imagemap=imagemap
        pathway_pdf.label_orthologs = orthologs
        pathway_pdf.label_compounds = compounds
        pathway_pdf.label_maps = maps
        pathway_pdf.label_reaction_entries = reactions

        if filename == "":
            pathway_pdf.draw(str(self.pathway_ID+".pdf"))
        else:
            pathway_pdf.draw(str(filename+".pdf"))

    def pathway_genes(self, KEGG_ID, color=[]):
        '''
        Given a list of kegg gene ids colors the pathway map with a specific color
        :param KEGG_ID: list of kegg gene ids following this sintax:
         speciescod:geneID
        :param color: hex color code
        '''
        self.KEGG_ID, ortholog_ID = self.convert_KEGGID_koID(KEGG_ID)
        color1 = self._set_colors(color)[0]
        for ortholog_rec in self.pathway.orthologs:
            for ortholog in ortholog_rec.name.split(" "):
                if ortholog in ortholog_ID:
                    self.ortholog_ID.append(ortholog)
                    self._set_bgcolor(ortholog_rec, color1)
                    self._set_fgcolor(ortholog_rec, color1)

    def pathway_organismo(self, organism_ID, color = ""):
        '''
        Given a list of kegg organism ids, colors the pathway map with
        a specific color
        :param organism_ID: lsit of species kegg id
        :param color: hex color code
        '''
        #TODO adaptar para recolher varios organismos

        color1 = self._set_colors(color)[0]

        pathway_organism = self.load_pathway(self.pathway_ID, organism_ID)
        organism_kegg_ID = []
        for gene_rec in pathway_organism.genes:
            for gene in gene_rec.name.split(" "):
                organism_kegg_ID.append(gene)
        organism_ortholog_ID = self.convert_KEGGID_koID(organism_kegg_ID)[1]

        for ortholog_rec in self.pathway.orthologs:
            for ortholog in ortholog_rec.name.split(" "):
                if ortholog in organism_ortholog_ID:
                    self.ortholog_ID.append(ortholog)
                    self._set_bgcolor(ortholog_rec, color1)
                    self._set_fgcolor(ortholog_rec, color1)

    def pathway_genes_organismo(self, KEGG_ID, maxshared = 5, color = []):
        '''
        Given a list of gene_kegg ids colors the the pathway map with
        a specific color for organism
        :param KEGG_ID: list of kegg gene id
        :param maxshared: naximum number of boxes to be drawed
        :param color: list of hex color codes to be used to color the pathway map
        '''
        coverted_KEGG_ID, orthologs_ID = self.convert_KEGGID_koID(KEGG_ID)
        organisms = self.organismo_genes(coverted_KEGG_ID)

        # Set colors
        organism_colors = {}
        if len(color) < len(organisms):
            colors = self._set_colors([], len(organisms))
        else:
            colors = self._set_colors(color)
        i = 0
        for organism in organisms:
            if organism not in organism_colors.keys():
                organism_colors[organism] = colors[i]
                i += 1

        # Dictionario Ortholog - [organismo]
        ortholog_dic = {}
        for i in range(len(orthologs_ID)):
            if orthologs_ID[i-1] in ortholog_dic.keys():
                ortholog_dic[orthologs_ID[i - 1]].append(organisms[i - 1])
            else:
                ortholog_dic[orthologs_ID[i - 1]] = [organisms[i - 1]]

        for ortholog_rec in self.pathway.orthologs:
            # organismos por ortholog_rec
            organisms = []
            for ortholog in ortholog_rec.name.split(" "):
                if ortholog in ortholog_dic.keys():
                    organisms += ortholog_dic[ortholog]
            # Se tiver organismos para o ortholog_rec
            if len(organisms) > 0:
                if len(organisms) <= maxshared:
                    organism_size = len(organisms)
                else:
                    organism_size = maxshared
                organisms.sort()
                if organism_size%2==0:
                    n = 0
                    for i in range(-organism_size+1, organism_size,2):
                        newrecord = self.create_pair_boxes(ortholog_rec, organism_size, i)
                        newrecord.bgcolor = organism_colors[organisms[n]]
                        ortholog_rec.graphics.append(newrecord)
                        n+=1
                    self.create_tile_box(ortholog_rec)

                if organism_size% 2 != 0:
                    n= 0
                    for i in range(int(-(organism_size-1)/2),int((organism_size+1)/2)):
                        newrecord = self.create_odd_boxes(ortholog_rec, organism_size, i)
                        newrecord.bgcolor = organism_colors[organisms[n]]
                        ortholog_rec.graphics.append(newrecord)
                        n += 1
                    self.create_tile_box(ortholog_rec)

    def pathway_genes_diferencial(self, dataframe, scale = 100, maxshared = 5, colormap = ""):
        '''
        Given a dataframe with expression values and index corresponding to
        ortholog kegg id, multiples the numbers by the scale number and
        represents in a heatmap format the expression levels of each sample
        in each ortholog representation in the pathway map.
        :param dataframe: pandas dataframe with expression values in each collumn
        and index corresponding to ortholog kegg id
        :param scale: int of number to be used to scale the values of dataframe
        :param maxshared: int with maximum number of boxes to draw in each
        ortholog representation in pathway map
        :param colormap: str with name of matplotlib colormap to be used
        :return:
        '''
        coverted_KEGG_ID, orthologs_ID = self.convert_KEGGID_koID(dataframe.index.tolist())
        dataframe = dataframe.ix[coverted_KEGG_ID] # filtrar kegg_IDs com KO
        orthologs_dic = self.ortholog_dic()

        new_index = []
        for ortholog in orthologs_ID:
            if ortholog in orthologs_dic.keys():
                new_index.append(orthologs_dic[ortholog])

        dataframe.index = new_index  # Passar os kegg_ID para Ko
        #TODO verificar este groupby com sum
        dataframe = dataframe.groupby(level=0).sum()

        value_range = (dataframe.max().max() - dataframe.min().min()) * scale
        #dataframe = dataframe.apply(self._scale_values, args=(scale,))

        if colormap == "":
            colormap = cm.get_cmap("coolwarm",value_range)
        else:
            try:
                colormap = cm.get_cmap(colormap, value_range)
            except:
                print("Colormap doesn't exist")
        print(dataframe)
        dataframe = dataframe.apply(self._conv_value_rgb, args=(colormap, scale))
        dataframe = dataframe.apply(self._conv_rgb_hex)

        nrboxes = len(dataframe.columns.tolist())
        for ortholog_rec in self.pathway.orthologs:
            ortholog_name = ortholog_rec.name.split(" ")[0]
            if ortholog_name in dataframe.index.tolist():
                colors = dataframe.ix[ortholog_name].tolist()
                if nrboxes > maxshared:
                    nrboxes = maxshared
                if nrboxes%2==0:
                    n = 0
                    for i in range(-nrboxes+1, nrboxes, 2):
                        newrecord = self.create_pair_boxes(ortholog_rec, nrboxes, i)
                        newrecord.bgcolor = colors[n]
                        ortholog_rec.graphics.append(newrecord)
                        n+=1
                    self.create_tile_box(ortholog_rec)

                if nrboxes% 2 != 0:
                    n= 0
                    for i in range(int(-(nrboxes-1)/2),int((nrboxes+1)/2)):
                        newrecord = self.create_odd_boxes(ortholog_rec, nrboxes, i)
                        newrecord.bgcolor = colors[n]
                        ortholog_rec.graphics.append(newrecord)
                        n += 1
                    self.create_tile_box(ortholog_rec)


    def pathway_box_list(self, dic_box_items, items, maxshared = 10, color = []):
        '''
        Represents items in the pathway map
        :param dic_box_items: dic with keys correspoding to int index of the list
        of ortholog elements in pathway and values correspoding to to elements
        of items
        :param items: list of items to be represented and given a specific color
        :param maxshared: int representing the maximum number of boxes to be
        represented in each element in the pathway map
        :param color: list of costum colors to be used to color the elements
        '''
        colors = self._set_colors(color, len(items))
        dic_colors = {}
        for i in range(0, len(items)):
            dic_colors[items[i]] = colors[i]

        for box in dic_box_items.keys():
            number_types = len(dic_box_items[box])
            if number_types > maxshared:
                number_types = maxshared
            if number_types % 2 == 0:
                n = 0
                for i in range(-number_types + 1, number_types, 2):
                    newrecord = self.create_pair_boxes(self.pathway.orthologs[box], number_types, i)
                    newrecord.bgcolor = dic_colors[dic_box_items[box][n]]
                    self.pathway.orthologs[box].graphics.append(newrecord)
                    n += 1
                self.create_tile_box(self.pathway.orthologs[box])

            if number_types % 2 != 0:
                n = 0
                for i in range(int(-(number_types - 1) / 2),
                               int((number_types + 1) / 2)):
                    newrecord = self.create_odd_boxes(self.pathway.orthologs[box], number_types, i)
                    newrecord.bgcolor = dic_colors[dic_box_items[box][n]]
                    self.pathway.orthologs[box].graphics.append(newrecord)
                    n += 1
                self.create_tile_box(self.pathway.orthologs[box])


    def pathway_boxes_diferencial(self, dataframe, log = False, colormap=""):
        '''
        Represents expression values present in a dataframe in the
        pathway map
        :param dataframe: pandas DataFrame with each column representing a sample
        and index corresponding to int list index of the ortholog elment in the
        pathway
        :param log: bol providing the option for a log normalization of data
        :param colormap: str representing a costum matplotlib colormap to be used
        '''
        dataframe = dataframe.groupby(level=0).sum()
        print("a")
        print(dataframe)
        print("b")
        print(dataframe.min().min())
        print(dataframe.max().max())

        if log:
            norm = cm.colors.LogNorm(vmin=dataframe.min().min(), vmax=dataframe.max().max())
        else:
            norm = cm.colors.Normalize(vmin=dataframe.min().min(), vmax=dataframe.max().max())

        print(norm(26739))
        print(norm(83316))
        print(norm(53316))
        print(norm(0))

        if colormap == "":
            colormap = cm.get_cmap("coolwarm")
        else:
            try:
                colormap = cm.get_cmap(colormap)
            except:
                print("Colormap doesn't exist")
        dataframe = dataframe.apply(self._conv_value_rgb, args=(colormap, norm))
        dataframe = dataframe.apply(self._conv_rgb_hex)



        nrboxes = len(dataframe.columns.tolist())
        for box in dataframe.index.tolist():
            colors = dataframe.ix[box].tolist()

            if nrboxes % 2 == 0:
                n = 0

                for i in range(-nrboxes + 1, nrboxes, 2):
                    newrecord = self.create_pair_boxes(self.pathway.orthologs[box], nrboxes, i)
                    newrecord.bgcolor = colors[n]

                    self.pathway.orthologs[box].graphics.append(newrecord)
                    n += 1
                self.create_tile_box(self.pathway.orthologs[box])

            if nrboxes % 2 != 0:
                n = 0
                for i in range(int(-(nrboxes - 1) / 2), int((nrboxes + 1) / 2)):
                    newrecord = self.create_odd_boxes(self.pathway.orthologs[box], nrboxes, i)
                    newrecord.bgcolor = colors[n]
                    self.pathway.orthologs[box].graphics.append(newrecord)
                    n += 1
                self.create_tile_box(self.pathway.orthologs[box])

################################################################################
################################################################################

class MoscaData():
    def __init__(self, file, pathway=[]):
        self.data = self.load_file(file)
        if len(pathway) > 0:
            self.pathway = pathway
        else:
            self.pathway = self.maps().keys()


    def load_file(self, file):
        '''
        Load MOSCA .csv file and prepares it for operations
        :param file: .csv file from MOSCA analysis
        return: pandas DataFrame object containing the formated .csv data
        '''
        df = pd.read_csv(file, delimiter=",", index_col=["Cross-reference (KEGG)"])
        df = df[df.index.notnull()]
        coverted_KEGG_ID, orthologs_ID = self.convert_KEGGID_koID(df.index.tolist())
        for i in range(0, len(coverted_KEGG_ID)):
            coverted_KEGG_ID[i] = str(coverted_KEGG_ID[i]) + ";"
        df = df.loc[coverted_KEGG_ID]
        df["orthologs"] = orthologs_ID
        df = df.set_index("orthologs")

        return df

    def maps(self):
        '''
        Creates a dic with all specific kegg pathway maps and their description
        :return: dic with all specific kegg pathway maps
        '''
        pathways = kegg_list("pathway")
        dic_pathway = {}
        generalmaps = ["01100", "01110", "01120", "01130", "01200", "01210",
                       "01212", "01212", "01220"]
        for i in pathways.read().split("\n")[:-1]:
            pathway = i.split("\t")[0][8:]
            desc = i.split("\t")[1]
            if pathway not in generalmaps:
                dic_pathway[pathway] = desc


        return dic_pathway

    def top_gemus(self, samples, n=10):
        '''
        Calculates top genus from samples
        :param samples: list of samples
        :param n: number of top genus to return
        :return: list of top genus
        '''
        data = self.data.groupby("Taxonomic lineage (GENUS)")[samples].sum()
        data["sums"] = data.sum(axis=1)
        data = data.sort_values(by=["sums"], ascending=False)
        if n > len(data.index.tolist()):
            n = len(data.index.tolist())
        return data.index.tolist()[:n - 1]


    ############################################################################
    ####                            Gets                                    ####
    ############################################################################

    def get_df(self):
        '''
        Gets loaded MOSCA data in pandas dataframe format
        :return: pandas dataframe with MOSCA data
        '''
        return self.data

    ############################################################################
    ####                            Sets                                    ####
    ############################################################################

    def set_pathway(self, pathway):
        '''
        Sets a pathway
        :param pathway:  KEGG pathway ID to be set
        '''
        self.pathway = pathway

    ############################################################################
    ####                              Helper                                ####
    ############################################################################

    def convert_KEGGID_koID(self, KEGG_ID):
        '''
        Converts KEGG_ID genes to Ortholog KO ID from KEGG
        :param KEGG_ID: (list) - KEGG ID genes
        :return: (tupple) - (list,list) - KEGG ID genes converted and ko IDs
        '''
        #TODO Adicionar excecoes para gerir KEGG IDs incorretos
        coverted_KEGG_ID = []
        orthologs_ID = []
        print(len(KEGG_ID))
        if len(KEGG_ID) <= 150:
            links = kegg_link("ko", KEGG_ID).read().split("\n")
            for link in links:
                if len(link) > 1:
                    coverted_KEGG_ID.append(link.split("\t")[0])
                    orthologs_ID.append(link.split("\t")[1][3:])
        else:
            i = 0
            while ((i+150) < len(KEGG_ID)):
                print(i)
                links = kegg_link("ko", KEGG_ID[i:i+150]).read().split("\n")
                for link in links:
                    if len(link) > 1:
                        coverted_KEGG_ID.append(link.split("\t")[0])
                        orthologs_ID.append(link.split("\t")[1][3:])
                i += 150
            links = kegg_link("ko", KEGG_ID[i:len(KEGG_ID)]).read().split("\n")
            for link in links:
                if len(link) > 1:
                    coverted_KEGG_ID.append(link.split("\t")[0])
                    orthologs_ID.append(link.split("\t")[1][3:])
        print(len(coverted_KEGG_ID), len(orthologs_ID))

        return (coverted_KEGG_ID, orthologs_ID)


    def set_colors(self, colors = [], ncolor = 1):
        '''
        Creates list of hex colors to be used,
        using matplotlib or using costum colors
        :param colors: list of hex colors
        :param ncolor: int indicating the ammount of hex colors should be created
        :return: returns list with hex color codes
        '''
        ret_color = []
        if len(colors) == 0:
            # if no colors are given creates a list of hex colors with ncolor
            # ammount from matplotlib discrete colormaps. If ncolor > 12
            # a continuous colormap is used instead.
            if ncolor <= 8:
                pastel2 = cm.get_cmap('Pastel2', 8)
                for i in range(ncolor):
                    ret_color.append(to_hex(pastel2(i)))
                return ret_color
            elif ncolor <= 12:
                set3 = cm.get_cmap("Set3", 12)
                for i in range(ncolor):
                    ret_color.append(to_hex(set3(i)))
                return ret_color
            else:
                rainbow = cm.get_cmap("rainbow",ncolor)
                for i in range(ncolor):
                    ret_color.append(to_hex(rainbow(i)))
                return ret_color
        else:
            # validates hex values and returns the original list
            isvalidhex = True
            for hexvalue in colors:
                match = re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', str(hexvalue))
                if not match:
                    isvalidhex = False
            if isvalidhex:
                return colors
            else:
                raise Exception("Colors aren't valid hex codes")

    def load_pathway(self, pathway):
        '''
        Loads a given pathway map in Pathway object
        :param pathway: str of kegg map pathway
        :return: Pathway map object
        '''
        return Pathway(pathway)


    def pathway_reset(self):
        '''
        resets the Pathway object
        '''
        self.pathway.reset_pathway()

    ############################################################################
    ####                          Operations                                ####
    ############################################################################

    def genomic_potential_genus(self, samples, genus, nr = 10):
        '''
        Represents the genomic potential of the dataset, by coloring each
        genus with an unique color. Only the orthologs with at least one
        sample expression level above the threshold will be represented
        :param samples:list of str column names of the dataset correspoding to
        expression values
        :param genus: list of genus to represent
        :param threshold: int representing the expression threshold
        '''
        #KO for each genus type
        if len(genus) > 0:
            genus = list(set(genus))
        else:
            genus = self.top_gemus(samples, nr)
        colors = self.set_colors([], len(genus))
        for map in self.pathway:
            pathway = self.load_pathway(map)
            dic_boxes = {}
            present_genus = []
            for type in genus:
                df = self.data[self.data["Taxonomic lineage (GENUS)"] == type][samples]
                df = df[df.any(axis=1)]
                orthologs = df.index.tolist()
                for ortholog in orthologs:
                    if ortholog in pathway.get_ko_boxes().keys():
                        if type not in present_genus:
                            present_genus.append(type)
                        for box in pathway.get_ko_boxes()[ortholog]:
                            if box in dic_boxes.keys():
                                dic_boxes[box].append(type)
                            else:
                                dic_boxes[box] = [type]
            pathway.pathway_box_list(dic_boxes, present_genus, 10, colors)
            name_pdf = str(map) + "_" + "genomic_potential"
            pathway.pathway_pdf(name_pdf)
            print("Map saved to " + name_pdf + ".pdf")


    def genomic_potential_sample(self, samples, threshold = 0):
        '''
        Represents the genomic potential of the dataset, by coloring each sample
        with an unique color. Only orthologs iwth at least one sample expression
        level above the theshold will be represented
        :param samples: ist of str column names of the dataset correspoding to
        expression values
        :param threshold: threshold: int representing the expression threshold
        '''
        dic_box_sample = {}
        df = self.data[samples]
        for sample in samples:
            orthologs = df.loc[df[sample] > threshold].index.tolist()
            for ortholog in orthologs:
                if ortholog in self.pathway.get_ko_boxes().keys():
                    for box in self.pathway.get_ko_boxes()[ortholog]:
                        if box in dic_box_sample.keys():
                            dic_box_sample[box].append(sample)
                            dic_box_sample[box] = list(set(dic_box_sample[box]))
                        else:
                            dic_box_sample[box] = [sample]

        self.pathway.pathway_box_list(dic_box_sample, samples)


    def diferential_expression_sample(self, samples, log=False):
        '''
        Represents in small heatmaps the expression levels of each sample on the
        dataset present in the given pathway map. The values can be transford to
        a log10 scale
        :param samples: ist of str column names of the dataset correspoding to
        expression values
        :param log: bol providing the option for a log normalization of data
        '''
        for map in self.pathway:
            pathway = self.load_pathway(map)
            data = []
            new_index = []
            df = self.data.groupby(level = 0)[samples].sum()
            df = df[df.any(axis=1)]
            orthologs = df.index.tolist()
            for ortholog in orthologs:
                if ortholog in pathway.get_ko_boxes().keys():
                    for box in pathway.get_ko_boxes()[ortholog]:
                        for ko in df.loc[orthologs].index.tolist():
                            data.append(df.loc[ko].tolist())
                            new_index.append(box)
            new_df = pd.DataFrame(data, columns=samples, index=new_index)
            new_df = new_df[new_df.any(axis=1)]
            pathway.pathway_boxes_diferencial(new_df, log)

            if log:
                name_pdf = str(map) + "_" + "differential_expression" + "_" + "log"
                pathway.pathway_pdf(name_pdf)
            else:
                name_pdf = str(map) + "_" + "differential_expression"
                pathway.pathway_pdf(name_pdf)
            print("Map saved to " + name_pdf + ".pdf")


############################################################################
####                            Test                                    ####
############################################################################

def test_clean_pathway():
    clean_pathway = Pathway("map00680")
    clean_pathway.pathway_pdf("test_clean_pathway")

def test_pathway_genes():
    test_pathway = Pathway("map00680")
    kegg_IDs = ["pmy:Pmen_3467", "pre:PCA10_14290", "mba:Mbar_A0841", "psub:Pelsub_P2997", "mby:MSBRM_0151", "mby:MSBRM_0158"]
    test_pathway.pathway_genes(kegg_IDs)
    test_pathway.pathway_pdf("test_pathway_genes")
    test_pathway.reset_pathway()
    kegg_IDs.append("ddh:Desde_0778")
    test_pathway.pathway_genes(kegg_IDs, ["#ef0000"])
    test_pathway.pathway_pdf("test_pathway_genes_vermelho")

def test_pathway_organismo():
    test_pathway = Pathway("map00680")
    test_pathway.pathway_organismo("eco")
    test_pathway.pathway_pdf("test_pathway_organismo")
    test_pathway.reset_pathway()
    test_pathway.pathway_organismo("eco",["#ef0000"])
    test_pathway.pathway_pdf("test_pathway_organismo_vermelho")

def test_pathway_genes_organismo():
    test_pathway = Pathway("map00680")
    kegg_IDs = ["ppg:PputGB1_2248","ppt:PPS_3154","ppx:T1E_2051","ppun:PP4_21640","ppud:DW66_3426","sme:SMc02610", "mby:MSBRM_0152"]
    test_pathway.pathway_genes_organismo(kegg_IDs)
    test_pathway.pathway_pdf("test_pathway_genes_organismo")

def test_pathway_genes_diferencial():
    df_pair = pd.DataFrame(
        [[1, 2, 3, 4, 5], [2, 3, 4, 5, 6], [3, 4, 5, 6, 7], [4, 5, 6, 7, 8], [5, 6, 7, 8, 9]],
        columns=['c1', 'c2', 'c3', 'c4', 'c5'],
        index=["pmy:Pmen_3467","pre:PCA10_14290","mba:Mbar_A0841","psub:Pelsub_P2997","mby:MSBRM_0151"])

    df_odd = pd.DataFrame(
        [[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6], [4, 5, 6, 7],
         [5, 6, 7, 8]],
        columns=['c1', 'c2', 'c3', 'c4'],
        index=["pmy:Pmen_3467", "pre:PCA10_14290", "mba:Mbar_A0841",
               "psub:Pelsub_P2997", "mby:MSBRM_0151"])

    def_nine = pd.DataFrame(
        [[1, 2, 5, 15, 17, 19, 45, 55, 50], [2, 4, 3, 20, 19, 22, 50, 52, 70], [2, 1, 2, 20, 15, 24, 50, 57, 48], [2, 5, 10, 22, 25, 20, 60, 61, 65],
         [1, 3, 4, 20, 15, 19, 50, 66, 48]],
        columns=['c1', 'c2', 'c3', 'c4', "c5", "c6", "c7", "c8", "c9"],
        index=["pmy:Pmen_3467", "pre:PCA10_14290", "mba:Mbar_A0841",
               "psub:Pelsub_P2997", "mby:MSBRM_0151"])

    def_nine = pd.DataFrame(
        [[1, 2, 5, 15, 17, 19, 45, 55, 50], [2, 4, 3, 20, 19, 22, 50, 52, 70],
         [2, 1, 2, 20, 15, 24, 50, 57, 48], [2, 5, 10, 22, 25, 20, 60, 61, 65]],
        columns=['c1', 'c2', 'c3', 'c4', "c5", "c6", "c7", "c8", "c9"],
        index=[ "pre:PCA10_14290", "mba:Mbar_A0841",
               "psub:Pelsub_P2997","mby:MSBRM_0151"])

    print(def_nine)
    test_pathway = Pathway("map00680")
    #test_pathway.pathway_genes_diferencial(df_pair)
    #test_pathway.pathway_genes_diferencial(df_odd)
    test_pathway.pathway_genes_diferencial(def_nine, 100, 9)
    test_pathway.pathway_pdf("test_pathway_genes_diferencial")

def test_pathway_genes_amostras():
    test_pathway = Pathway("map00680")
    kegg_IDs_matrix = [["ppg:PputGB1_2248", "ppt:PPS_3154", "ppx:T1E_2051",
                "ppun:PP4_21640", "ppud:DW66_3426", "sme:SMc02610",
                "mby:MSBRM_0152"]]
    test_pathway.pathway_genes_amostras(kegg_IDs_matrix)
    test_pathway.pathway_pdf("test_pathway_genes_amostras")


def test_genomic_potential_genus():
    data = MoscaData("all_info_normalized.csv",["map00680"])
    samples_mg = ["grinder-reads0.17a", "grinder-reads0.17b", "grinder-reads0.17c"]
    data.genomic_potential_genus(samples_mg, ["Pseudomonas", "Methanosarcina","Pelolinea","Thermoplasmatales"])


def test_genomic_potential_sample():
    data = MoscaData("all_info_normalized.csv", ["map00680"])
    samples_mg = ["grinder-reads0.17a", "grinder-reads0.17b","grinder-reads0.17c"]
    data.genomic_potential_sample(samples_mg)
    data.pathway_pdf("test_genomic_potential_sample")


def test_diferential_expression_sample():
    data = MoscaData("all_info_normalized_all.csv")
    #data.set_pathway(["map00680"])

    samples_mt = ["grinder-reads1a", "grinder-reads1b", "grinder-reads1c", "grinder-reads3a", "grinder-reads3b", "grinder-reads3c"]


    data.diferential_expression_sample(samples_mt, True)
    data.diferential_expression_sample(samples_mt, False)


############################################################################
####                            Main                                    ####
############################################################################

if __name__ == '__main__':
    ## Test Functions
    #test_clean_pathway()
    #test_pathway_genes()
    #test_pathway_organismo()
    #test_pathway_genes_organismo()
    #test_pathway_genes_diferencial()
    #test_pathway_genes_amostras()
    test_genomic_potential_genus()
    #test_genomic_potential_sample()
    #test_diferential_expression_sample()
    import matplotlib as mpl


    '''
    import matplotlib.pyplot as plt
    from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
                                      AnnotationBbox)
    from matplotlib.cbook import get_sample_data
    from PIL import Image
    
    im = plt.imread(get_sample_data('/home/nowinder/Projects/Uminho/Projeto/test_pathway.png'))

    fig, ax = plt.subplots()
    ax.plot(range(10))
    newax = fig.add_axes([0, 0, 1, 1], anchor='NE', zorder=1)
    newax.imshow(im)
    newax.axis('off')
    plt.savefig('myfig.png', dpi=5000
    plt.show()


    min, max = (-40, 30)
    step = 10
    mymap = cm.get_cmap("rainbow")
    Z = [[0, 0], [0, 0]]
    levels = range(min, max + step, step)

    CS3 = plt.contourf(Z, levels, cmap=mymap)
    plt.clf()
    plt.colorbar(CS3)  # using the colorbar info I got from contourf
    plt.axis('off')
    plt.show()
    # plt.figure(figsize=(1023.36 / 100, 670.8 / 100), dpi=10)

    plt.savefig('myfig.pdf')    
    '''


