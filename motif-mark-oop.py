#!/usr/bin/env python
 
import re
import argparse
import cairo
import math
import random
import bioinfo 


### Argparse
def get_args():
    parser = argparse.ArgumentParser(description="Program for Motif Mark")
    parser.add_argument('-f', help='fasta file')
    parser.add_argument('-m', help='motifs file')
    parser.add_argument('-o', help='output file')
    return parser.parse_args()

args = get_args()
f = args.f
m = args.m
o = args.o

# Global dictionary to store and reference read info {header:sequence}
sequence_dictionary:dict = {}


### Classes
class Sequence:
    ''' Reads a fasta file and returns a list of gene headers and sequences and extracts the genes and exons from the sequence'''
    def __init__(self, oneline_file:str, sequence_dictionary:dict):
       
        def _load_oneline():
            
            headers:list = []
            sequences:list = []

            with open(oneline_file, 'r') as oneline:
                for line in oneline.readlines():
                    line = line.strip('\n')
                    if line[0] == '>':
                        header = line.split()[0]
                        header = header.strip('>')
                        headers.append(header)
                    else:
                        seq = line
                        sequences.append(seq)
                        sequence_dictionary[header] = (seq)

            return headers, sequences

        self.oneline_file = oneline_file
        self.headers, self.sequences = _load_oneline()

    def get_gene(self, sequence_dictionary:dict) -> list:
        '''returns list of gene sequences from the dictionary'''
        gene_sequences = []
        for i in sequence_dictionary:
            gene_sequences.append(sequence_dictionary[i])
        
        return gene_sequences

    def get_exon(self, sequence_dictionary:dict) -> dict:
        '''returns exon dict containing exon begin and end loc {exon_seq:span}'''
        exon_dict = {}
        for key, val in sequence_dictionary.items():
            iterate = re.finditer("[A-Z]+", val)
            for match in iterate:
                exon_dict[key] = match.span()
               
        return exon_dict 
    

class Motif:
    def __init__(self, motif_file:str):
        '''Reads a motif file and returns a dictionary with searchable version of each motif {unambiguous:motif}'''
        def _load_motif():
            motifs:list = []

            with open(motif_file, 'r') as motif:
                for line in motif.readlines():
                    line = line.strip('\n').upper()
                    motifs.append(line)
                
            return motifs 
        
        self.motif_file = motif_file
        self.motifs = _load_motif()
    
    def parse_motifs(self) -> list:
        '''Takes the sequence_dictionary as an argument and returns a list of the locations of matching motifs'''
        motif_locs:list = []

        for val in sequence_dictionary.values():
            gene_seq = val
            converted_motifs:list = []
            tmp:list = []
            for line in self.motifs:
                conv_motif = ""     # empty string that will hold converted motifs
                for raw_motif in line:
                    conv_motif += str(bioinfo.IUPAC_bases[raw_motif])    # convert to unambiguouos IUPAC terms
                converted_motifs.append(conv_motif)
                for i in converted_motifs:
                    matches = re.finditer(i, gene_seq)   # find matches based off list of converted motif seqs   
                    match_locs:list = []
                    for match_seq in matches:
                        motif_loc = match_seq.span()    # find location of matching motif
                        match_locs.append(motif_loc)
                tmp.append(match_locs)
            motif_locs.append(tmp)

        return motif_locs
    

### Functions
def assign_colors():
    ''' Set values for color input'''
    motif_colors:dict = {}
    for i, motif in enumerate(range(50)):
        x,y,z = random.random(), random.random(), random.random()
        motif_colors[i] = [x,y,z]
    
    return motif_colors

def longest_gene() -> int:
    '''Find the longest gene length for the plot width'''
    longest_seq:int = 0
    for seq in sequence_dictionary.values():
        if len(seq) > longest_seq:
            longest_seq = len(seq)
    
    return longest_seq 

def oneline_fasta(file):
    '''Reads a fasta file and writes to a new fasta file with the sequence line from each record contained in one line'''
    output_file = ""
    with open(o, 'w') as out:                                    
        with open(file, 'r') as fa:
            for i, line in enumerate(fa):
                line = line.strip('\n')
                if i == 0:
                    print(line, file=out)
                elif line[0] == '>' and i != 0:
                    print(file=out)
                    print(line, file=out)
                else:
                    print(line, end='', file=out)
            print(file=out)


if __name__ == "__main__":
    oneline_fasta(f)
    print(f"\t> All fasta record sequences converted to oneline")

    # Instantiate classes
    sequences = Sequence('oneline.fasta', sequence_dictionary)
    genes = sequences.get_gene(sequence_dictionary)
    exons = sequences.get_exon(sequence_dictionary)
    
    motifs = Motif(m)
    motif_locations = motifs.parse_motifs()
    motifs_list = motifs.motifs
    
    colors = assign_colors()

    # PyCairo set up
    width_value = longest_gene()
    WIDTH = int((width_value)*1.25)
    HEIGHT = int(len(sequence_dictionary.values()) * 115) + 40

    # create the coordinates display graphic, designate output
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
    # create the coordinates you will be drawing on (like a transparency)
    image = cairo.Context(surface)

    #set vertical and horizontal positons 
    x_start = 25
    y_start = 40

    motif_count = 0
    gene_count = 0

    for gene_sequence in sequence_dictionary.keys():

        image.set_line_width(7)
        image.set_source_rgb(.3,.3,.3)      # set gray scale
        image.move_to(x_start,y_start+25)
        image.line_to(x_start+len(genes[gene_count]), y_start+25)       # use gene sequence length to define line
        image.stroke()

        head = list(sequence_dictionary.keys())[gene_count]     # pull gene name from headers in dictionary for labels
        image.move_to(x_start + 10, (y_start) - 5)      
        image.set_source_rgb(.8,.7,.8) 
        image.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        image.set_font_size(17)
        image.show_text(head) 
        image.stroke()

        exon_cord = exons[gene_sequence]        # get exon pos
        exon_cord = list(exon_cord)
        image.rectangle(exon_cord[0], y_start, exon_cord[1] - exon_cord[0],50) 
        image.set_source_rgb(0.8,0.8,0.8)       # set gray scale
        image.fill()

        color_count = 1
        for motif in motif_locations[gene_count]:
            for pos in motif:       # iterate through motif locations

                x,y,z = colors[color_count]     # random colors for motifs 

                image.set_line_width(5)
                image.set_source_rgb(x,y,z)     # set gray scale
                image.rectangle(pos[0] + 25, y_start,pos[1] - pos[0], 50) # use motif locations to place vertical motif markers
                image.fill()
            color_count += 1

        y_start += 100

        gene_count += 1
        motif_count += 1

    # Key/labeling
    motif_count = 1
    vert = 425
    for motif in motifs_list:

        image.set_source_rgb(.8,.7,.8)      # set gray scale
        image.move_to(25, vert + 5)
        image.show_text(motif)                                                        
        image.stroke()

        x,y,z = colors[motif_count]         # random colors for motifs

        image.set_source_rgb(x,y,z)         # set gray scale
        image.rectangle(WIDTH - 750, vert - 10 , 10, 10) 
        image.fill()

        motif_count += 1
        vert += 22

    surface.write_to_png("Figure_1.png") 
    surface.finish()
    print(f"\t> Image Completed")