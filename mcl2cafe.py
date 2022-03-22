#! /usr/bin/env python3
"""
mcl2cafe.py: Format mcl output to be used by cafe
Author: Zhipeng Qu
Data: 28/01/2022
"""

# import modules
import argparse
import re
import os

# Get arguments. Three arguments will be required:
# falist: list of fasta files
# mcl: input mcl output file
# cafe: raw table can be used as CAFE input

parser = argparse.ArgumentParser(description = 'Format mcl output to CAFE input')
parser.add_argument('falist', help = 'A file including all fasta files, with file name prefix as species name.')
parser.add_argument('mcl', help = 'Dump output of mcl.')
parser.add_argument('cafe', help = 'Table used as CAFE input.')
parser.add_argument('cafe_large', help = 'Table used as CAFE input with big gene count')
parser.add_argument('cafe_small', help = 'Table used as CAFE input with big gene count filtered')


args = parser.parse_args()
falist = args.falist
mcl_file = args.mcl
cafe_file = args.cafe
cafe_large_file = args.cafe_large
cafe_small_file = args.cafe_small

# Define fa_parser function
def fa_parser(fa):
    """
    Function fa_parser will read fasta file and store it in fasta class
    Args:
        fa (file name): fasta format file
    Returns:
        fa_cls (dict): fasta class
    """

    fa_dict = {} # use dictory to refer sequence based on header
    with open(fa, 'r') as fasta:
        header = ''
        sequence = ''

        for line in fasta:
            line = line.rstrip() # remove last character
            if line.startswith('>'):
                if header != "":
                    fa_dict[header] = sequence
                header = line.lstrip('>') # remove '>' in header
                sequence = ''
            else: # concatenate sequences if they are in multiple lines
                sequence += line

        # add last sequence into dictory
        fa_dict[header] = sequence
    return fa_dict

# get fasta files for all species into a table
fas_header = {} # initialize the dict to store all headers
sp_list = [] # initialize species list
with open(falist, 'r') as fas:
    for single_fa in fas: 
        single_fa = single_fa.rstrip()
        sp = os.path.splitext(single_fa)[0] # get sp name from file name
        sp_list.append(sp)

        fas_dict = fa_parser(single_fa)
        for header in fas_dict:
            header_s = header.split(" ", 1)[0]
            fas_header[header_s] = sp
            out_line = [header_s, sp]

# process mcl dump output
with open(mcl_file, 'r') as mcl, open(cafe_file, 'w') as cafe, open(cafe_large_file, 'w') as cafe_large, open(cafe_small_file, 'w') as cafe_small:
    n = 0

    sps_header = '\t'.join(sp_id for sp_id in sp_list)
    cafe.write('Desp\tFamily\t' + sps_header + '\n')
    cafe_large.write('Desp\tFamily\t' + sps_header + '\n')
    cafe_small.write('Desp\tFamily\t' + sps_header + '\n')
    for ortho in mcl:
        n += 1
        ortho = ortho.rstrip()
        ortho_genes = ortho.split('\t')
        sp_ids_dict = dict((sp_id, 0) for sp_id in sp_list) # init count = 0 for sp dict
        for gene_id in ortho_genes:
            for sp_id in sp_list:
                if fas_header[gene_id] == sp_id:
                    sp_ids_dict[sp_id] += 1

        # write raw cafe input and get count info
        cafe.write('ortho' + str(n) + '\t' + 'ortho' + str(n) + '\t')
        clade_count = 0
        size_filter = True
        for sp_id in sp_list:
            cafe.write(str(sp_ids_dict[sp_id]) + '\t')
            if sp_ids_dict[sp_id] > 0:
                clade_count += 1
            if sp_ids_dict[sp_id] >= 100:
                size_filter = False
        cafe.write("\n")

        # filter cafe input based on clade presence and gene numbers
        if clade_count >= 2 and size_filter:
            cafe_small.write('ortho' + str(n) + '\t' + 'ortho' + str(n) + '\t')
            for sp_id in sp_list:
                cafe_small.write(str(sp_ids_dict[sp_id]) + '\t')
            cafe_small.write('\n')
        else:
            cafe_large.write('ortho' + str(n) + '\t' + 'ortho' + str(n) + '\t')
            for sp_id in sp_list:
                cafe_large.write(str(sp_ids_dict[sp_id]) + '\t')
            cafe_large.write('\n')
