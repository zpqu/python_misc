#! /usr/bin/env python3
"""
format_links.py: format collinear links for WGD to bed format used by circlize
Author: Zhipeng Qu
Date: 21/03/2022
"""

# import modules
import argparse
import re
import os

# Get arguments, three arguments will be required:
# gff: pseodu gff gene list
# colli: collinear file
# bed1: output source bed file
# bed2: output desitination bed file

parser = argparse.ArgumentParser(description = 'Format collinear liks for WGD to bed format used by circlize')
parser.add_argument('gff', help = 'pseodu gff gene list')
parser.add_argument('colli', help = "collinear file from WGD analysis")
parser.add_argument('bed1', help = 'output bed1 file')
parser.add_argument('bed2', help = 'output bed2 file')

args = parser.parse_args()
gff = args.gff
colli = args.colli
bed1 = args.bed1
bed2 = args.bed2

# Get gene coordiantes
gene_chr = {} # init the dict to store gene chr
gene_start = {} # init the dict to store gene start
gene_end = {} # init the dict to store gene end
with open(gff, 'r') as gff_gene:
    for single_gff in gff_gene:
        single_gff = single_gff.rstrip()
        single_gene = single_gff.split('\t')
        single_chr = re.sub('Sfla-', '', single_gene[0]) 
        single_name = single_gene[1]
        single_start = int(single_gene[2])
        single_end = int(single_gene[3])
        gene_chr[single_name] = single_chr
        gene_start[single_name] = single_start
        gene_end[single_name] = single_end

# Process pseodu links
with open(colli, 'r') as colli_link, open(bed1, 'w') as bed1_out, open(bed2, 'w') as bed2_out:
    ct = 1
    n = 0
    source_start = 0
    source_end = 0
    target_start = 0
    target_end = 0
    block = ''
    strand = ''
    bed_list = []

    for line in colli_link:
        line = line.rstrip()
        if line.startswith('#') and (not line.startswith('## Alignment')):
            continue
        if line.startswith('## Alignment'):
            if ct > 1:
                bed_list.append(gene_chr[source] + '\t' + str(source_start) + '\t' + str(source_end) + '\t' + block + '\t' + gene_chr[target] + '\t' + str(target_start) + '\t' + str(target_end) + '\t' + block)
            block = re.sub(':.+', '', line)
            block = re.sub('## ', '', block)
            block = re.sub(' ', '', block)
            strand = re.sub('.+ ', '', line)
            ct = ct + 1
            n = 0
        else:
            elemts = line.split()
            source = elemts[-3]
            target = elemts[-2]
            #print(source + '\t' + target)
            if n == 0:
                source_start = gene_start[source]
                if strand == 'plus':
                    target_start = gene_start[target]
                else:
                    target_start = gene_end[target]
            source_end = gene_end[source]
            if strand == 'plus':
                target_end = gene_end[target]
            else:
                target_end = gene_start[target]
            n = n + 1
    bed_list.append(gene_chr[source] + '\t' + str(source_start) + '\t' + str(source_end) + '\t' + block + '\t' + gene_chr[target] + '\t' + str(target_start) + '\t' + str(target_end) + '\t' + block)

    for bed in bed_list:
        if bool(re.search('scaffold', bed)):
            continue
        #print(bed)
        beds = bed.split('\t')
        bed1_out.write(beds[0] + '\t' + beds[1] + '\t' + beds[2] + '\t' + beds[3] + '\n')
        bed2_out.write(beds[4] + '\t' + beds[5] + '\t' + beds[6] + '\t' + beds[7] + '\n')
