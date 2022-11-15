#! /usr/bin/env python3
"""
blastn_coverage.py: Calculate the coverage of blastn alignment in query and subject.
Author: Zhipeng Qu
Date: 30/11/2021
"""

# import modules
import argparse
import re

# Get arguments. four arguments will be requrired:
# query: the query fasta file used to calculate query length
# subject: the subject fasta file used to calcualte subject length
# input: tab delimited file from blast result using -outfmt 6
# output: tab delimited file with query and subject coverage information added
#         for all hits
parser = argparse.ArgumentParser(description = "Calculate the blast alignment coverage")
parser.add_argument('query', help = 'Fasta file including all query sequences.')
parser.add_argument('subject', help = 'Fasta file including all subject sequences.')
parser.add_argument('input', help = 'Tab delimited file from blast.')
parser.add_argument('output', help = 'Tab delimited file, with coverage added.')

args = parser.parse_args()
query_fa = args.query
subject_fa = args.subject
in_name = args.input
out_name = args.output


# Define fasta_parser function 
def fasta_parser(fasta):
    """
    Function fasta_parser will read fasta file and put it in Fasta class
    Args: 
        fasta (file name): fasta format file
    Returns:
        fa_cls (dict): fasta class
    """

    fa_dict = {} # use dictory to refer sequence based on header
    with open(fasta, 'r') as fa:
        header = ""
        sequence = ""

        for line in fa:
            line = line.rstrip() # remove last character
            if line.startswith('>'):
                if header != "":
                    header = header.split(' ', 1)[0] # capture non-space as header
                    fa_dict[header] = sequence
                header = line.lstrip('>') # remove '>' in header
                sequence = ""
            else: # concatenate sequences if they are in multiple lines
                sequence += line

        # add last sequence into dictory
        fa_dict[header] = sequence
    return fa_dict

# Import query and subject sequence files and put them in dictories
querys = fasta_parser(query_fa)
subjects = fasta_parser(subject_fa)

# Open input and out files for reading and writing. Two steps of processing 
# were included: Step 1, multiple hsps from the same subject hits will be 
# merged based on following criteria: gaps in query hsps < 50 (nt) and gaps 
# in subject hsps < 10 (aa). Each subject hit will have one final merged 
# hsps and stored in a list "merged_blast". Step 2, Calculate coverages of
# hits and put in output.
with open(in_name, 'r') as in_file, open(out_name, 'w') as out_file:

    # Step 1, merge hsps for the same subject
    merged_blast = [] # put merged subject into list
    qid_base = "" # initilize empty qid_base to make the test
    for line in in_file.readlines():
        line = line.rstrip()
        blast_list = line.split('\t')
        
        if qid_base == "": # get the 1st item
            qid_base = blast_list[0]
            sid_base = blast_list[1]

            # Get intervals in case they are not stranded
            qstart_base = min(int(blast_list[6]), int(blast_list[7]))
            qend_base = max(int(blast_list[6]), int(blast_list[7]))
            sstart_base = min(int(blast_list[8]), int(blast_list[9]))
            send_base = max(int(blast_list[8]), int(blast_list[9]))
            slen_base = send_base - sstart_base

        else: # processsing remaining items
            qid = blast_list[0]
            sid = blast_list[1]
            qstart = min(int(blast_list[6]), int(blast_list[7]))
            qend = max(int(blast_list[6]), int(blast_list[7]))
            sstart = min(int(blast_list[8]), int(blast_list[9]))
            send = max(int(blast_list[8]), int(blast_list[9]))
            slen = send - sstart
            
            # Merge hsps if their are multiple
            if qid == qid_base and sid == sid_base:
                if (qstart - qend_base) < 50 and (sstart - send_base) < 10:
                    qstart_base = min(qstart_base, qstart)
                    qend_base = max(qend_base, qend)
                    sstart_base = min(sstart_base, sstart)
                    send_base = max(send_base, send)
                
                # For hsps with large gaps, use the hsp with longest subject alignment
                else:
                    if slen > slen_base:
                        qstart_base = qstart
                        qend_base = qend
                        sstart_base = sstart
                        send_base = send

            else: # New item when cid changes
                line_items =[qid_base, sid_base, qstart_base, qend_base, sstart_base, send_base] 
                merged_blast.append(line_items)

                # Reset base values
                qid_base = qid
                sid_base = sid
                qstart_base = qstart
                qend_base = qend
                sstart_base = sstart
                send_base = send

    # Store the last item
    line_items = [qid_base, sid_base, qstart_base, qend_base, sstart_base, send_base] 
    merged_blast.append(line_items)
    
    # Step 2, find the longest subject alignment for each query
    qname_base = ""
    for hit_list in merged_blast:
        qname = hit_list[0]
        sname = hit_list[1]
        qStart = hit_list[2]
        qEnd = hit_list[3]
        qLen = qEnd - qStart
        sStart = hit_list[4]
        sEnd = hit_list[5]
        sLen = sEnd - sStart
        qCov = round(qLen/len(querys.get(qname))*100, 2)
        sCov = round(sLen/len(subjects.get(sname))*100, 2)

        out_line = [qname, sname, qStart, qEnd, qLen,
            qCov, sStart, sEnd, sLen, sCov]
        out_file.write('\t'.join(str(v) for v in out_line))
        out_file.write('\n')
