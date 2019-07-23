#!/usr/bin/env python

# sam2vcf.py ver. 1.1.1
# Copyright (C) 2019 Yuki Kato
# This script is used to convert SAM alignments for Bivartect into predicted VCF variants.
# The output will be shown in standard out.
# The variation class (VC) consists of:
# SNV: single nucleotide variant
# DIV: deletion/insertion variant
# SV: structural variant of length >= 50 bp
# BP: unassigned breakpoint
# Note: predicted variants on the reverse strand will appear on the forward strand in the VCF output.
# Usage: Type "./sam2vcf.py -h" in your terminal.

import argparse

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    usage='./sam2vcf.py [option]* <sam> <fasta>',
    description='sam2vcf.py ver. 1.1.1\n\nThis script is used to convert SAM alignments for Bivartect into predicted VCF variants.\nThe output will be shown in standard out.\n\nThe variation class (VC) consists of:\n\n  SNV: single nucleotide variant\n  DIV: deletion/insertion variant\n  SV: structural variant of length >= 50 bp\n  BP: unassigned breakpoint\n\nNote: predicted variants on the reverse strand will appear on the forward strand in the VCF output.'
)
parser.add_argument('sam', metavar='sam <STR>', type=str,
                    help='path to the SAM alignments for Bivartect')
parser.add_argument('fasta', metavar='fasta <STR>', type=str,
                    help='path to the (gzipped) reference FASTA sequence')
args = parser.parse_args()

import re
import gzip
import operator

def readFASTA(obj, num_lines):
    chrnum = 0
    seq = ""

    for i, line in enumerate(obj):
        if i != num_lines - 1:
            m2 = r2.match(line)

            if m2: # Header in FASTA
                if i != 0:
                    reference[chrnum] = seq # For the previous chromosome

                if m2.group(1) == "X" or m2.group(1) == "Y" or m2.group(1) == "MT":
                    chrnum = char_numbers[m2.group(1)] # Initialization
                
                else:
                    chrnum = int(m2.group(1)) # Initialization
                
                seq = "" # Initialization

            else:
                seq += line.strip()
        
        elif i == num_lines - 1:
            seq += line.strip()
            reference[chrnum] = seq # For the current chromosome

def revComp(seq):
    r_seq = seq[::-1]
    rc_seq = ""

    for i in range(len(r_seq)):
        if r_seq[i] == "A":
            rc_seq = rc_seq + "T"
        elif r_seq[i] == "C":
            rc_seq = rc_seq + "G"
        elif r_seq[i] == "G":
            rc_seq = rc_seq + "C"
        elif r_seq[i] == "T":
            rc_seq = rc_seq + "A"
        else:
            rc_seq = rc_seq + r_seq[i]
    
    return rc_seq

# Read the FASTA input
r1 = re.compile(r"\.gz$")
gz = False

if r1.search(args.fasta):
    gz = True

r2 = re.compile(r">(\w+)")
char_numbers = {"X":97, "Y":98, "MT":99} # Dictionary of {chars:int}
number_chars = {97:"X", 98:"Y", 99:"MT"} # Dictionary of {int:chars}
reference = {} # Dictionary of {chrnum:reference_seq}

if gz:
    num_lines = sum(1 for i in gzip.open(args.fasta, "rt"))

    with gzip.open(args.fasta, "rt") as f:
        readFASTA(f, num_lines)

else:
    num_lines = sum(1 for i in open(args.fasta, "r"))

    with open(args.fasta, "r") as f:
        readFASTA(f, num_lines)

# SAM file
# QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ
r3 = re.compile(r"@")
r4 = re.compile(r"(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\w+)")
r5 = re.compile(r"(\d+):(\w+):([\w\-\*]+):([\w\-\*]+)")
vcfs = {} # Dictonary of {(CHROM, POS, REF, ALT):[ID, QUAL, FILTER, INFO]}

# Read the SAM input
with open(args.sam, "r") as f:
    for line in f:
        m3 = r3.match(line)

        if not m3:
            m4 = r4.match(line)

            if m4:
                qname = m4.group(1)
                flag = int(m4.group(2))
                rname = m4.group(3)
                pos = int(m4.group(4))
                seq = m4.group(10)
                bp = 0
                id = 0
                ref = "."
                alt = "."
                info = "."
                
                # Case: unmapped
                if rname == "*":
                    continue
                
                # Case: unmapped
                #if flag == 4 or flag == 77 or flag == 141:
                    #bp = 0
                
                chrnum = 0

                if rname == "X" or rname == "Y" or rname == "MT":
                    chrnum = char_numbers[rname]
                
                else:
                    chrnum = int(rname)
                
                m5 = r5.match(qname)
                
                if m5:
                    # Parse a Bivartect output header
                    id = int(m5.group(1))
                    vc = m5.group(2)
                    ref = m5.group(3)
                    alt = m5.group(4)

                    # Set appropriate bp, ref, alt & info depending on VC & strand
                    # Case: SNV
                    if vc == "SNV":
                        info = "VC=SNV" # VC: variation class

                        if flag == 16: # Reverse strand
                            bp = pos - 1
                            ref = revComp(ref.upper())
                            alt = revComp(alt.upper())
                        
                        else: # Forward strand
                            bp = pos + len(seq)
                    
                    # Case: INS
                    elif ref == "-":
                        if len(alt) >= 50 or "*" in alt:
                            info = "VC=SV"
                        
                        else:
                            info = "VC=DIV" # DIV: deletion/insertion variant
                        
                        if flag == 16: # Reverse strand
                            bp = pos - 1 # pos is the last common base before the variant
                            ref = reference[chrnum][pos-2]
                            alt = reference[chrnum][pos-2] + revComp(alt.upper())

                        else: # Forward strand
                            bp = pos + len(seq) - 1 # '-1' means inclusion of the last common base before the variant
                            ref = seq[-1]
                            alt = seq[-1] + alt
                    
                    # Case: DEL
                    elif alt == "-":
                        if len(ref) >= 50 or "*" in ref:
                            info = "VC=SV"
                        
                        else:
                            info = "VC=DIV"
                        
                        if flag == 16: # Reverse strand
                            bp = pos - len(ref) - 1
                            start = pos - len(ref) - 2 # 0-based
                            ref = reference[chrnum][start] + revComp(ref.upper())
                            alt = reference[chrnum][start]

                            if "*" in ref:
                                info += ";IMPRECISE"
                        
                        else: # Forward strand
                            bp = pos + len(seq) - 1
                            ref = seq[-1] + ref
                            alt = seq[-1]

                    # Case: BP
                    else:
                        info = "VC=BP" # BP: unassigned breakpoint

                        if flag == 16: # Reverse strand
                            bp = pos - len(ref) - 1 # Best estimate
                            start = pos - len(ref) - 2 # 0-based
                            ref = reference[chrnum][start] + revComp(ref.upper())
                            alt = reference[chrnum][start] + revComp(alt.upper())
                            info += ";IMPRECISE"
                        
                        else: # Forward strand
                            bp = pos + len(seq) - 1
                            ref = seq[-1] + ref
                            alt = seq[-1] + alt

                    locus = (chrnum, bp, ref, alt)

                    if not locus in vcfs:
                        vcfs[locus] = [id, ".", ".", info]

# Output
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

for locus in sorted(vcfs, key=operator.itemgetter(0, 1)):
    chrchar = ""

    if locus[0] == 97 or locus[0] == 98 or locus[0] == 99:
        chrchar = number_chars[locus[0]]
    
    else:
        chrchar = str(locus[0])
    
    pos = locus[1]
    id = vcfs[locus][0]
    ref = locus[2]
    alt = locus[3]
    qual = vcfs[locus][1]
    filter = vcfs[locus][2]
    info = vcfs[locus][3]

    print(chrchar, pos, id, ref, alt, qual, filter, info, sep='\t')
