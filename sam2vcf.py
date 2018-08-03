#!/usr/bin/env python

# sam2vcf.py ver. 1.0.6
# Copyright (C) 2018 Yuki Kato
# This script is used to convert a mapping result SAM file for Bivartect into a sorted list of predicted variants in VCF format.
# The output will be shown in standard out.
# Note: structural variants on the negative strand will appear on the positive strand in the VCF output. The leftmost position (POS) with "IMPRECISE" INFO is a best approximate due to the Bivartect's prediction on the negative strand.
# Usage: Type "./sam2vcf.py -h" in your terminal.

import argparse

parser = argparse.ArgumentParser(
    usage = './sam2vcf.py [option]* <sam>',
    description = 'This script is used to convert a mapping result SAM file for Bivartect into a sorted list of predicted variants in VCF format. The output will be shown in standard out. Note: structural variants on the negative strand will appear on the positive strand in the VCF output. The leftmost position (POS) with "IMPRECISE" INFO is a best approximate due to the Bivartect\'s prediction on the negative strand.',
    epilog = 'sam2vcf.py ver. 1.0.6')
parser.add_argument('sam', metavar = 'sam <STRING>', type = str,
                    help = 'path to the mapping result SAM file for Bivartect')
args = parser.parse_args()

import re
import operator

def revComp(seq):
    if seq == "-":
        return "-"

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

# SAM file
# QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ
r1 = re.compile(r"@")
r2 = re.compile(r"(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\w+)")
r3 = re.compile(r"(\d+):(\w+):([\w\-\*]+):([\w\-\*]+)")
vcf = {} # Dictonary of {(CHROM, POS, REF, ALT):[ID, QUAL, FILTER, INFO]}

with open(args.sam, "r") as f:
    for line in f:
        m1 = r1.match(line)

        if not m1:
            m2 = r2.match(line)

            if m2:
                qname = m2.group(1)
                flag = int(m2.group(2))
                rname = m2.group(3)
                pos = int(m2.group(4))
                seq = m2.group(10)
                bp = 0
                num = qname
                ref = "."
                alt = "."
                info = "."
                
                # Case: unmapped
                if rname == "*":
                    continue
                
                # Case: unmapped
                #if flag == 4 or flag == 77 or flag == 141:
                    #bp = 0

                if rname == "X":
                    rname = 97
                elif rname == "Y":
                    rname = 98
                elif rname == "MT":
                    rname = 99
                
                m3 = r3.match(qname)
                
                if m3:
                    # Parse a Bivartect output header
                    num = m3.group(1)
                    ref = m3.group(3)
                    alt = m3.group(4)
                    info = "SVTYPE=" + m3.group(2)

                    # Set appropriate bp, ref, alt & info depending on SVTYPE & strand
                    # Case: SNV
                    if m3.group(2) == "SNV":
                        if flag == 16: # Negative strand
                            bp = pos - 1
                            ref = revComp(ref)
                            alt = revComp(alt)
                        else: # Positive strand
                            bp = pos + len(seq)
                    
                    # Case: INS
                    elif ref == "-":
                        if flag == 16: # Negative strand
                            bp = pos # Best estimate
                            ref = revComp(seq[-1])
                            alt = revComp(seq[-1] + alt)
                            info += ";IMPRECISE"
                        else: # Positive strand
                            bp = pos + len(seq) - 1
                            ref = seq[-1]
                            alt = seq[-1] + alt
                    
                    # Case: DEL
                    elif alt == "-":
                        if flag == 16: # Negative strand
                            bp = pos - len(ref) # Best estimate
                            ref = revComp(seq[-1] + ref)
                            alt = revComp(seq[-1])
                            info += ";IMPRECISE"
                        else: # Positive strand
                            bp = pos + len(seq) - 1
                            ref = seq[-1] + ref
                            alt = seq[-1]

                    # Case: BP
                    else:
                        if flag == 16: # Negative strand
                            bp = pos - len(ref) # Best estimate
                            ref = revComp(seq[-1] + ref)
                            alt = revComp(seq[-1] + alt)
                            info += ";IMPRECISE"
                        else: # Positive strand
                            bp = pos + len(seq) - 1
                            ref = seq[-1] + ref
                            alt = seq[-1] + alt
                
                tup = (int(rname), bp, ref, alt)
                
                if not tup in vcf:
                    vcf[tup] = [num, ".", ".", info]

# Output
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

for t in sorted(vcf, key=operator.itemgetter(0, 1)):
    if t[0] == 97:
        chrom = "X"
    elif t[0] == 98:
        chrom = "Y"
    elif t[0] == 99:
        chrom = "MT"
    else:
        chrom = t[0]
    
    print(chrom, t[1], vcf[t][0], t[2], t[3], vcf[t][1], vcf[t][2], sep = '\t', end = '\t')
    print(vcf[t][3])
