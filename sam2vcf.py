#!/usr/bin/env python

# sam2vcf.py ver. 1.0.7
# Copyright (C) 2019 Yuki Kato
# This script is used to convert SAM alignments for Bivartect into predicted VCF variants.
# The output will be shown in standard out.
# The variation class (VC) consists of:
# SNV: single nucleotide variant
# DIV: deletion/insertion variant
# SV: structural variant of length >= 50 bp
# BP: unassigned breakpoint
# Note: small indels and structural variants on the reverse strand will appear on the forward strand in the VCF output. The leftmost position (POS) with "IMPRECISE" INFO is a best approximate due to the Bivartect's prediction on the reverse strand.
# Usage: Type "./sam2vcf.py -h" in your terminal.

import argparse

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    usage='./sam2vcf.py [option]* <sam>',
    description='sam2vcf.py ver. 1.0.7\n\nThis script is used to convert SAM alignments for Bivartect into predicted VCF variants.\nThe output will be shown in standard out.\n\nThe variation class (VC) consists of:\n\n SNV: single nucleotide variant\n DIV: deletion/insertion variant\n SV: structural variant of length >= 50 bp\n BP: unassigned breakpoint\n\nNote: small indels and structural variants on the reverse strand will appear on the forward strand in the VCF output.\nThe leftmost position (POS) with "IMPRECISE" INFO is a best approximate due to the Bivartect\'s prediction on the reverse strand.'
)
parser.add_argument('sam', metavar='sam <STR>', type=str,
                    help='path to the SAM alignments for Bivartect')
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
vcfs = {} # Dictonary of {(CHROM, POS, REF, ALT):[ID, QUAL, FILTER, INFO]}
char_numbers = {"X":97, "Y":98, "MT":99} # Dictionary of {chars:int}
number_chars = {97:"X", 98:"Y", 99:"MT"} # Dictionary of {int:chars}

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
                
                m3 = r3.match(qname)
                
                if m3:
                    # Parse a Bivartect output header
                    id = int(m3.group(1))
                    vc = m3.group(2)
                    ref = m3.group(3)
                    alt = m3.group(4)

                    # Set appropriate bp, ref, alt & info depending on SVTYPE & strand
                    # Case: SNV
                    if vc == "SNV":
                        info = "VC=SNV" # VC: variation class

                        if flag == 16: # Reverse strand
                            bp = pos - 1
                            ref = revComp(ref)
                            alt = revComp(alt)
                        else: # Forward strand
                            bp = pos + len(seq)
                    
                    # Case: INS
                    elif ref == "-":
                        if "*" in alt:
                            info = "VC=SV"
                        else:
                            info = "VC=DIV" # DIV: deletion/insertion variant
                        
                        if flag == 16: # Reverse strand
                            bp = pos # Best estimate
                            ref = revComp(seq[-1])
                            alt = revComp(seq[-1] + alt)
                            info += ";IMPRECISE"
                        else: # Forward strand
                            bp = pos + len(seq) - 1
                            ref = seq[-1]
                            alt = seq[-1] + alt
                    
                    # Case: DEL
                    elif alt == "-":
                        if "*" in ref:
                            info = "VC=SV"
                        
                        else:
                            info = "VC=DIV"
                        
                        if flag == 16: # Reverse strand
                            bp = pos - len(ref) # Best estimate
                            ref = revComp(seq[-1] + ref)
                            alt = revComp(seq[-1])
                            info += ";IMPRECISE"
                        else: # Forward strand
                            bp = pos + len(seq) - 1
                            ref = seq[-1] + ref
                            alt = seq[-1]

                    # Case: BP
                    else:
                        info = "VC=BP" # BP: unassigned breakpoint

                        if flag == 16: # Reverse strand
                            bp = pos - len(ref) # Best estimate
                            ref = revComp(seq[-1] + ref)
                            alt = revComp(seq[-1] + alt)
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

    print(chrchar, pos, id, ref, alt, qual, filter, sep='\t', end='\t')
    print(info)
