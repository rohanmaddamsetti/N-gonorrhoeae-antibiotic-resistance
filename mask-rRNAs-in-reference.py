#!/usr/bin/env python

''' 
mask-rRNAs-in-reference.py by Rohan Maddamsetti.

This script takes in a reference genbank file,
finds all rRNA operons after the first one,
and masks them in the sequence with NNNNNN characters.

This preprocessing should help mutation-calling in the full genome.
''' 

from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq
from Bio.Alphabet import generic_nucleotide
import os

def main():
    
    rRNAs_to_keep = ['NGO_r01', 'NGO_r02', 'NGO_r03']
    rRNAs_to_mask = []
    ref = SeqIO.read("./gonorrhoeae-FA1090.gb.txt","genbank")
    rRNAlist = [feat for feat in ref.features if feat.type == 'rRNA']
    for feat in rRNAlist:
        if feat.qualifiers['locus_tag'][0] not in rRNAs_to_keep:
            rRNAs_to_mask.append(feat)

    mutable_ref = MutableSeq(str(ref.seq))

    for feat in rRNAs_to_mask:
        my_mask = 'N'*len(feat)
        my_start = feat.location.start
        my_end = feat.location.end
        mutable_ref[my_start:my_end] = my_mask

    ## now add the masked sequence to the reference genome object.
    ref.seq = Seq(str(mutable_ref),generic_nucleotide)

    with open("./masked-gonorrhoeae-FA1090.gbk", "w") as masked_handle:
        SeqIO.write([ref], masked_handle, "genbank")

main()
