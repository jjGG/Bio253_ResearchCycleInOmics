#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 21:58:28 2022

@author: nataliazajac
#"""
#"""Main function."""
#    
#    usage = """
#    %python extract_from_fasta.py [fasta] [start] [end]
#    Note:
#     list of samples has a single sample per line
# """

import sys
from Bio import SeqIO


fasta = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])

record = SeqIO.read(fasta, "fasta")

with open("my_region.fa", "w") as out:
	dna = record[start:end]
	SeqIO.write(dna, out, "fasta")

with open("my_region.protein.fa", "w") as out:
	protein = dna.translate()
	SeqIO.write(protein, out, "fasta")
