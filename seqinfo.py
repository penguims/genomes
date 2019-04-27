#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © %YEAR% %USER% <%MAIL%>
#
# Distributed under terms of the %LICENSE% license.

import sys;
import re;
import argparse;
from Bio import SeqIO;

"""
%HERE%
"""


parser = argparse.ArgumentParser();
parser.add_argument("-i", "--infile", help = "input file", required = "yes");
parser.add_argument("-l", "--len", help = "", type = int, default = 4000000);
parser.add_argument("-o", "--outfile", help = "output file");
args = parser.parse_args();

seqs = [];

for seq in SeqIO.parse(args.infile, "fasta"):
	if (len(seq)>=args.len):
		SeqIO.write(seq, seq.id+".fa", "fasta");
	else:
		seqs.append(seq);

SeqIO.write(seqs, "hg19_nochrs.fa", "fasta");

