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
args = parser.parse_args();

for seq in SeqIO.parse(args.infile, "fasta"):
	SeqIO.write(seq, seq.id+".fa", "fasta"); 
