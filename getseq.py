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
parser.add_argument("-b", "--begin", help = "seq start site", type = int, required = "yes");
parser.add_argument("-e", "--end", help = "seq end site", type = int, required = "yes");
parser.add_argument("-a", "--app", help = "appendix length", type = int, default = 10);
args = parser.parse_args();

for seq in SeqIO.parse(args.infile, "fasta"):
	nseq=seq.seq;
	print("{}|{}|{}".format(
		nseq[args.begin-args.app:args.begin], 
		nseq[args.begin:args.end], 
		nseq[args.end+1:args.end+1+args.app])
	);

