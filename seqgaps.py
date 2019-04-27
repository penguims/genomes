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

c = 0;
for seq in SeqIO.parse(args.infile, "fasta"):
	b, l, i = 0, 0, 0;
	print("seq name: " + seq.id);
	for nc in seq.seq:
		if nc == "N" :
			if l == 0:
				b = i;
			l += 1;
		else:
			if l > 0:
				print("(%d, %d) %30s...%d...%30s" % (b, b+l-1, seq.seq[b-20:b+10], l, seq.seq[b+l-1-10:b+l-1+20]));
				l = 0;
				c += 1;
		i += 1;
	if l > 0:
		print("(%d, %d) %30s...%d...%30s" % (b, b+l-1, seq.seq[b-20:b+10], l, seq.seq[b+l-1-10:b+l-1+20]));
		c += 1;
print 'total: ', c;

