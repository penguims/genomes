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
import Bio;

"""
%HERE%
"""


parser = argparse.ArgumentParser();
parser.add_argument("-i", "--infile", help = "input file", required = "yes");
parser.add_argument("-g", "--gaps", help = "gap list", required = "yes");
parser.add_argument("-w", "--window", type = int, help = "overlap range", default = 10);
parser.add_argument("-o", "--outfile", help = "output file");
args = parser.parse_args();

rg = [];
ct = 0;
with open(args.gaps, "r") as fh:
	for ln in fh:
		m = re.match("^\((\d+)\,\s+(\d+)\)", ln);
		rg.append([]);
		rg[ct].append(int(m.group(1)));
		rg[ct].append(int(m.group(2)));
		ct += 1;

with open(args.infile, "r") as fh:
	for ln in fh:
		cls = ln.split("\t");
		s = int(cls[8]);
		e = int(cls[9]);
		for r in rg:
			if e >= r[0]-args.window and e <= r[0] or s >= r[1] and s <= r[1]+args.window:
				print(ln);

