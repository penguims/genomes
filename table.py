#!/usr/local/bin/python
# -*- coding: UTF-8 -*-

import sys;
import re;
import argparse;


parser = argparse.ArgumentParser();
parser.add_argument("-i", "--infile", help = "input file", required = "yes");
parser.add_argument("-o", "--outfile", help = "output file");
parser.add_argument("-m", "--comment", help = "comment line");
parser.add_argument("-c", "--columns", help = "columns to left");
parser.add_argument("-r", "--rows", help = "rows to left");
parser.add_argument("-s", "--split", help = "spliter of columns", default = "\t");
parser.add_argument("-e", "--eol", help = "end of line", default = "\n");
args = parser.parse_args();

rows=[];
columns=[];

lc=0;
with open(args.infile, "r") as fh:
	for ln in fh:
		cells = ln.strip(args.eol).split(args.split);
		if args.comment and cells[0][0] == args.comment:
			continue;
		

def get_list(l):
	ret = [];
	t = l.split(",");
	for c in t:
		r = c.split(":");
		if len(r) > 1:
			ret.extend(range(int(r[0]), int(r[1])));
		else:
			ret.append(int(r[0]));
	return ret;
