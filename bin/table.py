#!/usr/local/bin/python
# -*- coding: UTF-8 -*-

import sys;
import re;
import argparse;


def get_list(l):
	ret = [];
	t = l.split(",");
	for c in t:
		r = c.split(":");
		if len(r) > 1:
			ret.extend(range(int(r[0]), int(r[1])+1));
		else:
			ret.append(int(r[0]));
	return ret;

parser = argparse.ArgumentParser();
parser.add_argument("-i", "--infile", help = "input file", required = "yes");
parser.add_argument("-o", "--outfile", help = "output file");
parser.add_argument("-m", "--comment", help = "comment line", action = "store_true");
parser.add_argument("-n", "--symbol", help = "comment start up symbol", default ="#");
parser.add_argument("-c", "--columns", help = "columns to left");
parser.add_argument("-r", "--rows", help = "rows to left");
parser.add_argument("-s", "--split", help = "spliter of columns", default = "\t");
parser.add_argument("-e", "--eol", help = "end of line", default = "\n");
parser.add_argument("-f", "--filter", help = "record filter: cell[column]=='xxx'");
args = parser.parse_args();

rows=[];
cols=[];
if args.rows:
	rows = get_list(args.rows);
if args.columns:
	cols = get_list(args.columns);

with open(args.infile, "r") as fh:
	lc = -1;
	for ln in fh:
		lc += 1;
		cells = ln.strip(args.eol).split(args.split);
		if not args.comment and ln[0] == args.symbol:
			continue;
		if rows and (lc not in rows):
			continue;
		if ln[0] != args.symbol:
			if args.filter and (not eval(args.filter)):
				continue;
		cc=-1;
		for cell in cells:
			cc += 1;
			if cols and (cc not in cols):
				continue;
			print(cell, args.split, end="");
		print();
		

