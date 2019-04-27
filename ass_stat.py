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
parser.add_argument("infile", help = "input file", nargs="+");
args = parser.parse_args();

stat={};
ktup=();
for fn in args.infile:
	mat = re.match(r'GCA\_\d+\.\d+\_([a-z|A-Z|0-9]+)[\.\_]', fn);
	if mat:
		ass = mat.group(1); 
		stat[ass]={};
		fh = open(fn, "r");
		for ln in fh:
			mat = re.match(r'all\tall\tall\tall\t([\w\-]+)\t(\d+)', ln);
			if mat:
				stat[ass][mat.group(1)]=mat.group(2);
				if not mat.group(1) in ktup:
					ktup+=tuple([mat.group(1)]);

l=" \t";
for n in ktup:
	l+=n+"\t";
print l;
l="";
for k in stat:
	l=k+"\t";
	for n in ktup:
		if n in stat[k]:
			l+=str(stat[k][n])+"\t";
		else:
			l+=" \t";
	print l;
	l="";
				
