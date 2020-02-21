#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © %YEAR% %USER% <%MAIL%>
#
# Distributed under terms of the %LICENSE% license.

import sys;
import os;
import re;
import argparse;
from Bio import SeqIO;
from Bio import SeqUtils;
from Bio.Seq import Seq;
from Bio.SeqRecord import SeqRecord;
from Bio.SeqFeature import SeqFeature, FeatureLocation;
from Bio.Alphabet import IUPAC;

"""
%HERE%
"""

comm = ["list", "get", "stat", "split", "cut", "extract"];

parser = argparse.ArgumentParser();
parser.add_argument("-i", "--infile", help = "input file", required = "yes");
parser.add_argument("-c", "--command", help = "functions", default = "list", choices = comm);
parser.add_argument("-f", "--format", help = "sequence file format", default = "fasta");
parser.add_argument("-n", "--name", help = "sequence to get");
parser.add_argument("-l", "--list", help = "sequence list");
parser.add_argument("-b", "--begin", help = "sequence start", type = int);
parser.add_argument("-e", "--end", help = "sequence end", type = int);
parser.add_argument("-o", "--output", help = "output file name");
parser.add_argument("-s", "--nodes", help = "no description when stat", action = "store_true");
args = parser.parse_args();

total = 0;

nlist = set();
seqs = [];

if (args.command in {"split", "extract", "get"} and not args.output):
	parser.error("Extract or get command should have output file!");

if args.command == "split" and not os.path.exists(args.output):
	os.mkdir(args.output);

if (args.command == "extract"):
	fh = open(args.list, mode = "r");
	for ln in fh.readlines():
		ln = ln.strip();
		nlist.add(ln);
	fh.close();

for seq in SeqIO.parse(args.infile, args.format):
	total += 1;
	if (args.command == "list"):
		print(seq.id);
	elif (args.command == "stat"):
		if (args.nodes):
			print("{0}\t{1}\t{2:2.2f}%".format(seq.id, len(seq), SeqUtils.GC(seq.seq)));
		else:
			print("{0}\t{1}\t{2:2.2f}%\t{3}".format(seq.id, len(seq), SeqUtils.GC(seq.seq), seq.description));
	elif (args.command == "get"):
		if (seq.id == args.name):
			seqs.append(seq);
	elif (args.command == "split"):
		ofn = seq.id+"."+args.format;
		if (args.output):
			ofn = args.output+"/"+ofn;
		SeqIO.write(seq, ofn, args.format);
	elif (args.command == "cut"):
		if (seq.id == args.name):
			subid = seq.id+"_"+str(args.begin)+"_"+str(args.end);
			ofn = subid+"."+args.format;
			nseq = SeqRecord(Seq(seq[args.begin:args.end], IUPAC.IUPACAmbiguousDNA), id = subid);
			SeqIO.write(nseq, ofn, args.format);
	elif (args.command == "extract" and seq.id in nlist):
		seqs.append(seq);		

if (seqs):
	SeqIO.write(seqs, args.output, args.format);

if (args.command == "list"):
	print("=====\nTotal: {0}".format(total));
