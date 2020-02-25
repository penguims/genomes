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

parser = argparse.ArgumentParser();

subpar = parser.add_subparsers(dest = "command");

list_parser = subpar.add_parser("list", help = "list sequences");
list_parser.add_argument("-i", "--infile", help = "input file", required = "yes");
list_parser.add_argument("-f", "--format", help = "sequence file format", default = "fasta");

stat_parser = subpar.add_parser("stat", help = "list sequences stat information");
stat_parser.add_argument("-s", "--nodes", help = "no description when stat", action = "store_true");
stat_parser.add_argument("-i", "--infile", help = "input file", required = "yes");
stat_parser.add_argument("-f", "--format", help = "sequence file format", default = "fasta");

get_parser = subpar.add_parser("get", help = "get a sequence by name");
get_parser.add_argument("-i", "--infile", help = "input file", required = "yes");
get_parser.add_argument("-f", "--format", help = "sequence file format", default = "fasta");
get_parser.add_argument("-n", "--name", help = "sequence to get", required = "yes");
get_parser.add_argument("-o", "--output", help = "output file name", required = "yes");

extr_parser = subpar.add_parser("extract", help = "extract a list of seqeunce");
extr_parser.add_argument("-i", "--infile", help = "input file", required = "yes");
extr_parser.add_argument("-f", "--format", help = "sequence file format", default = "fasta");
extr_parser.add_argument("-l", "--list", help = "sequence list", required = "yes");
extr_parser.add_argument("-o", "--output", help = "output file name", required = "yes");

split_parser = subpar.add_parser("split", help = "split a bundle squence to a single sequence file");
split_parser.add_argument("-i", "--infile", help = "input file", required = "yes");
split_parser.add_argument("-f", "--format", help = "sequence file format", default = "fasta");
split_parser.add_argument("-d", "--dir", help = "split sequence file into the dir", default = ".");

cut_parser = subpar.add_parser("cut", help = "cut a sequence");
cut_parser.add_argument("-i", "--infile", help = "input file", required = "yes");
cut_parser.add_argument("-f", "--format", help = "sequence file format", default = "fasta");
cut_parser.add_argument("-b", "--begin", help = "sequence start", type = int);
cut_parser.add_argument("-e", "--end", help = "sequence end", type = int);

args = parser.parse_args();

total = 0;

if args.command == "split" and not os.path.exists(args.dir):
	os.mkdir(args.dir);

fh = None;
nlist = set();
if args.command == "extract":
	if args.list:
		fh = open(args.list, mode = "r");
		for ln in fh.readlines():
			ln = ln.strip();
			nlist.add(ln);
		fh.close();
	else:
		parser.error("Extract command need a sequence id list file!");
	fh = open(args.output, mode = "w");

for seq in SeqIO.parse(args.infile, args.format):
	total += 1;
	if args.command == "list":
		print(seq.id);
	elif args.command == "stat":
		if args.nodes:
			print("{0}\t{1}\t{2:2.2f}%".format(seq.id, len(seq), SeqUtils.GC(seq.seq)));
		else:
			print("{0}\t{1}\t{2:2.2f}%\t{3}".format(seq.id, len(seq), SeqUtils.GC(seq.seq), seq.description));
	elif args.command == "get":
		if seq.id == args.name:
			SeqIO.write(seq, args.output, args.format);
			break;
	elif args.command == "split":
		ofn = seq.id+"."+args.format;
		ofn = os.path.join(args.dir, ofn);
		SeqIO.write(seq, ofn, args.format);
	elif args.command == "cut":
		if seq.id == args.name:
			subid = seq.id+"_"+str(args.begin)+"_"+str(args.end);
			ofn = subid+"."+args.format;
			nseq = SeqRecord(Seq(seq[args.begin:args.end], IUPAC.IUPACAmbiguousDNA), id = subid);
			SeqIO.write(nseq, ofn, args.format);
			break;
	elif args.command == "extract" and seq.id in nlist:
		SeqIO.write(seq, fh, args.format);	
		nlist.remove(seq.id);
		if not len(nlist):
			break;

if args.command == "extract":
	fh.close();

if args.command == "list":
	print("=====\nTotal: {0}".format(total));
