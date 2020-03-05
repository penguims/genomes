#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020,  Magic Fang, magicfang@gmail.com
#
# Distributed under terms of the GPL-3 license.

import sys;
sys.path.append("/Users/magic/bin/seqman");
from Bio import SeqIO;
from seqman.consts import *;
from seqman.commands import *;
from seqman import argparser;

args = argparser();
cmd = commands(args.command, args);
if not cmd:
	print("command not found!", file = sys.stderr);
	exit(1);

# ifh: Input file handle, ofh: Output file handle
(ifh, ofh) = (args.infile, args.output);
if args.infile != sys.stdin:
	ifh = open(args.infile, "r");
elif sys.stdin.isatty():
	cparser.error("{}: {}".format(args.command, "stdin without input content!"));
if args.output != sys.stdout:
	ofh = open(args.output, mode = "w");

#Initial output file format, default is input file format
if not args.oformat:
	args.oformat = args.format;


# Main program
for seq in SeqIO.parse(ifh, args.format, alphabet = ALPHABET[args.alphabet]):
	cseq = cmd.run(seq);
	if cseq:
		if args.command == "split":
			SeqIO.write(seq, cseq, args.oformat);
		elif args.command == "stat":
			ofh.write(cseq);
		else:
			SeqIO.write(cseq, ofh, args.oformat);

# Input and output file handle destroy
if args.infile != sys.stdin:
	ifh.close();
if args.output != sys.stdout:
	ofh.close();

