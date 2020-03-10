#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020,  Magic Fang, magicfang@gmail.com
#
# Distributed under terms of the GPL-3 license.

import sys;
from Bio import SeqIO;
from seqman import argparser, getfh, close;
from seqman.consts import *;
from seqman.commands import *;
from seqman.seqfilter import seqfilter;

parser = argparser(True);
args = parser.parse_args();
cmd = commands(args.command, args);
fil = seqfilter(args);
if not cmd:
	parser.error("command not found!");

# ifh: Input file handle, ofh: Output file handle
(ifh, ofh) = (args.infile, args.output);
ifh = getfh(args.infile, "rt");
ofh = getfh(args.output, "wt");
if not ifh or not ofh:
	parser.error("input or output file error");

#Initial output file format, default is input file format
if not args.oformat:
	args.oformat = args.format;

for seq in SeqIO.parse(ifh, args.format, alphabet = ALPHABET[args.alphabet]):
	if args.upper or args.lower:
		seq = fil.upperlower(seq);
	if not fil.check(seq):
		continue;
	cseq = cmd.run(seq);
	if cseq:
		if args.command == "split":
			SeqIO.write(seq, cseq, args.oformat);
		elif args.command == "stat" or args.command == "translate" and args.six:
			ofh.write(cseq);
		else:
			SeqIO.write(cseq, ofh, args.oformat);

# Input and output file handle destroy
close(ifh);
close(ofh);
