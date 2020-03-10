#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020,  Magic Fang, magicfang@gmail.com
#
# Distributed under terms of the GPL-3 license.

import sys;
sys.path.append("/Users/magic/Downloads/genomes-master");
from Bio import SeqIO;
from seqman import argparser, getfh, close;
from seqman.consts import *;
from seqman.commands import *;
from seqman.seqfilter import seqfilter;

parser = argparser(True);
args = parser.parse_args();
cmd = commands(args.command, args);
if not cmd:
	parser.print_help();
	parser.exit("command not found!");
fil = seqfilter(args);
if args.show:
	showitem(args);
	exit(0);

#Initial output file format, default is input file format
if not args.oformat:
	args.oformat = args.format;

# ifh: Input file handle, ofh: Output file handle
(ifh, ofh) = (args.infile, args.output);
ifh = getfh(args.infile, "rt");
ofh = getfh(args.output, "wt");
if not ifh or not ofh:
	parser.parse_args([args.command, "-h"]);
#	parser.exit("no input file!");


ocseq = "";
for seq in SeqIO.parse(ifh, args.format, alphabet = ALPHABET[args.alphabet]):
	if args.upper or args.lower:
		seq = fil.upperlower(seq);
	if not fil.check(seq):
		continue;
	cseq = cmd.run(seq);
	if cseq:
		if args.command == "split":
			if not ocseq or ocseq != cseq:
				close(ofh);
				ocseq = cseq;
				ofh = getfh(ocseq, "wt");
			SeqIO.write(seq, ofh, args.oformat);
		elif args.command == "stat" or args.command == "translate" and args.six:
			ofh.write(cseq);
		else:
			SeqIO.write(cseq, ofh, args.oformat);

# Input and output file handle destroy
close(ifh);
close(ofh);
