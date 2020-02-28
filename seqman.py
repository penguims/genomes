#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020,  Magic Fang, magicfang@gmail.com
#
# Distributed under terms of the GPL-3 license.

import sys;
import os;
import re;
import argparse;
from Bio import SeqIO;
from Bio import SeqUtils;
from Bio.Seq import Seq;
from Bio.SeqRecord import SeqRecord;
from Bio.SeqFeature import SeqFeature, FeatureLocation;
from Bio.Alphabet import generic_dna, generic_rna, generic_protein, IUPAC;
from random import randint;

"""
%HERE%
"""
# Sequence alphabet dataset
alphabet = {
	"dna":generic_dna,
	"rna":generic_rna,
	"pro":generic_protein,
	};

# NCBI codon tables
codon = {
	"Standard":1,
	"Vertebrate Mitochondrial":2,
	"Yeast Mitochondrial":3,
	"Mold Mitochondrial":4,
	"Protozoan Mitochondrial":4,
	"Coelenterate Mitochondrial":4,
	"Mycoplasma; Spiroplasma":4,
	"Invertebrate Mitochondrial":5,
	"Ciliate Nuclear":6,
	"Dasycladacean Nuclear":6,
	"Hexamita Nuclear":6,
	"Echinoderm Mitochondrial":9,
	"Flatworm Mitochondrial":9,
	"Euplotid Nuclear":10,
	"Bacterial, Archaeal and Plant Plastid":11,
	"Alternative Yeast Nuclear":12,
	"Ascidian Mitochondrial":13,
	"Alternative Flatworm Mitochondrial":14,
	"Blepharisma Macronuclear":15,
	"Chlorophycean Mitochondrial":16,
	"Trematode Mitochondrial":21,
	"Scenedesmus obliquus Mitochondrial":22,
	"Thraustochytrium Mitochondrial":23,
	"Pterobranchia Mitochondrial":24,
	"Candidate Division SR1 and Gracilibacteria":25,
	"Pachysolen tannophilus Nuclear":26,
	"Karyorelict Nuclear":27,
	"Condylostoma Nuclear":28,
	"Mesodinium Nuclear":29,
	"Peritrich Nuclear":30,
	"Blastocrithidia Nuclear":31
	};

# Check filter and return if pass
# seq: Seq Object
# fil: len or gc
# opr: <=, <, =, >, >=
# val: sequence length or gc content
def checkFilter(seq, fil, opr, val):
	tmp = 0;
	if fil == "len":
		tmp = len(seq);
	elif fil == "gc":
		tmp = SeqUtils.GC(seq.seq);
	if opr == ">":
		return tmp > val;
	elif opr == ">=":
		return tmp >= val;
	elif opr == "=":
		return tmp == val;
	elif opr == "<=":
		return tmp <= val;
	elif opr == "<":
		return tmp < val;
	else:
		return False;

# Copy feature informtion from SeqRecodr to Seq
# seq: Seq object
# seqrec: SeqRecord
# return SeqRecodr object
def copySeqRecord(seq, seqrec):
	nseqrec = SeqRecord(seq, id = seqrec.id, description = seqrec.description);
	return nseqrec;

pparser = argparse.ArgumentParser(add_help=False);
pparser.add_argument("-in", "--infile", help = "input file", required = "yes");
pparser.add_argument("-if", "--format", help = "sequence file format", default = "fasta");
pparser.add_argument("-le", "--len", help = "length filter");
pparser.add_argument("-gc", "--gc", help = "gc filter");
pparser.add_argument("-of", "--oformat", help = "output sequence format");
pparser.add_argument("-vo", "--verbose", action = "store_true", help = "verbose output to stderr");
pparser.add_argument("-al", "--alphabet", choices = alphabet.keys(), help = "sequence alphabet", default = "dna");
pparser.add_argument("-lo", "--lower", action = "store_true", help = "get lower case sequence");
pparser.add_argument("-up", "--upper", action = "store_true", help = "get upper case sequence");

cparser = argparse.ArgumentParser();
subpar = cparser.add_subparsers(dest = "command");

list_parser = subpar.add_parser("list", help = "list sequences", aliases=["li"], parents=[pparser]);

stat_parser = subpar.add_parser("stat", help = "list sequences stat information", parents=[pparser]);
stat_parser.add_argument("-no", "--nodes", help = "no description when stat", action = "store_true");

get_parser = subpar.add_parser("get", help = "get a sequence by name", parents=[pparser]);
get_parser.add_argument("-nm", "--name", help = "sequence to get", required = "yes");
get_parser.add_argument("-ou", "--output", help = "output file name", required = "yes");

extr_parser = subpar.add_parser("extract", help = "extract a list of seqeunce", parents=[pparser]);
extr_parser.add_argument("-li", "--list", help = "sequence list", required = "yes");
extr_parser.add_argument("-ou", "--output", help = "output file name", required = "yes");

split_parser = subpar.add_parser("split", help = "split a bundle squence to a single sequence file", parents=[pparser]);
split_parser.add_argument("-di", "--dir", help = "split sequence file into the dir", default = ".");

cut_parser = subpar.add_parser("cut", help = "cut a sequence", parents=[pparser]);
cut_parser.add_argument("-st", "--start", help = "sequence start", required = "yes", type = int, default = 1);
cut_group = cut_parser.add_mutually_exclusive_group();
cut_group.add_argument("-en", "--end", help = "sequence end", type = int);
cut_group.add_argument("-ln", "--length", help = "sequence end", type = int);
cut_parser.add_argument("-ra", "--random", help = "random cut a sequence by a length", action = "store_true");
cut_parser.add_argument("-ou", "--output", help = "output file name", required="yes");
cut_parser.add_argument("-fr", "--frequency", help = "cut times", type = int, default = 1);
cut_parser.add_argument("-ml", "--minlen", help = "mix cut length", type = int, default = 10);
cut_parser.add_argument("-xl", "--maxlen", help = "max cut length", type = int, default = 100);

tran_parser = subpar.add_parser("translate", help = "translate dna to protein, or in 6 frames", parents = [pparser]);
tran_parser.add_argument("-f6", "--six", help = "translate in 6 frames", action = "store_true");
tran_parser.add_argument("-ct", "--codon", help = "codon table for translation", choices = codon.keys(), default = "Standard");
tran_parser.add_argument("-lc", "--list", help = "list codon table", action = "store_true");
tran_parser.add_argument("-ou", "--output", help = "output file name", required="yes");

tran_parser = subpar.add_parser("transcribe", help = "transcribe dna to rna or backward", parents = [pparser]);
tran_parser.add_argument("-ou", "--output", help = "output file name", required="yes");

rev_parser = subpar.add_parser("transform", help = "reverse and complimant sequence", parents = [pparser]);
rev_parser.add_argument("-ou", "--output", help = "output file name", required="yes");
rev_parser.add_argument("-rv", "--reverse", help = "reverse sequence", action = "store_true");
rev_parser.add_argument("-co", "--complement", help = "get complememt sequence", action = "store_true");

conv_parser = subpar.add_parser("convert", help = "convert from one format to another", parents = [pparser]);
conv_parser.add_argument("-ou", "--output", help = "output file name", required = "yes");

args = cparser.parse_args();


# Defined the global vars
# Sequence number counter
total = 0;
# Command = stat, the sequnce id list
nlist = set();
# fil: {len, gc}; opr:{<, <=, =, >=, >}; val: float
(fil, opr, val) = ("", "", 0);
# ifh: Input file handle, ofh: Output file handle
(ifh, ofh) = (None, None);

# List the NCBI codon tables
if args.command == "translate" and args.list:
	for sp in codon.keys():
		print("{}: {}".format(codon[sp], sp));
	exit();

# Check if the directory exists, mkdir if not
if args.command == "split" and not os.path.exists(args.dir):
	os.mkdir(args.dir);

# Get sequence id list
# File format:
# NM_004483.5
# NM_139022.3
# NM_001319155.2
# NM_001319162.2
# NM_000779.4
# NM_001319161.2
# ...
if args.command == "extract":
	if args.list:
		fh = open(args.list, mode = "r");
		for ln in fh.readlines():
			ln = ln.strip();
			nlist.add(ln);
		fh.close();
	else:
		parser.error("Extract command need a sequence id list file!");
# Parse filter options
if args.len or args.gc:
	if args.len:
		valstr = args.len;
		fil = "len";
	else:
		valstr = args.gc;
		fil = "gc";
	mat = re.match(r'^([\<\>\=]+)\s*(\d+)', valstr);
	if mat:
		opr = mat.group(1);
		val = float(mat.group(2));
# Initial output file format, default is input file format
if not args.oformat:
	args.oformat = args.format;

# Initial input and output file handle, stdin and stdout 
try:
	if args.output == "stdout":
		ofh = sys.stdout;
	else:
		ofh = open(args.output, mode = "w");
except AttributeError:
	pass;

if args.infile == "stdin":
	ifh = sys.stdin;
else:
	ifh = open(args.infile, "r");

# Main program
for seq in SeqIO.parse(ifh, args.format, alphabet = alphabet[args.alphabet]):
	total += 1;
	# Setting sequence alphabet in lower or uper case, or no changed
	if args.lower:
		seq = seq.lower();
	elif args.upper:
		seq = seq.upper();
	# Show verbose debugging info.
	if args.verbose:
		print("current process sequence: ", seq.id, file = sys.stderr);
	# Check if the sequences pass the filter
	if (args.len or args.gc) and not checkFilter(seq, fil, opr, val):
		continue;
	# Just list the sequences id
	if args.command in {"list", "li"}:
		print(seq.id);
	# Get a list of the sequences stat info. with description or not
	elif args.command == "stat":
		if args.nodes:
			print("{0}\t{1}\t{2:2.2f}%".format(seq.id, len(seq), SeqUtils.GC(seq.seq)));
		else:
			print("{0}\t{1}\t{2:2.2f}%\t{3}".format(seq.id, len(seq), SeqUtils.GC(seq.seq), seq.description));
	# Get a sequence by id
	elif args.command == "get":
		if seq.id == args.name:
			SeqIO.write(seq, ofh, args.oformat);
			break;
	# Split a multi-seqeunces file into single sequence file in to the given directory
	elif args.command == "split":
		ofn = seq.id+"."+args.oformat;
		ofn = os.path.join(args.dir, ofn);
		SeqIO.write(seq, ofn, args.oformat);
	# Cut a sequence by location,
	# If random option is presented, program will check if the length var is presented, 
	# will generate the number of random seqnence with fixed lenth, otherwise will
	# generate the sequence by random start position and length.
	# If random option is not presented, it will generate a sequence with fixed length
	# or end position.
	elif args.command == "cut":
		if args.random:
			for i in range(args.frequency):
				(start, end) = (0, 0);
				if args.length:
					start = randint(0, len(seq)-args.length);
					end = start + args.length;
				else:
					ln = randint(args.minlen, args.maxlen);
					start = randint(0, len(seq)-ln);
					end = start+ln;
				seqid = "{}_{}_{}".format(seq.id, start, end);
				SeqIO.write(SeqRecord(seq.seq[start:end], id = seqid, description = seq.description), ofh, args.oformat);
		else:
			start = args.start;
			end = 0;
			if args.end:
				end = args.end;
			else:
				end = start+args.length;
			if end > len(seq):
				end = len(seq);
			seqid = "{}_{}_{}".format(seq.id, start, end);
			SeqIO.write(SeqRecord(seq.seq[start:end], id = seqid, description = seq.description), ofh, args.oformat);
	# Extract sequence(s) by list
	elif args.command == "extract" and seq.id in nlist:
		SeqIO.write(seq, ofh, args.oformat);	
		nlist.remove(seq.id);
		if len(nlist) == 0:
			break;
	# Translate the sequence to protein in six frames or just translate
	# Codon table: NCBI Codon Table Codes
	elif args.command == "translate":
		if args.six:
			ofh.write(SeqUtils.six_frame_translations(seq.seq, codon[args.codon]).tostring());
		else:
			SeqIO.write(seq.translate(table = codon[args.codon], id = True, description = True), ofh, args.oformat);
	# Transcribe DNA to RNA or RNA to DNA based on sequence alphabet
	elif args.command == "transcribe":
		if args.alphabet == "dna":
			SeqIO.write(copySeqRecord(seq.seq.transcribe(), seq), ofh, args.oformat);
		elif args.alphabet == "rna":
			SeqIO.write(copySeqRecord(seq.seq.back_transcribe(), seq), ofh, args.oformat);
		else:
			print("Can not transcribe sequence {}".format(seq.id), file = sys.stderr);
	# Transform sequence: reverse, complement or reverse and complement
	elif args.command == "transform":
		if args.complement:
			SeqIO.write(copySeqRecord(seq.seq.complement(), seq), ofh, args.oformat);
		elif args.reverse:
			tseq = seq.reverse_complement().seq.complement();
			SeqIO.write(copySeqRecord(tseq, seq), ofh, args.oformat);
		else:
			SeqIO.write(seq.reverse_complement(id = True, description = True), ofh, args.oformat);
	# Convert sequence file format
	elif args.command == "convert":
		SeqIO.write(seq, ofh, args.oformat);

# Input and output file handle destroy
if args.infile != "stdin":
	ifh.close();

try:
	if args.output and not args.output == "stdout":
		ofh.close();
except AttributeError:
	pass;

if args.command == "list":
	print("=====\nTotal: {0}".format(total));
