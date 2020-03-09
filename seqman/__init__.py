#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020,  Magic Fang, magicfang@gmail.com
#
# Distributed under terms of the GPL-3 license.


import argparse, sys, gzip, bz2;
from random import randint;
from seqman.consts import *;

__all__ = ["commands", "consts", "seqfilter", "argparser", "getfh", "close", "genseq"];

def genseq(size):
	seq = "";
	for i in range(size):
		seq += NUCLS[randint(0, len(NUCLS)-1)];
	return seq;

def getfh(fn, mode):
	if fn == "stdout":
		return sys.stdout;
	elif fn == "stdin":
		if sys.stdin.isatty():
			return None;
		else:
			return sys.stdin;
	elif fn == "stderr":
		return sys.stderr;
	else:
		try:
			if fn[-3:] == ".gz":
				return gzip.open(fn, mode);
			elif fn[-4:] == ".bz2":
				return bz2.open(fn, mode);
			else:
				return open(fn, mode);
		except:
			return None;

def close(fh):
    if fh == sys.stderr or fh == sys.stdin or sys.stdout:
        return True;
    else:
        try:
            fh.close();
            return True;
        except:
            return False;

def argparser(parser = False):
	pparser = argparse.ArgumentParser(add_help=False);
	pparser.add_argument("-in", "--infile", help = "input file", default = "stdin");
	pparser.add_argument("-ou", "--output", help = "output file name", default = "stdout");
	pparser.add_argument("-if", "--format", help = "sequence file format", default = "fasta");
	pparser.add_argument("-of", "--oformat", help = "output sequence format");
	pparser.add_argument("-fi", "--filter", help = "sequence filter", nargs = "*");
	pparser.add_argument("-vo", "--verbose", action = "store_true", help = "verbose output to stderr");
	pparser.add_argument("-al", "--alphabet", choices = ALPHABET.keys(), help = "sequence alphabet", default = "dna");
	pparser.add_argument("-lo", "--lower", action = "store_true", help = "get lower case sequence");
	pparser.add_argument("-up", "--upper", action = "store_true", help = "get upper case sequence");
	
	cparser = argparse.ArgumentParser();
	subpar = cparser.add_subparsers(dest = "command");
	
	stat_parser = subpar.add_parser("stat", help = "list sequences stat information", parents=[pparser]);
	stat_parser.add_argument("-ds", "--desc", help = "no description when stat", action = "store_true");
	stat_parser.add_argument("-mw", "--weight", help = "molecular weight", action = "store_true");
	stat_parser.add_argument("-3g", "--gc123", help = "3 position gc content", action = "store_true");
	stat_parser.add_argument("-ck", "--checksum", help = "seqeucne checksum", choices=CHECKSUM.keys(), default = "seguid");
	stat_parser.add_argument("-ip", "--isopoint", help = "isoelectric point", action = "store_true");
	stat_parser.add_argument("-mp", "--meltpoint", help = "melting point", action = "store_true")
	stat_parser.add_argument("-ic", "--icc", help = "local component complecity", action = "store_true");
	
	get_parser = subpar.add_parser("get", help = "get a sequence by a name or a list", parents=[pparser]);
	get_parser.add_argument("-nm", "--name", help = "sequence to get");
	get_parser.add_argument("-li", "--list", help = "sequence list");
	get_parser.add_argument("-ft", "--features", help = "get sequence by features", nargs= "+");
	
	split_parser = subpar.add_parser("split", help = "split a bundle squence to a single sequence file", parents=[pparser]);
	split_parser.add_argument("-di", "--dir", help = "split sequence into dir", default = ".");
	
	cut_parser = subpar.add_parser("cut", help = "cut a sequence", parents=[pparser]);
	cut_parser.add_argument("-st", "--start", help = "sequence start", type = int, default = 1);
	cut_group = cut_parser.add_mutually_exclusive_group();
	cut_group.add_argument("-en", "--end", help = "sequence end", type = int);
	cut_group.add_argument("-ln", "--length", help = "sequence end", type = int);
	cut_parser.add_argument("-ra", "--random", help = "random cut a sequence by a length", action = "store_true");
	cut_parser.add_argument("-fr", "--frequency", help = "cut times", type = int, default = 1);
	cut_parser.add_argument("-ml", "--minlen", help = "mix cut length", type = int, default = 10);
	cut_parser.add_argument("-xl", "--maxlen", help = "max cut length", type = int, default = 100);
	cut_parser.add_argument("-re", "--reverse", help = "cut and concate sequence", action = "store_true");
	
	tran_parser = subpar.add_parser("translate", help = "translate dna to protein, or in 6 frames", parents = [pparser]);
	tran_parser.add_argument("-f6", "--six", help = "translate in 6 frames", action = "store_true");
	tran_parser.add_argument("-ct", "--codon", help = "codon table for translation", choices = CODON.keys(), default = "Standard");
	tran_parser.add_argument("-lc", "--list", help = "list codon table", action = "store_true");
	tran_parser.add_argument("-st", "--start", help = "start translate", type = int);
	tran_parser.add_argument("-en", "--end", help = "start translate", type = int);
	tran_parser.add_argument("-fr", "--frame", help = "start translate", type = int);


	trab_parser = subpar.add_parser("transcribe", help = "transcribe dna to rna or backward", parents = [pparser]);
	
	rev_parser = subpar.add_parser("transform", help = "reverse and complimant sequence", parents = [pparser]);
	rev_parser.add_argument("-rv", "--reverse", help = "reverse sequence", action = "store_true");
	rev_parser.add_argument("-co", "--complement", help = "get complememt sequence", action = "store_true");
	
	conv_parser = subpar.add_parser("convert", help = "convert from one format to another", parents = [pparser]);
	
	mut_parser = subpar.add_parser("mutation", help = "make mutation on sequences", parents = [pparser]);
	mut_parser.add_argument("-tp", "--type", help = "mutation type", nargs = "+", choices = MUTYPE);
	mut_parser.add_argument("-fr", "--frequency", help = "mutation frequency", default = 1, type = int);
	mut_parser.add_argument("-fc", "--force", help = "force snp to be not same mutation", action = "store_true");
	mut_parser.add_argument("-ln", "--length", help = "common length to INS, DEL, SV and CNV", default = 1, type = int);
	mut_parser.add_argument("-sp", "--sv_span", help = "SV span length", default = 1, type = int);
	mut_parser.add_argument("-ra", "--random", help = "random mutation length", action = "store_true");
	mut_parser.add_argument("-rf", "--ran_freq", help = "random mutation frequency", action = "store_true");
	mut_parser.add_argument("-ir", "--min_random", help = "min mutation length", default = 1, type = int);
	mut_parser.add_argument("-xr", "--max_random", help = "max mutation length", default = 10, type = int);
	mut_parser.add_argument("-si", "--single", help = "output single mutation sequence", action = "store_true");
	mut_parser.add_argument("-sy", "--stype", help = "output single mutation type sequence", action = "store_true");
	mut_parser.add_argument("-sq", "--seq", help = "seq to insert, sv or cnv");
	mut_parser.add_argument("-cf", "--cnv_freq", help = "cnv repeat times", default = 10, type = int);
	mut_parser.add_argument("-st", "--start", help = "start location for ins, del, sv and cnv", type = int);
	mut_parser.add_argument("-en", "--end", help = "end location for del, sv and cnv", type = int);
	mut_parser.add_argument("-rv", "--reverse", help = "reverse sequence", action = "store_true");
	mut_parser.add_argument("-co", "--complement", help = "get complememt sequence", action = "store_true");
	mut_parser.add_argument("-rc", "--rev_com", help = "get reverse complememt sequence", action = "store_true");

	cnm_parser = subpar.add_parser("chgname", help = "change sequence name", parents = [pparser]);
	cnm_parser.add_argument("-or", "--origin", help = "seq id to be changed");
	cnm_parser.add_argument("-to", "--to", help = "seq changed to name");
	cnm_parser.add_argument("-li", "--list", help = "sequnce list: from_id to_id");
	cnm_parser.add_argument("-kp", "--keep", help = "keep origin seq id", action = "store_true");
	
	if parser:
		return cparser;
	else:
		return	cparser.parse_args();

