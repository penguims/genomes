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
from seqman.commands import *;

__all__ = ["commands", "consts", "seqfilter", "seqman"];


"""
	Seqman class for init program
"""
class seqman:
	'seqman class'
	# Internal method for get file handle
	# Will check the file name is stdout, stdin or stderr
	# Supprt bzip2 and gzip file format
	# fn:		filename, string
	# mode:		open file mode
	# return None if faild get file handle
	def _getfh(self, fn, mode):
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
	# argconf:		module name include a var named "parser"
	# 				the module sould be at the command dir well run program
	def __init__(self, argconf = "args"):
		exec("import "+argconf)
		self.parser = eval(argconf+".parser");
	# Check and validate the command line options
	def check(self):
		# ArgumentParser
		self.args = self.parser.parse_args();
		# Command string
		self.command = self.args.command;
		# seqman.command module
		self.cmd = commands(self.command, self.args);
		if not self.cmd:
			self.parser.print_help();
			self.parser.exit("command not found!");
		self.format = self.args.format;
		if not self.args.oformat:
			self.args.oformat = self.args.format;
		self.oformat = self.format;
		self.alphabet = self.args.alphabet;
		self.infile = self.args.infile;
		self.output = self.args.output;
		self.ifh = self._getfh(self.infile, "rt");
		self.ofh = self._getfh(self.output, "wt");
		# If command with not input and output option
		# Invoked when not stdin input content`
		if not (self.ifh and self.ofh):
			self.parser.parse_args([self.args.command, "-h"]);
	# Destroy class object
	def __del__(self):
		if self.ifh != sys.stdin:
			self.ifh.close();
		if self.ofh != sys.stdout or self.ofh != sys.stderr:
			self.ofh.close();
