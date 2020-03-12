#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020,  Magic Fang, magicfang@gmail.com
#
# Distributed under terms of the GPL-3 license.

import re;
from Bio import SeqUtils;

"""
Sequence fileter class
"""

class seqfilter:
	''
	# args.upper:		set sequence to upper case, not effect when genbank file format
	# args.lower:		set sequence to lower case
	# args.filter:		filter option, "len>1000" "gc<50"
	def __init__(self, args):
		self.upper = args.upper;
		self.lower = args.lower;
		self.filter = args.filter;
		self.items = [];
		if self.filter:
			for f in self.filter:
				res = re.match(r'([a-zA-Z]+)([\=\<\>]+)(\d+)$', f, re.I);
				if res:
					self.items.append([res.group(1), res.group(2), float(res.group(3))]);
	# Get upper case sequence
	# seq:			Bio.SeqRecord
	# return Bio.SeqRecord
	def upperseq(self, seq):
		return seq.upper();
	# Get lower case sequence
	# seq:			Bio.SeqRecord
	# return Bio.SeqRecord
	def lowerseq(self, seq):
		return seq.lower();
	# return filter result 
	# val1:			filter name string: len or gc
	# opr:			compare operator, "<", "<=", "=", ">=" or "<"
	# val2:			value
	# return true or false
	def _compare(self, val1, opr, val2):
		if opr == "<":
			return (val1-val2)<0;
		elif opr == "<=":
			return (val1-val2)<=0;
		elif opr == "=":
			return (val1-val2)<0.0000001;
		elif opr == ">=":
			return (val1-val2)>=0;
		elif opr == ">":
			return (val1-val2)>0;
		else:
			return False;
	# Check len filter
	# seq:		Bio.SeqRecord
	# return filter result
	def checklen(self, seq):
		res = True;
		for it in self.items:
			if re.match(r'^len$', it[0], re.I):
				res = self._compare(len(seq), it[1], it[2]);
				if not res:
					break;
		return res;
	# Cehck gc filter
	# seq:		Bio.SeqRecord
	# return filter result
	def checkgc(self, seq):
		res = True;
		for it in self.items:
			if re.match(r'^gc$', it[0], re.I):
				res = self._compare(SeqUtils.GC(seq.seq), it[1], it[2]);
				if not res:
					break;
		return res;
	# Short method to check upper or lower case sequence
	# seq:		Bio.SeqRecord
	# return Bio.SeqRecord
	def upperlower(self, seq):
		if self.upper:
			return self.upperseq(seq);
		elif self.lower:
			return self.lowerseq(seq);
		return None;
	# Short method to check filter result
	# seq:		Bio.SeqRecord
	# return filter result of len and gc filter
	def check(self, seq):
		if len(self.items) > 0:
			return self.checklen(seq) and self.checkgc(seq);
		else:
			return True;
