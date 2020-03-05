#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © %YEAR% %USER% <%MAIL%>
#
# Distributed under terms of the %LICENSE% license.

import sys;
import os;
import argparse;
from random import randint;
from Bio.Seq import Seq;
from Bio.SeqRecord import SeqRecord;
from Bio.SeqFeature import SeqFeature, FeatureLocation;
from Bio import SeqUtils;
from seqman.consts import CHECKSUM, ALPHABET;

CMDS = ("get", "cut", "split", "stat", "chgname");

def commands(cmdstr, args):
	if cmdstr in CMDS:
		return globals()[cmdstr](args);
	else:
		return None;

class chgname:
	def __init__(self, args):
		self.origin = args.origin;
		self.to = args.to;
		self.list = args.list;
		self.keep = args.keep;
		self.ndict = {};
		try:
			fh = open(self.list, "r");
			for ln in fh.readlines():
				ln = ln.strip();
				tmp = ln.split("\t");
				self.ndict[tmp[0]] = tmp[1];
			fh.close();
		except:
			print("Can not open list file: "+self.list, file = sys.stderr);
		if self.origin and self.to:
			self.ndict[self.origin] = self.to;

	def run(self, seq):
		if seq.id in self.ndict.keys():
			tmp = seq.id;
			seq.id = self.ndict[seq.id];
			if self.keep:
				seq.description = " ".join([tmp, seq.description]);
			self.ndict.pop(tmp);
			return seq;
		else:
			return None;

	def done(self):
		return len(self.ndict) == 0;


class cut:
	def __init__(self, args):
		self.start = args.start;
		self.end = args.end;
		self.length = args.length;
		self.random = args.random;
		self.freq = args.frequency;
		self.minlen = args.minlen;
		self.maxlen = args.maxlen;
	
	def cutseq(self, seq):
		(start, end) = (self.start, self.start+self.length);
		if self.end:
			end = self.end;
		if end > len(seq):
			end = len(seq);
		seqid = "{}_{}_{}".format(seq.id, start, end);
		return(SeqRecord(seq.seq[start:end], id = seqid, description = seq.description));
	
	def randomcut(self, seq):
		subseqs = [];
		for i in range(self.freq):
			(start, end) = (0, 0);
			if self.length:
				start = randint(0, len(seq)-self.length);
				end = start + self.length;
			else:
				ln = randint(self.minlen, self.maxlen);
				start = randint(0, len(seq)-ln);
				end = start + ln;
			seqid = "{}_{}_{}".format(seq.id, start, end);
			subseqs.append(SeqRecord(seq.seq[start:end], id = seqid, description = seq.description));
		return subseqs;
	
	def run(self, seq):
		if self.random:
			return self.randomcut(seq);
		else:
			return self.cutseq(seq);

class get:
	def __init__(self, args):
		self.name = args.name;
		self.nlist = [];
		if args.list:
			try:
				fh = open(args.list, "r");
				for ln in fh.readlines():
					ln = ln.strip();
					self.nlist.append(ln);
				fh.close();
			except:
				print("File not found!", file = sys.stderr);
				exit(1);
		if self.name:
			self.nlist.append(self.name);
		self.feat = args.features;
		self.alphabet = args.alphabet;
	
	def getseq(self, seq):
		if seq.id in self.nlist:
			self.nlist.remove(seq.id);
			return seq;
		else:
			return None;

	def getfeature(self, seq):
		subseqs = [];
		for ft in seq.features:
			if not ft.type in self.feat:
				continue;
			loc = ft.location;
			sid = "{}_{}_{}".format(seq.id, loc.start, loc.end);
			sds = "{} {} {}".format(ft.type, ft.strand, seq.description);
			fseq = Seq(str(ft.extract(seq).seq), alphabet = ALPHABET[self.alphabet]);
			subseqs.append(SeqRecord(fseq, id = sid, description = sds));
		return subseqs;

	def run(self, seq):
		if self.feat:
			return self.getFeatureseq(seq);
		else:
			return self.getseq(seq);

	def done(self):
		if self.nlist:
			return False;
		else:
			return True;
class split:
	def __init__(self, args):
		self.dir = args.dir;
		self.oformat = args.oformat;
		if not os.path.exists(self.dir):
			os.mkdir(self.dir);
	
	def run(self, seq):
		return os.path.join(self.dir, ".".join([seq.id, self.oformat]));

class stat:
	def __init__(self, args):
		self.gc123 = args.gc123;
		self.weight = args.weight;
		self.checksum = args.checksum;
		self.desc = args.desc;
		self.isopoint = args.isopoint;
		self.meltpoint = args.meltpoint;
		self.icc = args.icc;
	
	def run(self, seq):
		items = [seq.id, str(len(seq))];
		if self.gc123:
			for gc in SeqUtils.GC123(seq.seq):
				items.append("{:2.2f}".format(gc));
		else:
			items.append("{:2.2f}".format(SeqUtils.GC(seq.seq)));
		if self.weight:
			items.append("{:.2f}".format(SeqUtils.molecular_weight(seq.seq)));
		if self.checksum:
			cksum = CHECKSUM[self.checksum](seq.seq);
			if type(cksum) is int:
				items.append(str(cksum));
			else:
				items.append(cksum);
		if self.desc:
			items.append(seq.description);
		return "\t".join(items)+"\n";
