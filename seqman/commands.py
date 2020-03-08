#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020,  Magic Fang, magicfang@gmail.com
#
# Distributed under terms of the GPL-3 license.

import sys;
import os;
import argparse;
from random import randint;
from Bio.Seq import Seq;
from Bio.SeqRecord import SeqRecord;
from Bio.SeqFeature import SeqFeature, FeatureLocation;
from Bio import SeqUtils;
from seqman.consts import CHECKSUM, ALPHABET, CODON, NUCLS;
from seqman import genseq;

CMDS = ("get", "cut", "split", "stat", "chgname", "translate", "transcribe", "transform", "convert", "mutation");

def commands(cmdstr, args):
	if cmdstr in CMDS:
		return globals()[cmdstr](args);
	else:
		return None;

class mutation:
	def __init__(self, args):
		self.type = args.type;
		self.freq = args.frequency;
		self.force = args.force;
		self.length = args.length;
		self.random = args.random;
		self.sv_span = args.sv_span;
		self.ran_freq = args.ran_freq;
		self.min_random = args.min_random;
		self.max_random = args.max_random;
		self.single = args.single;
		self.verbose = args.verbose;
		self.stype = args.stype;
		self.seq = args.seq;
		self.cnv_freq = args.cnv_freq;
		self.start = args.start;
		self.end = args.end;
		self.reverse = args.reverse;
		self.complement = args.complement;
		self.rev_com = args.rev_com;
	def _insert(self, museq, p, iseq):
		if p == 0:
			return iseq+museq;
		else:
			return museq[0:p]+iseq+museq[p:];
	def _delete(self, museq, p, rlen):
		if p == 0:
			return museq[rlen:];
		else:
			return museq[0:p]+museq[p+rlen:];
	def snp(self, seq):
		museq = seq.seq.tomutable();
		p = randint(0, len(museq)-1);
		if getattr(self, "start", True):
			p = self.start;
		ns = NUCLS[0:];
		if self.force:
			ns.remove(museq[p].upper());
		tmp = museq[p];
		museq[p]=ns[randint(0, len(ns)-1)];
		vstr = "{}:{}->{}".format(p, tmp, museq[p]);
		if self.verbose:
			print(vstr, file = sys.stderr);
		return SeqRecord(museq.toseq(), id = seq.id, description = vstr+" "+seq.description);
	def ins(self, seq):
		(rlen, museq) = (self.length, seq.seq.tomutable());
		p = randint(0, len(seq)-1);
		if getattr(self, "start", True):
			p = self.start;
		tmp = self.seq;
		if not getattr(self, "seq", True):
			tmp = genseq(rlen);
		museq = self._insert(museq, p, tmp);
		vstr = ("{}:{}").format(tmp, p);
		if self.verbose:
			print(seq.seq.tomutable(), file = sys.stderr);
			print(("{:>"+str(p+rlen)+"}").format(vstr), file = sys.stderr);
			print(museq, file = sys.stderr);
		return SeqRecord(museq.toseq(), id = seq.id, description = vstr+" "+seq.description);
	def delete(self, seq):
		(p, rlen, museq) = (0, self.length, seq.seq.tomutable());
		if getattr(self, "start", True):
			p = self.start;
			if getattr(self, "end", True):
				rlen = self.end - self.start + 1;
		else:
			p = randint(0, len(museq)-1-rlen);
		if len(museq) <= rlen:
			print("Sequence {} length deleted to be zero!".format(seq.id), file = sys.stderr);
			exit(1);
		tmp = str(museq[p:p+rlen]);
		museq = self._delete(museq, p, rlen);
		vstr = ("{}:{}").format(tmp, p)
		if self.verbose:
			print(seq.seq.tomutable(), file = sys.stderr);
			print(("{:>"+str(p+rlen)+"}").format(vstr), file = sys.stderr);
			print(museq, file = sys.stderr);
		return SeqRecord(museq.toseq(), id = seq.id, description = vstr+" "+seq.description);
	def sv(self, seq):
		(p, rlen, museq) = (0, self.length, seq.seq.tomutable());
		if getattr(self, "start", True):
			p = self.start;
			if getattr(self, "end", True):
				rlen = self.end - self.start + 1;
		else:
			p = randint(0, len(museq)-1-rlen);
		tmp = museq[p:p+rlen];
		otmp = tmp[0:];
		if self.reverse:
			tmp.reverse();
		elif self.complement:
			tmp.complement();
		elif self.rev_com:
			tmp.reverse_complement();
		t = p + rlen + self.sv_span;
		if t > len(museq)-1:
			t = len(museq);
		museq = self._insert(museq, t, tmp);
		museq = self._delete(museq, p, rlen);
		vstr = "{}:{}->{}".format(tmp, p, t);
		if otmp != tmp:
			vstr = "{}->{}:{}->{}".format(otmp, tmp, p, t);
		if self.verbose:
			pass;
		return SeqRecord(museq.toseq(), id = seq.id, description = vstr+" "+seq.description);
			
	def cnv(self, seq):
		museq = seq.seq.tomutable();
		fr = self.cnv_freq;
		rlen = self.length;
		if getattr(self, "start", True):
			p = self.start;
			if getattr(self, "end", True):
				rlen = self-end - self.start + 1;
		else:
			p = randint(0, len(seq)-1-rlen);
		tmp = self.seq;
		if not self.seq:
			tmp = museq[p:p+rlen];
			fr -= 1;
		vstr = "{}:{}*{}".format(p, tmp, fr);
		museq = self._insert(museq, p, tmp*fr);
		return SeqRecord(museq.toseq(), id = seq.id, description = vstr+" "+seq.description);
	def run(self, seq):
		seqs = [];
		mseq = seq[0:];
		for t in self.type:
			fr = self.freq;
			if self.ran_freq:
				fr = randint(1, self.ran_freq);
			for f in range(fr):
				if self.random:
					rlen = randint(self.min_random, self.max_random);
				if t == "snp":
					mseq = self.snp(mseq);
				elif t == "ins":
					mseq = self.ins(mseq);
				elif t == "del":
					mseq = self.delete(mseq);
				elif t == "sv":
					mseq = self.sv(mseq);
				elif t == "cnv":
					mseq = self.cnv(mseq);
				if self.single:
					seqs.append(mseq);
					mseq = seq[0:];
			if self.stype and not self.single:
				seqs.append(mseq);
				mseq = seq[0:];
		if not (self.single or self.stype):
			seqs.append(mseq);
		return seqs;

class translate:
	def __init__(self, args):
		self.six = args.six;
		self.codon = args.codon;
		self.start = args.start;
		self.end = args.end;
		self.frame = args.frame;
	def run(self, seq):
		if self.six:
			return "\n".join([
				seq.id, 
				seq.description, 
				SeqUtils.six_frame_translations(seq.seq, CODON[self.codon]).tostring()]);
		else:
			(start, end) = (0, len(seq));
			museq = seq;
			if self.frame:
				if self.frame < 4:
					start = self.frame-1;
					end = len(seq);
				else:
					museq = seq.reverse_complement();
					museq.id = seq.id;
					museq.description = seq.description;
					start = abs(self.frame-6);
					end = len(museq);
				museq.description = "frame: {} {}".format(self.frame, museq.description);
			elif self.start >=self.end:
				start = self.start;
				end = self.end;
				if end > len(seq):
					end = len(seq);
				museq.description = "start: {} end: {} {}".format(start, end, museq.description);
			museq = museq[start:end];
			return museq.translate(table = CODON[self.codon], id = True, description = True);

class trascribe:
	def __init__(self, args):
		self.alphabet = args.alphabet;
	def run(self, seq):
		if self.alphabet == "dna":
			return SeqRecord(seq.seq.transcribe(), id = seq.id, description = seq.description);
		else:
			return SeqRecord(seq.seq.back_transcribe(), id = seq.id, description = seq.description);

class transform:
	def __init__(self, args):
		self.complement = args.complement;
		self.reverse = args.reverse;
	def run(self, seq):
		museq = seq.seq.tomutable();
		if args.complement:
			museq = museq.complement();
		elif args.reverse:
			museq = museq.reverse();
		else:
			museq = museq.reverse_complement();
		return SeqRecord(museq.toseq(), id = seq.id, description = seq.description);

class convert:
	def __init__(self, args):
		pass;
	def run(self, seq):
		return seq;

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
		self.reverse = args.reverse;
	def cutseq(self, seq):
		(start, end) = (self.start, self.start+self.length);
		if self.end:
			end = self.end;
		if end > len(seq):
			end = len(seq);
		seqid = "{}_{}_{}".format(seq.id, start, end);
		mseq = seq.seq[start:end];
		if self.reverse:
			mseq = seq.seq[0:start]+seq.seq[end:];
			seqid = "{}-{}-{}".format(seq.id, start, end);
		return(SeqRecord(mseq, id = seqid, description = seq.description));
	
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
			mseq = seq.seq[start:end];
			if self.reverse:
				mseq = seq.seq[0:start]+seq.seq[end:];
				seqid = "{}-{}-{}".format(seq.id, start, end);
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
