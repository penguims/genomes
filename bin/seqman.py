#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020,  Magic Fang, magicfang@gmail.com
#
# Distributed under terms of the GPL-3 license.


from Bio import SeqIO;
from seqman import seqman;
from seqman.consts import ALPHABET;
from seqman.seqfilter import seqfilter;

if __name__ == "__main__":
	sm = seqman("conf");
	sm.check();
	fil = seqfilter(sm.args);
	seqitr = SeqIO.parse(sm.ifh, sm.format, alphabet = ALPHABET[sm.alphabet]);
	if sm.command == "split":
		sm.cmd.run(seqitr, sm.oformat);
	else:
		for seq in seqitr: 
			if sm.args.upper or sm.args.lower:
				seq = fil.upperlower(seq);
			if not fil.check(seq):
				continue;
			cseq = sm.cmd.run(seq);
			if cseq:
				if sm.command == "stat" or sm.command == "translate" and sm.args.six:
					sm.ofh.write(cseq);
				else:
					SeqIO.write(cseq, sm.ofh, sm.oformat);
