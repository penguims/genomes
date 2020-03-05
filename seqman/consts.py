#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020,  Magic Fang, magicfang@gmail.com
#
# Distributed under terms of the GPL-3 license.

from Bio.SeqUtils.CheckSum import crc32, crc64, seguid, gcg;
from Bio.Alphabet import generic_dna, generic_rna, generic_protein;

"""
%HERE%
"""
# DNA nucleics
NUCLS = ['A', 'C', 'G', 'T'];

# Sequence alphabet dataset
ALPHABET = {
	"dna":generic_dna,
	"rna":generic_rna,
	"pro":generic_protein,
	};

# NCBI codon tables
CODON = {
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

MUTYPE = ["snp", "del", "ins", "sv", "cnv"];

CHECKSUM = {
		"seguid":seguid,
		"gcg":gcg,
		"crc32":crc32,
		"crc64":crc64
	};

