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

#Support file format
FORMAT = {
		'abi':'Applied Biosystem\'s sequencing trace format',
		'abi-trim':'Same as "abi" but with quality trimming with Mott\'s algorithm',
		'ace':'Reads the contig sequences from an ACE assembly file.',
		'embl':'The EMBL flat file format. Uses Bio.GenBank internally.',
		'fasta':'The generic sequence file format where each record starts withan identifer line starting with a ">" character, followed bylines of sequence.',
		'fasta-2line':'Stricter interpretation of the FASTA format using exactlytwo lines per record (no line wrapping).',
		'fastq':'A "FASTA like" format used by Sanger which also stores PHREDsequence quality values (with an ASCII offset of 33).',
		'fastq-sanger':'An alias for "fastq" for consistency with BioPerl and EMBOSS',
		'fastq-solexa':'Original Solexa/Illumnia variant of the FASTQ format whichencodes Solexa quality scores (not PHRED quality scores) with anASCII offset of 64.',
		'fastq-illumina':'Solexa/Illumina 1.3 to 1.7 variant of the FASTQ formatwhich encodes PHRED quality scores with an ASCII offset of 64(not 33). Note as of version 1.8 of the CASAVA pipeline Illuminawill produce FASTQ files using the standard Sanger encoding.',
		'genbank':'The GenBank or GenPept flat file format.',
		'gb':'An alias for "genbank", for consistency with NCBI Entrez Utilities',
		'ig':'The IntelliGenetics file format, apparently the same as theMASE alignment format.',
		'imgt':'An EMBL like format from IMGT where the feature tables are moreindented to allow for longer feature types.',
		'pdb-seqres':'Reads a Protein Data Bank (PDB) file to determine thecomplete protein sequence as it appears in the header (no dependencies).',
		'pdb-atom':'Uses Bio.PDB to determine the (partial) protein sequence asit appears in the structure based on the atom coordinate section of thefile (requires NumPy for Bio.PDB).',
		'phd':'Output from PHRED, used by PHRAP and CONSED for input.',
		'pir':'A "FASTA like" format introduced by the National BiomedicalResearch Foundation (NBRF) for the Protein Information Resource(PIR) database, now part of UniProt.',
		'seqxml':'SeqXML, simple XML format described in Schmitt et al (2011).',
		'sff':'Standard Flowgram Format (SFF), typical output from Roche 454.',
		'sff-trim':'Standard Flowgram Format (SFF) with given trimming applied.',
		'swiss':'Plain text Swiss-Prot aka UniProt format.',
		'tab':'Simple two column tab separated sequence files, where eachline holds a record\'s identifier and sequence. For example,this is used as by Aligent\'s eArray software when savingmicroarray probes in a minimal tab delimited text file.',
		'qual':'A "FASTA like" format holding PHRED quality values fromsequencing DNA, but no actual sequences (usually providedin separate FASTA files).',
		'uniprot-xml':'The UniProt XML format (replacement for the SwissProt plaintext format which we call "swiss")'
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

#Show items
SHOWITEM = {
		"codon":CODON,
		"format":FORMAT,
	};

#Show ites
def showitem(args):
	for itn in args.show:
		print("{} list: ".format(itn));
		for it in SHOWITEM[itn].keys():
			print("\t{}:\t{}".format(it, SHOWITEM[itn][it]));
