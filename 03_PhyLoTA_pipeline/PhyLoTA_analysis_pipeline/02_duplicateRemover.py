#!/usr/bin/python

from Bio import SeqIO
import os

os.chdir("../01_FASTA")

def duplicate_remover(filename):
	taxon_list = []
	unique_sequences = []
	new_file = open(filename + ".output", 'w')
	for seq_record in SeqIO.parse(filename, "fasta"):
		if seq_record.description.split("_")[1] not in taxon_list:
			taxon_list.append(seq_record.description.split("_")[1])
			unique_sequences.append(seq_record)
	SeqIO.write(unique_sequences, new_file, "fasta")

def ls(path):
	files = os.listdir(path)
	return files

files = ls(".")
print files
for file in files:
	if file.endswith(".fasta"):
		duplicate_remover(str(file))