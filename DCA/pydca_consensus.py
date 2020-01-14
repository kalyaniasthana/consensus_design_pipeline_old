import sys
import os
sys.path.insert(1, '../')
from consensus import *
import pandas as pd
from collections import defaultdict
import numpy as np

mappings = {'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17, 'V': 18, 'W': 19, 'Y': 20, '-': 21}

def mfdca_compute_params(filename):
	cwd = 'mfdca compute_params protein ../temp_files/train_file.fasta --verbose'
	os.system(cwd)

def read_couplings():
	couplings_filename = 'DCA_output_train_file/couplings_train_file.txt'
	couplings = {}

	with open(couplings_filename, 'r') as fin:
		for line in fin:
			if line.startswith('#      Length of sequences in alignment data: '):
				line = line.strip('\n')
				line = line.split(':')
				loa = int(line[1])
			elif line.startswith('#'):
				continue
			else:
				line = line.strip('\n')
				line = line.split(',')
				site_1, site_2, residue_1, residue_2, coupling_value = int(line[0]), int(line[1]), int(line[2]), int(line[3]), float(line[4])
				tup = (site_1, site_2, residue_1, residue_2)
				couplings[tup] = coupling_value

	return couplings

def read_fields():
	fields_filename = 'DCA_output_train_file/fields_train_file.txt'
	fields = {}

	with open(fields_filename, 'r') as fin:
		for line in fin:
			if line.startswith('#'):
				continue
			else:
				line = line.strip('\n')
				line = line.split(',')
				site, residue, field_value = int(line[0]), int(line[1]), float(line[2])
				fields[(site, residue)] = field_value

	return fields

def train_test_partition(train_file, test_file, main_file):
	sequences, headers = fasta_to_list(main_file)
	nos = len(sequences)
	no_train = int(0.3*nos)
	no_test = nos - no_train
	counter = 0
	with open(main_file, 'r') as fin:
		with open(train_file, 'w') as fout_1:
			with open(test_file, 'w') as fout_2:
				for line in fin:
					if line.startswith('>'):
						counter += 1
					if counter <= no_train:
						fout_1.write(line)
					else:
						fout_2.write(line)

def energy_function(sequence, couplings, fields):
	negative_bh = 0
	sequence_length = len(sequence)
	for i in range(sequence_length):
		aa_i = sequence[i].upper()
		if aa_i == '-':
			continue
		aa_i_no = mappings[aa_i]
		negative_bh += fields[(i+1, aa_i_no)]
		for j in range(i+1, sequence_length):
			aa_j = sequence[j].upper()
			if aa_j == '-':
				continue
			aa_j_no = mappings[aa_j]
			negative_bh += couplings[(i+1, j+1, aa_i_no, aa_j_no)]
	return negative_bh

def main():
	filename = 'PF00131'
	main_file = '../refined_alignments/' + filename + '_refined.fasta'
	train_file = '../temp_files/train_file.fasta'
	test_file = '../temp_files/test_file.fasta'
	#train_test_partition(train_file, test_file, main_file)
	#mfdca_compute_params(train_file)
	couplings = read_couplings()
	#print(couplings)
	fields = read_fields()
	test_sequences, test_headers = fasta_to_list(test_file)
	sequence_energies_refined_alignment = []
	for sequence in test_sequences:
		e = energy_function(sequence, couplings, fields)
		sequence_energies_refined_alignment.append(e)

	print(sequence_energies_refined_alignment)

if __name__ == '__main__':
	main()


