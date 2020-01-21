import sys
import os
sys.path.insert(1, '../')
from consensus import *
import pandas as pd
from collections import defaultdict
import numpy as np
from Bio import SeqIO
from matplotlib import pyplot as plt

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
	no_train = int(0.2*nos)
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
		if aa_i not in mappings:
			continue
		aa_i_no = mappings[aa_i]
		negative_bh += fields[(i+1, aa_i_no)]
		for j in range(i+1, sequence_length):
			aa_j = sequence[j].upper()
			if aa_j == '-':
				continue
			if aa_j not in mappings:
				continue
			aa_j_no = mappings[aa_j]
			negative_bh += couplings[(i+1, j+1, aa_i_no, aa_j_no)]

	return negative_bh

def split_combined_alignment(combined_alignment, only_refined, only_hmm):
	refined = []
	hmm = []
	with open(only_refined, 'w') as fout_1:
		with open(only_hmm, 'w') as fout_2:
			for record in SeqIO.parse(combined_alignment, 'fasta'):
				if 'refined' in record.id:
					hmm.append(record)
				else:
					refined.append(record)

	SeqIO.write(refined, only_refined, 'fasta')
	SeqIO.write(hmm, only_hmm, 'fasta')

def main():
	filename = 'PF00167'
	combined_file = '../combined_alignments/' + filename + '_combined.fasta'
	train_file = '../temp_files/train_file.fasta'
	test_file = '../temp_files/test_file.fasta'
	only_refined = '../temp_files/only_refined_sequences.fasta'
	only_hmm = '../temp_files/only_hmm_emitted.fasta'
	dca_energy_plot = '../dca_energy_plots/' + filename + '_dca_energies.png'

	split_combined_alignment(combined_file, only_refined, only_hmm)
	train_test_partition(train_file, test_file, only_refined)
	
	mfdca_compute_params(train_file)
	couplings = read_couplings()
	fields = read_fields()

	test_sequences, test_headers = fasta_to_list(test_file)
	test_sequence_energies_refined_alignment = []
	for sequence in test_sequences:
		e = energy_function(sequence, couplings, fields)
		test_sequence_energies_refined_alignment.append(e)

	print(test_sequence_energies_refined_alignment)
	print('\n\n')
	hmm_sequences, hmm_headers = fasta_to_list(only_hmm)
	sequence_energies_from_hmm_alignment = []
	for sequence in hmm_sequences:
		e = energy_function(sequence, couplings, fields)
		sequence_energies_from_hmm_alignment.append(e)

	print(sequence_energies_from_hmm_alignment)
	print('\n\n')

	bins = np.linspace(-5000, 1000)
	plt.hist(test_sequence_energies_refined_alignment, bins, alpha = 0.5, label = 'test sequences from refined MSA')
	plt.hist(sequence_energies_from_hmm_alignment, bins, alpha = 0.5, label = 'sequences emitted from profile hmm')
	plt.legend(loc = 'upper right')
	plt.savefig(dca_energy_plot)
	plt.show()

if __name__ == '__main__':
	main()


