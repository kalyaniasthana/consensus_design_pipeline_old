import sys
import os
sys.path.insert(1, '../')
from consensus import *
from ugly_strings import *
import pandas as pd
from collections import defaultdict
import numpy as np
from Bio import SeqIO
from matplotlib import pyplot as plt
from random import randint

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

	return couplings, loa

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
	no_train = int(1*nos)
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
	energy = 0
	sequence_length = len(sequence)

	for i in range(sequence_length - 1):
		aa_i = sequence[i].upper()
		if aa_i == '-':
			continue
		if aa_i not in mappings:
			continue
		aa_i_no = mappings[aa_i]
		val = fields[(i+1, aa_i_no)]
		energy += val

		for j in range(i+1, sequence_length):
			aa_j = sequence[j].upper()
			if aa_j == '-':
				continue
			if aa_j not in mappings:
				continue

			aa_j_no = mappings[aa_j]
			val = couplings[(i+1, j+1, aa_i_no, aa_j_no)]
			energy -= val
	
	last_aa = sequence[sequence_length - 1].upper()
	if last_aa != '-':
		if last_aa in mappings:
			last_aa_no = mappings[last_aa]
			val = fields[(sequence_length, last_aa_no)]

	energy += val
	return energy

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

def sequence_energies_loop(sequences, couplings, fields):
	energy_list = []
	for sequence in sequences:
		e = energy_function(sequence, couplings, fields)
		energy_list.append(e)

	return energy_list

def percentage_identity(refined_consensus, hmm_consensus):
	matches = 0
	seq_length = len(refined_consensus)
	for i in range(seq_length):
		if refined_consensus[i] == hmm_consensus[i]:
			matches += 1

	pi = (matches*100)/seq_length
	return pi

def analyse_couplings(loa, couplings):

	for i in range(1, loa):
		values = []
		for j in range(i+1, loa + 1):
			for k in range(1, 21):
				for l in range(1, 21):
					values.append(couplings[(i, j, k , l)])
		max_value = max(values)
		min_value = min(values)
		print(i, max_value, min_value)

def analyse_fields(loa, fields):

	for i in range(1, loa + 1):
		values = []
		for j in range(1, 21):
			values.append(fields[(i, j)])
		max_value = max(values)
		min_value = min(values)
		print(i, max_value, min_value)

def fisher_yates_shuffling(sequences):
	seq_length = len(sequences[0])
	shuffled = []
	for sequence in sequences:
		seq = list(sequence)
		for i in range(seq_length):
			j = randint(0, i)
			temp = seq[i]
			seq[i] = seq[j]
			seq[j] = temp

		shuffled.append(''.join(seq))

	return shuffled

def main():
	filename = 'PF00131'
	combined_file, train_file, test_file, only_refined, only_hmm, dca_energy_plot, consensus_file, combined_with_consensus = pydca_strings(filename)

	split_combined_alignment(combined_file, only_refined, only_hmm)
	train_test_partition(train_file, test_file, only_refined)
	mfdca_compute_params(train_file)

	couplings, loa = read_couplings()
	fields = read_fields()
	#analyse_fields(loa, fields)
	#sys.exit()

	#analyse_couplings(loa, couplings)
	#sys.exit()

	hmm_sequences, hmm_headers = fasta_to_list(only_hmm)
	sequence_energies_from_hmm_alignment = sequence_energies_loop(hmm_sequences, couplings, fields)
	print(sequence_energies_from_hmm_alignment)
	print('\n\n')

	training_sequences, training_headers = fasta_to_list(train_file)
	sequence_energies_from_training_sequences = sequence_energies_loop(training_sequences, couplings, fields)
	print(sequence_energies_from_training_sequences)
	print('\n\n')

	consensus_seq, consensus_header = fasta_to_list(consensus_file)
	consensus_seq = consensus_seq[0]
	#consensus_energy = energy_function(consensus_seq, couplings, fields)

	hmm_pm = profile_matrix(hmm_sequences)
	hmm_consensus = consensus_sequence(hmm_sequences, hmm_pm)
	#hmm_consensus_energy = energy_function(hmm_consensus, couplings, fields)
	hmm_header = '>consensus-from-hmm-emitted-sequences'

	cons_seqs, cons_headers = fasta_to_list(consensus_file)
	if hmm_header[1: ] not in cons_headers:
		with open(consensus_file, 'a') as fin:
			fin.write(hmm_header)
			fin.write('\n')
			fin.write(hmm_consensus)
			fin.write('\n')

	option = '2'
	realign(option, combined_file, consensus_file, combined_with_consensus)

	for record in SeqIO.parse(combined_with_consensus, 'fasta'):
		#print(record.id)
		if record.id == 'consensus-from-hmm-emitted-sequences':
			hmm_consensus_aligned = record.seq
		elif record.id == 'consensus-from-refined-alignment':
			consensus_seq_aligned = record.seq

	consensus_energy = energy_function(consensus_seq_aligned, couplings, fields)
	hmm_consensus_energy = energy_function(hmm_consensus_aligned, couplings, fields)

	'''
	random_sequences_file = '../temp_files/random_sequences.fasta'
	random_seq, random_headers = fasta_to_list(random_sequences_file)
	sequence_energies_random_sequences = sequence_energies_loop(random_seq, couplings, fields)
	print(sequence_energies_random_sequences)
	print('\n\n')
	'''

	shuffled_sequences = fisher_yates_shuffling(training_sequences)
	sequence_energies_from_shuffled_sequences = sequence_energies_loop(shuffled_sequences, couplings, fields)
	print(sequence_energies_from_shuffled_sequences)
	print('\n\n')

	minimum = min([min(sequence_energies_from_hmm_alignment), min(sequence_energies_from_training_sequences), min(sequence_energies_from_shuffled_sequences)]) - 1000
	maximum = max([max(sequence_energies_from_hmm_alignment), max(sequence_energies_from_training_sequences), max(sequence_energies_from_shuffled_sequences)]) + 1000
	bins = np.linspace(minimum, maximum)
	pi = percentage_identity(consensus_seq_aligned, hmm_consensus_aligned)
	print('Percentage Identity of the two consensus sequences: ', pi, '\n')

	plt.hist(sequence_energies_from_training_sequences, alpha = 0.5, edgecolor = 'black', label = 'sequences from refined MSA')
	plt.hist(sequence_energies_from_hmm_alignment, bins, alpha = 0.5, edgecolor = 'black', label = 'sequences emitted from profile hmm')
	#plt.hist(sequence_energies_random_sequences, alpha = 0.5, edgecolor = 'black', label = 'sequence energies from random sequences')
	plt.hist(sequence_energies_from_shuffled_sequences, alpha = 0.5, edgecolor = 'black', label = 'shuffled sequences')
	plt.axvline(x = consensus_energy, color = 'red', label = 'consensus sequence from refined MSA')
	plt.axvline(x = hmm_consensus_energy, color = 'blue', label = 'consensus sequence from hmm sequences')
	plt.legend(loc = 'upper right')
	plt.xlabel('DCA energies')
	plt.savefig(dca_energy_plot)

if __name__ == '__main__':
	main()

