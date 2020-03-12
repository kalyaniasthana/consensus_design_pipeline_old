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
from random import randint, choice
import signal
import subprocess
from os import listdir
from os.path import isfile, join
import time
from datetime import datetime

start = time.time()

mappings = {'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,
'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17, 'V': 18, 'W': 19, 'Y': 20, '-': 21}

#mean field dca (pydca)
def mfdca_compute_params(filename):
	cwd = 'mfdca compute_params protein ../temp_files/train_file.fasta --verbose'
	cwd = cwd.split(' ')

	p = subprocess.Popen(cwd, stdout= subprocess.PIPE, stderr = subprocess.STDOUT)

	while True:
		line = p.stdout.readline()
		if not line:
			break
		line = str(line, 'utf-8')
		line = line.strip('\n')
		print(line)


		if 'effective number of sequences' in line:
			line = line.split('effective number of sequences: ')
			neffective = float(line[1])
			if neffective < 150:
				print('Effective number of sequences is too low')
				time.sleep(5)
				break

#reading couplings file
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

#reading fields file
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

#train test partition of the refined alignment 
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

#energy calculation using couplings and fields from mfdca (pydca)
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
			energy -= val/2
	
	last_aa = sequence[sequence_length - 1].upper()
	if last_aa != '-':
		if last_aa in mappings:
			last_aa_no = mappings[last_aa]
			val = fields[(sequence_length, last_aa_no)]

	energy += val

	'''
	for i in range(sequence_length):
		aa_i = sequence[i].upper()
		if aa_i == '-':
			continue
		if aa_i not in mappings:
			continue
		aa_i_no = mappings[aa_i]
		try:
			val = fields[(i+1, aa_i_no)]
			energy -= val
		except:
			continue

		for j in range(sequence_length):
			aa_j = sequence[j].upper()
			if i != j:
				if aa_j == '-':
					continue

				if aa_j not in mappings:
					continue

				aa_j_no = mappings[aa_j]
				try:
					val = couplings[(i+1, j+1, aa_i_no, aa_j_no)]
					energy -= val
				except:
					continue
		'''

	return energy

'''
def score_sequence(energy):
	score = -energy
	return score
'''
#split combined alignment to refined and hmm

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
'''
def score_loop(energies):
	scores = []
	for e in energy:
		scores.append(score_sequence(e))

	scores = np.array(scores)

	return scores
'''

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

#shuffing sequences
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

#shuffling sequences w/o moving dashes
def fisher_yates_without_dashes(sequences):
	seq_length = len(sequences[0])
	shuffled = []
	for sequence in sequences:
		seq = list(sequence)
		rand_list = []
		for i in range(seq_length):
			if seq[i] == '-':
				pass
			else:
				rand_list.append(i)
		for i in rand_list:
			j = choice(rand_list)
			temp = seq[i]
			seq[j] = seq[i]
			seq[j] = temp

		shuffled.append(''.join(seq))

	return shuffled

#plot energies
def plot_energies(
        sequence_energies_from_training_sequences, sequence_energies_from_hmm_alignment,
        consensus_energy, hmm_consensus_energy, bins_refined, bins_hmm, dca_energy_plot
    ):

	plt.hist(sequence_energies_from_training_sequences, alpha = 0.5, edgecolor = 'black', label = 'sequences from refined MSA', density = True)
	plt.hist(sequence_energies_from_hmm_alignment, alpha = 0.5, edgecolor = 'black', label = 'sequences emitted from profile hmm', density = True)
	#plt.hist(sequence_energies_from_shuffled_sequences, alpha = 0.5, edgecolor = 'black', label = 'shuffled sequences')
	plt.axvline(x = consensus_energy, color = 'red', label = 'consensus sequence from refined MSA')
	plt.axvline(x = hmm_consensus_energy, color = 'blue', label = 'consensus sequence from hmm sequences')
	plt.legend(loc = 'upper right')
	plt.xlabel('DCA energies')
	plt.savefig(dca_energy_plot)
	plt.clf()
	plt.cla()
	plt.close()

#for martin's program
def write_matlab_script(filename):

	dca_calculation_script = '../martin_dca/dca_energy.m'
	with open(dca_calculation_script, 'w') as f:
		f.write("try\n")
		f.write("\tfasta_test = '/media/Data/consensus/temp_files/only_hmm_emitted.fasta';\n")
		f.write("\tfasta_train = '/media/Data/consensus/temp_files/train_file.fasta';\n")
		f.write("\taccession = " + "'" + filename + "';\n")
		f.write("\t[score_train, score_test] = calculate_dca_scores(fasta_train,fasta_test,accession)\n")
		f.write("catch\n")
		f.write("\texit\n")
		f.write("end")

def call_matlab_script():
	os.chdir('/usr/local/MATLAB/R2016a/bin')
	cwd = './matlab -softwareopengl -nodesktop -r "run(' + "'/media/Data/consensus/martin_dca/dca_energy.m')" + ';exit;"'
	os.system(cwd)
	#os.chdir('../DCA')

def main_pydca(filename, accession_file):

	combined_file, train_file, test_file, only_refined, only_hmm, dca_energy_plot, consensus_file, combined_with_consensus = pydca_strings(filename)
	try:
		if os.stat(consensus_file).st_size == 0 or os.stat(combined_file) == 0:
			remove_accession(accession_file, filename)
			print('Some files are missing, so skipping this family for now!')
			return

		if path.exists(dca_energy_plot) and os.stat(dca_energy_plot).st_size != 0:
			remove_accession(accession_file, filename)
			print('Already calculated!')
			return

	except Exception as e:
		print(e)
		print('Skipping this family for now')
		return

	split_combined_alignment(combined_file, only_refined, only_hmm)
	train_test_partition(train_file, test_file, only_refined)


	#try:
	#	mfdca_compute_params(train_file)
	#except Exception as e:
	#	print(e)
	#	print('excuse me miss, we have a problem here.')
	#	return

	#couplings, loa = read_couplings()
	#fields = read_fields()


	try:

		hmm_sequences, hmm_headers = fasta_to_list(only_hmm)
		#sequence_energies_from_hmm_alignment = sequence_energies_loop(hmm_sequences, couplings, fields)
		#print(sequence_energies_from_hmm_alignment)
		#time.sleep(3)
		#print('\n')

		training_sequences, training_headers = fasta_to_list(train_file)
		#sequence_energies_from_training_sequences = sequence_energies_loop(training_sequences, couplings, fields)
		#print(sequence_energies_from_training_sequences)
		#time.sleep(3)

	except Exception as e:
		print(e)
		return

	hmm_pm = profile_matrix(hmm_sequences)
	hmm_consensus = consensus_sequence(hmm_sequences, hmm_pm)
	hmm_header = '>consensus-from-hmm-emitted-sequences'

	cons_seqs, cons_headers = fasta_to_list(consensus_file)
	if hmm_header[1: ] not in cons_headers:
		with open(consensus_file, 'a') as fin:
			fin.write(hmm_header)
			fin.write('\n')
			fin.write(hmm_consensus)
			fin.write('\n')

	consensus_seq, consensus_header = fasta_to_list(consensus_file)
	consensus_seq = consensus_seq[0]

	option = '2'
	realign(option, combined_file, consensus_file, combined_with_consensus)

	for record in SeqIO.parse(combined_with_consensus, 'fasta'):
		#print(record.id)
		if record.id == 'consensus-from-hmm-emitted-sequences':
			hmm_consensus_aligned = str(record.seq)
		elif record.id == 'consensus-from-refined-alignment':
			consensus_seq_aligned = str(record.seq)

	#consensus_energy = energy_function(consensus_seq_aligned, couplings, fields)
	#hmm_consensus_energy = energy_function(hmm_consensus_aligned, couplings, fields)
	#print('\n')
	#print(consensus_energy, hmm_consensus_energy)

	#print(consensus_seq_aligned, type(consensus_seq_aligned))
	#print('\n')
	#print(hmm_consensus_aligned, type(hmm_consensus_aligned))
	#rewrite both consensus in the consensus file in the aligned form
	with open(consensus_file, 'w') as fin:
		fin.write('>consensus-from-refined-alignment\n')
		fin.write(consensus_seq_aligned + '\n')
		fin.write('>consensus-from-hmm-emitted-sequences\n')
		fin.write(hmm_consensus_aligned)
	#shuffled_sequences = fisher_yates_without_dashes(training_sequences)
	#sequence_energies_from_shuffled_sequences = sequence_energies_loop(shuffled_sequences, couplings, fields)

	#minimum = min([min(sequence_energies_from_hmm_alignment), min(sequence_energies_from_training_sequences)]) - 1000
	#maximum = max([max(sequence_energies_from_hmm_alignment), max(sequence_energies_from_training_sequences)]) + 1000
	#bins_refined = np.linspace(min(sequence_energies_from_training_sequences), max(sequence_energies_from_training_sequences))
	#bins_hmm = np.linspace(min(sequence_energies_from_hmm_alignment), max(sequence_energies_from_hmm_alignment))
	pi = percentage_identity(consensus_seq_aligned, hmm_consensus_aligned)
	print('Percentage Identity of the two consensus sequences: ', pi, '\n')

	#plot_energies(sequence_energies_from_training_sequences, sequence_energies_from_hmm_alignment, 
	#	consensus_energy, hmm_consensus_energy, bins_refined, bins_hmm, dca_energy_plot)

	try:
		write_matlab_script(filename)
		call_matlab_script()
	except:
		return

	end = time.time() - start
	done_file = '/media/Data/consensus/temp_files/done_files.txt'
	print('It took ' + str(end) + ' seconds to do DCA calculation for: ' + filename)
	now = datetime.now()
	# dd/mm/YY H:M:S
	dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
	with open(done_file, 'a') as f:
		f.write(filename + ' ' + dt_string + '\n')
	print('\n\n\n')
	time.sleep(2)

