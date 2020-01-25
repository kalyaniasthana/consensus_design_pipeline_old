import gzip
import os
from Bio import SeqIO, AlignIO
import sys
from collections import Counter, OrderedDict
import matplotlib.pyplot as plt
import time
import copy
from shutil import copyfile
import json
from os import path
import pandas as pd
from ugly_strings import *
import subprocess

start = time.time()

amino_acids = 'ACDEFGHIKLMNPQRSTVWXY'
amino_acids = ['-'] + list(amino_acids)

#finding second largest number in a list
def second_largest(numbers):
	count = 0
	m1 = m2 = float('-inf')
	for x in numbers:
		count += 1
		if x > m2:
			if x >= m1:
				m1, m2 = x, m1
			else:
				m2 = x

	return m2 if count >= 2 else None

#fasta file to clustalo CLI using os module
def fasta_to_clustalo(in_file, out_file):
	cmd = 'clustalo -i ' + in_file + ' -o ' + out_file + ' --force -v'
	os.system(cmd)

#fasta file to MAFFT
def fasta_to_mafft(in_file, out_file):
	cmd = 'mafft ' + in_file + ' > ' + out_file
	os.system(cmd)

#fasta file to MUSCLE
def fasta_to_muscle(in_file, out_file):
	cmd = 'muscle -in ' + in_file + ' -out ' + out_file
	os.system(cmd)

#converting fasta file into two lists (sequence list and name/header list)
def fasta_to_list(out_file):
	sequences = []
	name_list = []

	for record in SeqIO.parse(out_file, 'fasta'):
		name_list.append(record.id)
		sequences.append(str(record.seq).upper())

	return sequences, name_list

#finding profile matrix for sequences in list format
def profile_matrix(sequences):
	sequence_length = len(sequences[0])
	profile_matrix = {}
	for acid in amino_acids:
		profile_matrix[acid] = [float(0) for i in range(sequence_length)]

	for i in range(len(sequences)):
		seq = sequences[i].upper()
		for j in range(len(seq)):
			profile_matrix[seq[j]][j] += float(1)

	for aa in profile_matrix:
		l = profile_matrix[aa]
		for i in range(len(l)):
			l[i] /= float(len(sequences))

	pm = OrderedDict([(x, profile_matrix[x]) for x in amino_acids])

	return pm

#finding index of bad sequence numbers in the sequence list
def find_bad_sequences(profile_matrix, sequences, name_list):
	max_value = max(profile_matrix['-'])

	if max_value == 1:
		max_value = second_largest(profile_matrix['-'])

	positions = []
	for i in range(len(profile_matrix['-'])):
		if profile_matrix['-'][i] == max_value:
			positions.append(i)

	bad_sequence_numbers = []
	for i in range(len(sequences)):
		for position in positions:
			if sequences[i][position] != '-':
				if i not in bad_sequence_numbers:
					bad_sequence_numbers.append(i)

	return bad_sequence_numbers

#removing bad sequence numbers and returning new sequence list and name list
def remove_bad_sequences(sequences, name_list, bad_sequence_numbers):

	sequences = [x for i, x in enumerate(sequences) if i not in bad_sequence_numbers]
	name_list = [x for i, x in enumerate(name_list) if i not in bad_sequence_numbers]

	return sequences, name_list

#write sequence list and main list to fasta file
def list_to_fasta(sequences, name_list, fasta_file):
	file = open(fasta_file, 'w')
	for i in range(len(sequences)):
		file.write('>' + name_list[i] + '\n' + sequences[i].upper() + '\n')

	file.close()

#remove dashes from sequences in fasta file and write to another fasta file
def remove_dashes(fasta_file_from, fasta_file_to):
	with open(fasta_file_from) as fin, open(fasta_file_to, 'w') as fout:
		for line in fin:
			if line.startswith('>'):
				fout.write(line)
			else:
				fout.write(line.translate(str.maketrans('', '', '-')))

#find consensus sequence from sequences in list format
def consensus_sequence(sequences, pm):
	consensus_seq = ''
	#pm = profile_matrix(sequences)
	sequence_length = len(sequences[0])
	for i in range(sequence_length):
		l = []
		for aa in pm:
			l.append(pm[aa][i])
		max_value = max(l)
		indices = get_all_indices(l, max_value)
		index = indices[0]
		if amino_acids[index] == '-':
			if l[index] < 0.5:
				second_largest_value = second_largest(l)
				if second_largest_value == max_value:
					index = indices[1]
				else:
					index = l.index(second_largest_value)
			else:
				continue
		consensus_seq += amino_acids[index]

	return consensus_seq

# returning a list of sequence lengths
def sequence_length_list(read_file):
	sequences, name_list = fasta_to_list(read_file)
	sequence_lengths = []
	#print(len(sequences))
	for seq in sequences:
		sequence_lengths.append(len(seq))
	return sequence_lengths

def mode_of_list(sequence_lengths):
	n = len(sequence_lengths)
	data = Counter(sequence_lengths) 
	get_mode = dict(data) 
	mode = [k for k, v in get_mode.items() if v == max(list(data.values()))]
	if n == len(mode):
		return None
	else:
		return mode

def selex_to_fasta(in_file, out_file):
	with open(in_file) as fin, open(out_file, 'w') as fout:
		headers = []
		sequences = []
		for line in fin:
			fout.write('>' + line[0:30].upper() + '\n')
			fout.write(line[30: ].upper())

def cdhit(in_file, out_file):
	cmd = 'cd-hit -i ' + in_file + ' -o ' + out_file + ' -T 1 -c 0.90'
	os.system(cmd)

def get_all_indices(l, value):

	return [i for i, val in enumerate(l) if val == value]

def stockholm_to_fasta(ifile, ofile):
	with open(ifile, 'r') as fin:
		with open(ofile, 'w') as fout:
			sequences = SeqIO.parse(ifile, 'stockholm')
			SeqIO.write(sequences, ofile, 'fasta')
	os.system('rm -rf ' + ifile)

def fasta_to_plain(accession, filename):
	alignment = AlignIO.read(open(filename, 'fasta'))
	sequences = [record.seq for record in alignment]
	plain_file = 'temp_files/' + accession + '_refined_noheader.txt'
	with open(plain_file, 'w') as f:
		for seq in sequences:
			f.write(str(seq))
			f.write('\n')

#Percent Identity = (Matches x 100)/Length of aligned region (with gaps)
def percentage_identity(consensus_fasta):
	seqs, head = fasta_to_list(consensus_fasta)
	matches = 0
	seq_length = len(seqs[0])
	for i in range(seq_length):
		if seqs[0][i] == seqs[1][i]:
			matches += 1
	pi = (matches*100)/seq_length
	return pi

def store_retrieve_identity_dict(accession, pi, filename):
	my_dict = {}
	with open(filename, 'a') as f:
		f.write(accession + ':' + str(pi) + '\n')

	with open(filename) as f:
		for line in f:
			line = line.strip('\n')
			line = line.split(':')
			my_dict[line[0]] = float(line[1])

	return my_dict

def plot_dict_key_and_value(my_dict):
	accessions = list(my_dict.keys())
	pi_values = list(my_dict.values())
	pi_values = [int(i) for i in pi_values]
	df = pd.DataFrame({'Accession' : accessions, '% Identity' : pi_values})
	ax = df.plot.bar(x = 'Accession', y = '% Identity', rot = 0,  stacked = True, colormap = 'Paired')
	fig = ax.get_figure()
	fig.savefig('temp_files/pi_plot.png')

def alignment(option, in_file, out_file):
	if option == '1':
		fasta_to_clustalo(in_file, out_file)
	elif option == '2':
		fasta_to_mafft(in_file, out_file)
	elif option == '3':
		fasta_to_muscle(in_file, out_file)
	else:
		print('Invalid Option')
		sys.exit()

def realign(option, original_alignment, hmm_sequences, out_file):
	if option == '2':
		cwd = 'mafft --add ' + hmm_sequences + ' --reorder --keeplength ' + original_alignment + ' > ' + out_file
		#mafft --add new_sequences --reorder existing_alignment > output
		os.system(cwd)
		#op = subprocess.check_output(cwd, shell=True)
		#with open(out_file, 'w') as fin:
		#	for line in op:
		#		fin.write(line)
	else:
		print('Invalid input')
		sys.exit()

def refine_filename(ip):
	ip = str(ip, 'utf-8')
	ip = ip.strip('\n')
	ip = ip.replace('./', '')
	return ip

def main():

	#accession_list = []
	#with open('temp_files/accession_list.txt', 'r') as f:
	#	for line in f:
	#		accession_list.append(line.strip('\n'))

	accession_list = ['PF00167']

	print('1. Clustal Omega 2. MAFFT 3. MUSCLE')
	#option = input()
	option = '2' #using only MAFFT for now
	write_file, out_file, temp_file, perc_idens = common_files()

	for accession in accession_list:
		plot = 'dca_energies/' + accession + '_dca_energies.png'
		if path.exists(plot):
			continue
		my_file = 'pfam_entries/' + accession + '.fasta'
		if path.exists(my_file):

			#0th iteration 
			filename = accession
			print(filename, '%'*30)
			file = 'pfam_entries/' + filename + '.fasta'
			copyfile(file, temp_file)
			remove_dashes(temp_file, write_file)
			my_dict = {}

			#####################################
			#break condition for testing purposes
			#test_seq, test_head = fasta_to_list(write_file)
			#if len(test_seq) > 3000:
			#	continue
			#####################################

			try:
				cdhit(write_file, out_file)
			except Exception as e:
				print('Exception: ' + str(e))
				continue

			refined_alignment, plot, final_consensus, profile_hmm, hmm_emitted_sequences, combined_alignment = specific_files(filename)

			sequence_lengths = sequence_length_list(out_file)
			x = [i for i in range(len(sequence_lengths))]
			plt.scatter(x, sequence_lengths)
			plt.savefig(plot)

			mode = mode_of_list(sequence_lengths)[0]
			try:
				alignment(option, out_file, write_file)
			except Exception as e:
				print('Exception: ' + str(e))
				continue

			iteration = 1
			#exit conditions
			#if number of sequences < 100
			#if length of alignment does not change in subsequent iterations
			#if length of alignment becomes to small i.e -15 the desired length (mode length)
			loa = 0
			try:

				while True:
					print("Iteration Number: " + str(iteration) + '*'*30)
					sequences, name_list = fasta_to_list(write_file)
					number_of_sequences = len(sequences)
					length_of_alignment = len(sequences[0])
					print('Length of Alignment = ', length_of_alignment)
					print('Alignment length (previous iteration): ', loa, 'Alignment length (current iteration): ', length_of_alignment)

					if number_of_sequences < 100 or length_of_alignment < mode - 15 or loa == length_of_alignment:
						copyfile(write_file, refined_alignment)
						f = open(final_consensus, 'w')
						f.write('>consensus-from-refined-alignment' + '\n')
						f.write(cs + '\n')
						cwd = 'hmmbuild ' + profile_hmm + ' ' + refined_alignment
						os.system(cwd)
						N = number_of_sequences
						L = length_of_alignment
						cwd = 'hmmemit -N ' + str(N) + ' -o ' + hmm_emitted_sequences + '-L' + str(L) + ' ' + profile_hmm
						os.system(cwd)
						os.chdir('/media/Data/consensus/hmm_emitted_sequences')
						cwd = 'find | grep ' + filename + '_hmmsequences.fasta'
						ip = subprocess.check_output(cwd, shell=True)
						ip = refine_filename(ip)
						ip = 'hmm_emitted_sequences/' + ip
						os.chdir('/media/Data/consensus')
						realign(option, refined_alignment, ip, combined_alignment)
						break

					pm = profile_matrix(sequences)
					cs = consensus_sequence(sequences, pm)

					print(cs, len(cs), 'Consensus from refined alignment')

					bad_sequence_numbers = find_bad_sequences(pm, sequences, name_list)
					sequences, name_list = remove_bad_sequences(sequences, name_list, bad_sequence_numbers)
					list_to_fasta(sequences, name_list, temp_file)

					#remove dashes to make file ready for new alignment
					remove_dashes(temp_file, out_file)
					alignment(option, out_file, write_file)
					loa = copy.deepcopy(length_of_alignment)
					iteration += 1

			except Exception as e:
				print('Exception: ' + str(e))
				time.sleep(5)
				continue

			print('***********Final Consensus Sequence from refined alignment: ')
			print(cs)
		
			end = time.time() - start
			print('It took ' + str(end) + ' seconds to run the script')
			print('\n\n\n')
			print('Time for next protein family/domain/motif\n')
			time.sleep(5)

if __name__ == '__main__':
    main()
