import gzip
import os
from Bio import SeqIO
import sys
from collections import Counter

amino_acids = 'ACDEFGHIKLMNPQRSTVWXY'
amino_acids = list(amino_acids) + ['-']

FGF_consensus_pdb = 'MRLRRLYCRTGGFHLQILPDGRVDGTREDNSPYSLLEIRAVEVGVVAIKGVKSGRYLAMNKKGRLYGSKHFTDECKFKERLLENGYNTYSSAKYRRGWYVALNKNGRPKKGNRTRRTQKATHFLPLPVSG'

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
'''
def gzip_to_fasta(gzfile, wf):
	write_file = open(wf, 'w')
	with gzip.open(gzfile, 'rb') as f:
		for line in f:
			line = line.decode('utf-8')
			write_file.write(line)
'''
#fasta file to clustalo CLI using os module
def fasta_to_clustalo(in_file, out_file):
	cmd = 'clustalo -i ' + in_file + ' -o ' + out_file + ' --force -v'
	os.system(cmd)

#converting fasta file into two lists (sequence list and name/header list)
def fasta_to_list(out_file):
	sequences = []
	name_list = []

	for record in SeqIO.parse(out_file, 'fasta'):
		name_list.append(record.id)
		sequences.append(str(record.seq))

	return sequences, name_list

#finding profile matrix for sequences in list format
def profile_matrix(sequences, pseudocount = 1):
	sequence_length = len(sequences[0])
	profile_matrix = {}
	for acid in amino_acids:
		profile_matrix[acid] = [pseudocount for i in range(sequence_length)]

	for i in range(len(sequences)):
		seq = sequences[i]
		for j in range(len(seq)):
			profile_matrix[seq[j]][j] += 1

	for aa in profile_matrix:
		l = profile_matrix[aa]
		for i in range(len(l)):
			if pseudocount > 0:
				l[i] /= (len(sequences)*(1 + pseudocount))
			else:
				l[i] /= len(sequences)

	return profile_matrix

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

	bad_sequence_set = set(bad_sequence_numbers)
	sequences = [x for i, x in enumerate(sequences) if i not in bad_sequence_set]
	name_list = [x for i, x in enumerate(name_list) if i not in bad_sequence_set]

	return sequences, name_list

#write sequence list and main list to fasta file
def list_to_fasta(sequences, name_list, fasta_file):
	file = open(fasta_file, 'w')
	for i in range(len(sequences)):
		file.write('>' + name_list[i] + '\n' + sequences[i] + '\n')

	file.close()

#remove dashes from sequences in fasta file and write to another fasta file
def remove_dashes(fasta_file_from, fasta_file_to):
	with open(fasta_file_from) as fin, open(fasta_file_to, 'w') as fout:
		for line in fin:
			fout.write(line.translate(str.maketrans('', '', '-')))

#find consensus sequence from sequences in list format
def consensus_sequence(sequences):
	consensus_seq = ''
	pm = profile_matrix(sequences, 1)
	sequence_length = len(sequences[0])
	for i in range(sequence_length):
		l = []
		for aa in pm:
			l.append(pm[aa][i])
		index = l.index(max(l))
		if amino_acids[index] == '-':
			if l[index] < 0.5:
				index = l.index(second_largest(l))
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

def median_of_list(sequence_lengths):
	n = len(sequence_lengths) 
	sequence_lengths.sort() 
	if n % 2 == 0:
		median1 = sequence_lengths[n//2]
		median2 = sequence_lengths[n//2 - 1] 
		median = (median1 + median2)/2
	else:
		median = sequence_lengths[n//2] 

	return median

def mean_of_list(sequence_lengths):
	n = len(sequence_lengths) 
	get_sum = sum(sequence_lengths)
	mean = get_sum / n
	return mean

#returns True only if all sequence length are within +/-20% of the initial mode length
def cut_off(sequence_lengths, mode):
	high = mode*1.2
	low = mode*0.8
	flag = True
	for length in sequence_lengths:
		if length > high or length < low:
			flag = False
			break
		
	return flag

def copy_file(in_file, out_file):
	with open(in_file) as fin, open(out_file, 'w') as fout:
		for line in fin:
			fout.write(line)
'''
def main():

	write_file = 'write.fasta'
	out_file = 'output.fasta'
	temp_file = 'temp.fasta'

	remove_dashes(write_file, temp_file)
	pm = {}
	sequence_lengths = sequence_length_list(temp_file)

	mode = mode_of_list(sequence_lengths)[0]
	median = median_of_list(sequence_lengths)
	mean = mean_of_list(sequence_lengths)

	flag = cut_off(sequence_lengths, mode)

	print('MODE | MEDIAN | MEAN' + '#'*100)

	print(str(mode) + '|' + str(median) + '|' + str(mean))
	sequences = []
	name_list = []
	iteration = 1

	while flag is False:

		print("ITERATION NUMBER: " + str(iteration) + '#'*100)

		fasta_to_clustalo(temp_file, out_file)

		sequences, name_list = fasta_to_list(out_file)
		pm = profile_matrix(sequences, 1)
		bad_sequence_numbers = find_bad_sequences(pm, sequences, name_list)
		sequences, name_list = remove_bad_sequences(sequences, name_list, bad_sequence_numbers)

		list_to_fasta(sequences, name_list, write_file)
		remove_dashes(write_file, temp_file)

		sequence_lengths = sequence_length_list(temp_file)
		print(sequence_lengths, '%'*100)
		print(len(sequence_lengths), 'number of sequences !!!!!!!!!')
		flag = cut_off(sequence_lengths, mode)
		print(flag, '   FLAG'*20)	
		print('MEDIAN | MEAN' + '#'*100)
		median = median_of_list(sequence_lengths)
		mean = mean_of_list(sequence_lengths)
		print(str(median) + '|' + str(mean))
		#print(consensus_sequence(sequences))
		print(len(consensus_sequence(sequences)), '$'*100)

		iteration += 1

	print(consensus_sequence(sequences))
'''
def main():

	write_file = 'write.fasta'
	copy_file('PF00167_full.fasta', write_file)
	#write_file is already aligned
	out_file = 'output.fasta'
	temp_file = 'temp.fasta'

	while True:
		#convert aligned write_file to list
		sequences, name_list = fasta_to_list(write_file)
		print(len(sequences), '#'*50)
		#profile matrix 
		pm = profile_matrix(sequences, 1)
		#find bad sequence indices in sequences list
		bad_sequence_numbers = find_bad_sequences(pm, sequences, name_list)
		#remove bad sequences
		sequences, name_list = remove_bad_sequences(sequences, name_list, bad_sequence_numbers)
		print(len(sequences), '$'*50)
		#convert new/smaller sequence and name list to fasta file
		list_to_fasta(sequences, name_list, temp_file)
		#remove dashes to make file ready for new alignment
		remove_dashes(temp_file, out_file)
		#new alignment 
		fasta_to_clustalo(out_file, write_file)



if __name__ == '__main__':
    main()