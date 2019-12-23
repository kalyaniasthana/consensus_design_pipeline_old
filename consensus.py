import gzip
import os
from Bio import SeqIO
import sys
from collections import Counter, OrderedDict
import matplotlib.pyplot as plt

amino_acids = 'ACDEFGHIKLMNPQRSTVWXY'
amino_acids = ['-'] + list(amino_acids)

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
def profile_matrix(sequences):
	sequence_length = len(sequences[0])
	profile_matrix = {}
	for acid in amino_acids:
		profile_matrix[acid] = [float(0) for i in range(sequence_length)]

	for i in range(len(sequences)):
		seq = sequences[i]
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
		file.write('>' + name_list[i] + '\n' + sequences[i] + '\n')

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
def consensus_sequence(sequences):
	consensus_seq = ''
	pm = profile_matrix(sequences)
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
'''
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
def sequence_length_without_dashes(sequence):
	new_seq = ''
	for i in range(len(sequence)):
		if sequence[i] == '-':
			continue
		else:
			new_seq += sequence[i]

	return len(new_seq)
'''
def consensus_length_cut_off(sequences, mode):
	flag = False
	cs = len(consensus_sequence(sequences))
	high = mode*1.05
	low = mode*0.95
	if cs > low and cs < high:
		flag = True

	return flag
'''
def selex_to_fasta(in_file, out_file):
	with open(in_file) as fin, open(out_file, 'w') as fout:
		headers = []
		sequences = []
		for line in fin:
			fout.write('>' + line[0:30] + '\n')
			fout.write(line[30: ])

def cdhit(in_file, out_file):
	cmd = 'cd-hit -i ' + in_file + ' -o ' + out_file + ' -T 1 -c 0.90'
	os.system(cmd)

def check_fasta(file):
	with open(file, "r") as handle:
		fasta = SeqIO.parse(handle, "fasta")
		return any(fasta)

def get_all_indices(l, value):

	return [i for i, val in enumerate(l) if val == value]

def main():

	#0th iteratin 

	write_file = 'write.fasta'
	out_file = 'output.fasta'
	temp_file = 'temp.fasta'
	bad_sequences = 'bad_sequences.fasta'

	selex_to_fasta('PF00167_full.txt', temp_file)
	remove_dashes(temp_file, write_file)
	cdhit(write_file, out_file)

	sequence_lengths = sequence_length_list(out_file)
	x = [i for i in range(len(sequence_lengths))]
	plt.scatter(x, sequence_lengths)
	plt.savefig('length_distribution.png')

	mode = mode_of_list(sequence_lengths)[0]
	fasta_to_clustalo(out_file, write_file)
	
	iteration = 1

	while True:
		print("ITERATION: " + str(iteration) + '*'*30)
		#convert aligned write_file to list
		sequences, name_list = fasta_to_list(write_file)
		pm = profile_matrix(sequences)
		#print(pm) #error here because these sequences are not aligned, FIX THIS
		print(consensus_sequence(sequences), len(consensus_sequence(sequences)), 'CONSENSUS!!!!!')
		#profile matrix 
		#pm = profile_matrix(sequences)
		#find bad sequence indices in sequences list
		bad_sequence_numbers = find_bad_sequences(pm, sequences, name_list)
		#remove bad sequences
		sequences, name_list = remove_bad_sequences(sequences, name_list, bad_sequence_numbers)
		#convert new/smaller sequence and name list to fasta file
		#print(consensus_sequence(sequences), 'CONSENSUS!!!!!')
		list_to_fasta(sequences, name_list, temp_file)
		#remove dashes to make file ready for new alignment
		remove_dashes(temp_file, out_file)
		#find new cut off
		#sequence_lengths = sequence_length_list(write_file)

		#flag = consensus_length_cut_off(sequences, mode)
		#new alignment
		fasta_to_clustalo(out_file, write_file)
		iteration += 1

if __name__ == '__main__':
    main()