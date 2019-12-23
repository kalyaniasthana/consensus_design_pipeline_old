import gzip
import os
from Bio import SeqIO
import sys
from collections import Counter, OrderedDict
import matplotlib.pyplot as plt
import time

start = time.time()

amino_acids = 'ACDEFGHIKLMNPQRSTVWXY'
amino_acids = ['-'] + list(amino_acids)

FGF_consensus_pdb = 'MRLRRLYCRTGGFHLQILPDGRVDGTREDNSPYSLLEIRAVEVGVVAIKGVKSGRYLAMNKKGRLYGSKHFTDECKFKERLLENGYNTYSSAKYRRGWYVALNKNGRPKKGNRTRRTQKATHFLPLPVSG'

def Nmaxelements(input_list, N): 
	final_list = []
	input_list.sort(reverse = True)
	for i in range(N):
		final_list.append(input_list[i])
	return final_list

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

	#0th iteration 

	write_file = 'write.fasta'
	out_file = 'output.fasta'
	temp_file = 'temp.fasta'
	bad_sequences = 'bad_sequences.fasta'
	#pfam 30.0 FGF 
	#selex_to_fasta('PF00167_full.txt', temp_file)
	#pfam 32.0 FGF
	#selex_to_fasta('PF00167_latest.txt', temp_file)
	#SH3 family
	selex_to_fasta('PF00018_full.txt', temp_file)
	remove_dashes(temp_file, write_file)
	cdhit(write_file, out_file)

	sequence_lengths = sequence_length_list(out_file)
	x = [i for i in range(len(sequence_lengths))]
	plt.scatter(x, sequence_lengths)
	plt.savefig('length_distribution.png')

	mode = mode_of_list(sequence_lengths)[0]
	fasta_to_clustalo(out_file, write_file)
	
	iteration = 1
	#exit conditions
	#if number of sequences < 100
	#if length of alignment does not change in subsequent iterations
	#if length of alignment becomes to small i.e -15 the desired length (mode length)

	while True:
		print("ITERATION: " + str(iteration) + '*'*30)
		#convert aligned write_file to list
		sequences, name_list = fasta_to_list(write_file)
		number_of_sequences = len(sequences)
		length_of_alignment = len(sequences[0])
		print('LENGTH OF ALIGNMENT = ', length_of_alignment)
		if number_of_sequences < 100 or length_of_alignment < mode - 15:
			break
		pm = profile_matrix(sequences)
		#print(pm) #error here because these sequences are not aligned, FIX THIS
		cs = consensus_sequence(sequences)
		print(cs, len(cs), 'CONSENSUS!!')
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
		temp_seqs, temp_names = fasta_to_list(write_file)
		loa = len(temp_seqs[0])
		if loa == length_of_alignment:
			break
		iteration += 1

	print('***********Final Consensus Sequence: ' + '\n')
	print(cs)
	end = time.time() - start
	print('It took ' + str(end) + ' seconds to run the script' )

'''
def main_re():
	my_list = [2, 6, 41, 85, 0, 3, 7, 6, 10, 211, 64, 77, 5, 99, 10]
	res = Nmaxelements(my_list, 10)
	print(res)
'''
if __name__ == '__main__':
    main()