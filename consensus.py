import gzip
import os
from Bio import SeqIO

amino_acids = 'ACDEFGHIKLMNPQRSTVWXY'
amino_acids = list(amino_acids) + ['-']

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

def gzip_to_fasta(gzfile, wf):
	write_file = open(wf, 'w')
	with gzip.open(gzfile, 'rb') as f:
		for line in f:
			line = line.decode('utf-8')
			write_file.write(line)

def fasta_to_clustalo(in_file, out_file):
	cmd = 'clustalo -i ' + in_file + ' -o ' + out_file + ' --force -v'
	os.system(cmd)

def clustalo_output_to_list(out_file):
	sequences = []
	name_list = []

	for record in SeqIO.parse(out_file, 'fasta'):
		name_list.append(record.id)
		sequences.append(str(record.seq))

	return sequences, name_list

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

	#bad_sequences = []
	#for number in bad_sequence_numbers:
	#		bad_sequences.append(name_list[number])

	return bad_sequence_numbers

def remove_bad_sequences(sequences, name_list, bad_sequence_numbers):

	bad_sequence_set = set(bad_sequence_numbers)
	sequences = [x for i, x in enumerate(sequences) if i not in bad_sequence_set]
	name_list = [x for i, x in enumerate(name_list) if i not in bad_sequence_set]

	return sequences, name_list

def main():
	cut_off = 500
	gzfile = 'PF00167_full_length_sequences.fasta.gz'
	write_file = 'write.fasta'
	out_file = 'output.fasta'
	#gzip_to_fasta(gzfile, write_file)
	#fasta_to_clustalo(write_file, out_file)
	sequences, name_list = clustalo_output_to_list(out_file)
	pm = profile_matrix(sequences, 1)
	bad_sequence_numbers = find_bad_sequences(pm, sequences, name_list)

	print(remove_bad_sequences(sequences, name_list, bad_sequences, bad_sequence_numbers))
if __name__ == '__main__':
    main()