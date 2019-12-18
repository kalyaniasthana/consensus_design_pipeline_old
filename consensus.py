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

def fasta_to_list(out_file):
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

def list_to_fasta(sequences, name_list, fasta_file):
	file = open(fasta_file, 'w')
	for i in range(len(sequences)):
		file.write('>' + name_list[i] + '\n' + sequences[i] + '\n')

	file.close()

def remove_dashes(fasta_file_from, fasta_file_to):
	with open(fasta_file_from) as fin, open(fasta_file_to, 'w') as fout:
		for line in fin:
			fout.write(line.translate(str.maketrans('', '', '-')))

def consensus_sequence(sequences):
	consensus_seq = ''
	pm = profile_matrix(sequences, 1)
	sequence_length = len(sequences[0])
	for i in range(sequence_length):
		l = []
		for aa in pm:
			l.append(pm[aa][i])
		index = l.index(max(l))
		consensus_seq += aa[index]

	return consensus_seq


def main():
	'''
	cut_off = 2800
	gzfile = 'PF00167_full_length_sequences.fasta.gz'
	write_file = 'write.fasta'
	out_file = 'output.fasta'
	gzip_to_fasta(gzfile, write_file)
	fasta_to_clustalo(write_file, out_file)
	sequences, name_list = fasta_to_list(out_file)
	#print(len(sequences), '###########################')
	pm = profile_matrix(sequences, 1)
	bad_sequence_numbers = find_bad_sequences(pm, sequences, name_list)
	sequences, name_list = remove_bad_sequences(sequences, name_list, bad_sequence_numbers)
	list_to_fasta(sequences, name_list, 'test.fasta')
	remove_dashes('test.fasta', write_file)
	fasta_to_clustalo(write_file, out_file)
	sequences, name_list = fasta_to_list(out_file)
	print(len(sequences), '###########################')
	pm = profile_matrix(sequences, 1)
	bad_sequence_numbers = find_bad_sequences(pm, sequences, name_list)
	sequences, name_list = remove_bad_sequences(sequences, name_list, bad_sequence_numbers)
	print(len(sequences), '###########################')

	'''
	cut_off = 3200
	gzfile = 'PF00167_full_length_sequences.fasta.gz'
	write_file = 'write.fasta'
	out_file = 'output.fasta'
	gzip_to_fasta(gzfile, write_file)

	seq, names = fasta_to_list(write_file)
	num = len(seq)
	iteration = 1
	pm = {}


	while num > cut_off:
		print("ITERATION NUMBER " + str(iteration) + '#'*100)
		fasta_to_clustalo(write_file, out_file)
		sequences, name_list = fasta_to_list(out_file)

		pm = profile_matrix(sequences, 1)
		bad_sequence_numbers = find_bad_sequences(pm, sequences, name_list)
		sequences, name_list = remove_bad_sequences(sequences, name_list, bad_sequence_numbers)
		list_to_fasta(sequences, name_list, 'temp.fasta')
		remove_dashes('temp.fasta', write_file)
		num = len(sequences)
		iteration += 1

	print(consensus_sequence(sequences))


if __name__ == '__main__':
    main()