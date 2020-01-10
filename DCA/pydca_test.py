import sys
import os
sys.path.insert(1, '../')
from consensus import *

mappings = {'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17, 'V': 18, 'W': 19, 'Y': 20, '-': 21}

def mfdca_compute_params(filename):
	cwd = 'mfdca compute_params protein ../refined_alignments/' + filename + '_refined.fasta --verbose'
	os.system(cwd)

def read_couplings_or_fields(filename, couplings_or_fields):
	couplings_filename = 'DCA_output_' + filename + '_refined/' + couplings_or_fields + '_' + filename + '_refined.txt'
	values = []
	with open(couplings_filename, 'r') as fin:
		for line in fin:
			if line.startswith('#'):
				continue
			line = line.strip('\n')
			line = line.split(',')
			values.append(line)

	return values

def main():
	filename = 'PF00131'
	#mfdca_compute_params(filename)
	couplings_or_fields = 'fields'
	values = read_couplings_or_fields(filename, couplings_or_fields)
	#print(values)
	file = '../refined_alignments/' + filename + '_refined.fasta'
	sequences, headers = fasta_to_list(file)

if __name__ == '__main__':
	main()


