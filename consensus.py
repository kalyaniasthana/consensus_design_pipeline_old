import gzip
import os
from Bio import SeqIO

amino_acids = 'ACDEFG'

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

def profile_matrix(sequences):
	sequence_length = len(sequences[0])
	profile_matrix = {}


def main():
    gzfile = 'PF00167_full_length_sequences.fasta.gz'
    write_file = 'write.fasta'
    out_file = 'output.fasta'
    gzip_to_fasta(gzfile, write_file)
    fasta_to_clustalo(write_file, out_file)
    sequences, name_list = clustalo_output_to_list(out_file)
    print(sequences[0], name_list[0:5])

if __name__ == '__main__':
    main()