import gzip
import os
from Bio import SeqIO, AlignIO
import sys
from collections import Counter, OrderedDict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import copy
from shutil import copyfile
import json
from os import path
import pandas as pd
from ugly_strings import *
import subprocess
sys.path.insert(1, 'DCA/')
from pydca_consensus_copy import *
import threading

start = time.time()

#list of all amino acid letter plus dash
amino_acids = 'ACDEFGHIKLMNPQRSTVWXY'
amino_acids = ['-'] + list(amino_acids)

#finding second largest number in a list(found this on stack overflow)
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

#function for calling clustalo for aligning sequences in fasta format
def fasta_to_clustalo(in_file, out_file):
        cmd = 'clustalo -i ' + in_file + ' -o ' + out_file + ' --force -v'
        os.system(cmd)

#function for calling mafft for aligning sequences in fasta format
def fasta_to_mafft(in_file, out_file):
        cmd = 'mafft ' + in_file + ' > ' + out_file
        os.system(cmd)

#function for calling muscle for aligning sequences in fasta format
def fasta_to_muscle(in_file, out_file):
        cmd = 'muscle -in ' + in_file + ' -out ' + out_file
        os.system(cmd)

#converting fasta file into two lists (sequence list and name/header list)
def fasta_to_list(out_file):
        sequences = []
        name_list = []
        #SeqIO.parse in a function from the biopython module
        for record in SeqIO.parse(out_file, 'fasta'):
                name_list.append(record.id)
                sequences.append(str(record.seq).upper())

        return sequences, name_list

#finding profile matrix for sequences in list format
def profile_matrix(sequences):
        sequence_length = len(sequences[0]) #length of first sequences (length of all sequences is the same after alignment)
        profile_matrix = {} #profile matrix in dictionary format
        for acid in amino_acids:
                profile_matrix[acid] = [float(0) for i in range(sequence_length)] #initialise all entries in profile matrix to zero

        for i in range(len(sequences)):
                seq = sequences[i].upper() #convert sequence to upper case, just in case it isn't
                for j in range(len(seq)): #for each letter in the sequence
                    profile_matrix[seq[j]][j] += float(1) #increase frequency of the letter (seq[j]) at position j

        for aa in profile_matrix: #for amino acid in profile matrix
                l = profile_matrix[aa] #l i sthe list of frequencies associated with that amino acid
                for i in range(len(l)): #for position i in l
                        l[i] /= float(len(sequences)) #divide frequency at i by the length of the list l

        pm = OrderedDict([(x, profile_matrix[x]) for x in amino_acids])

        return pm

#finding index of bad sequence numbers in the sequence list
def find_bad_sequences(profile_matrix, sequences, name_list):
        max_value = max(profile_matrix['-']) #max probability of finding a dash
        if max_value == 1: #just in case there are dashes in every sequence at that position
                max_value = second_largest(profile_matrix['-'])

        positions = [] 
        for i in range(len(profile_matrix['-'])):
                if profile_matrix['-'][i] == max_value: #all positions at which probability of finding a dash is maximum
                        positions.append(i)

        bad_sequence_numbers = []
        for i in range(len(sequences)):
                for position in positions:
                        if sequences[i][position] != '-': #if sequence does not have a dash at the position where probability of finding a dash is the maximum
                                if i not in bad_sequence_numbers: #if sequence is not already in the list
                                        bad_sequence_numbers.append(i)

        return bad_sequence_numbers

#removing bad sequence numbers and returning new sequence list and name/header list
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

#find consensus sequence from sequences in list format (with dashes)
def consensus_sequence(sequences, pm):
        consensus_seq = ''
        #pm = profile_matrix(sequences)
        sequence_length = len(sequences[0]) #length of any sequence
        for i in range(sequence_length):
                l = []
                for aa in pm:
                        l.append(pm[aa][i]) #list of probabilities of amino acid 'aa' at every position
                max_value = max(l) #find maximum value in the above list
                indices = get_all_indices(l, max_value) #get all indices in the list which have the above maximum value
                index = indices[0] #get first index
            
                if amino_acids[index] == '-': #if amino acid at that index is a dash
                        if l[index] < 0.5: #if probability of occurence of dash is less than 0.5
                                second_largest_value = second_largest(l) #then find the second largest value
                                if second_largest_value == max_value: #if second largest and largest and largest values are equal, then get the second index from the list of max values
                                        index = indices[1]
                                else:
                                    index = l.index(second_largest_value) #get index of amino acid with second largest probability of occurence (after dash)
                        #else:
                        #        continue #if probability of occurence of dash is greater than 0.5 then skip adding an amino acid at that position

                consensus_seq += amino_acids[index] #append amino acid to consensus sequence

        return consensus_seq

#find consensus sequence from sequences in list format(without dashes)
def consensus_sequence_nd(sequences, pm):
    consensus_seq = ''
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

#mode of a list
def mode_of_list(sequence_lengths):
        n = len(sequence_lengths)
        data = Counter(sequence_lengths) 
        get_mode = dict(data) 
        mode = [k for k, v in get_mode.items() if v == max(list(data.values()))]
        if n == len(mode):
                return None
        else:
                return mode

#selex format to fasta format conversion (not using this function anywhere as of now)
def selex_to_fasta(in_file, out_file):
        with open(in_file) as fin, open(out_file, 'w') as fout:
                headers = []
                sequences = []
                for line in fin:
                        fout.write('>' + line[0:30].upper() + '\n')
                        fout.write(line[30: ].upper())

#call cd-hit for clustering and removing similar sequences
def cdhit(in_file, out_file):
        cmd = 'cd-hit -i ' + in_file + ' -o ' + out_file + ' -T 1 -c 0.90'
        os.system(cmd)

#get ann indices of a value in a list
def get_all_indices(l, value):

        return [i for i, val in enumerate(l) if val == value]

#stockholm format to fasta format conversion - not using
def stockholm_to_fasta(ifile, ofile):
        with open(ifile, 'r') as fin:
                with open(ofile, 'w') as fout:
                        sequences = SeqIO.parse(ifile, 'stockholm')
                        SeqIO.write(sequences, ofile, 'fasta')
        os.system('rm -rf ' + ifile)
#fasta format to plain format conversion
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
'''
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
'''

#different alignment options
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

#realigning sequences to an existing alignment using mafft
def realign(option, original_alignment, hmm_sequences, out_file):
        if option == '2':
                cwd = 'mafft --add ' + hmm_sequences + ' --reorder --keeplength ' + original_alignment + ' > ' + out_file
                #mafft --add new_sequences --reorder existing_alignment > output
                os.system(cwd)
                #op = subprocess.check_output(cwd, shell=True)
                #with open(out_file, 'w') as fin:
                #       for line in op:
                #               fin.write(line)
        else:
                print('Invalid input')
                sys.exit()

#removing unwanted characters from a filename
def refine_filename(ip):
        ip = str(ip, 'utf-8')
        ip = ip.strip('\n')
        ip = ip.replace('./', '')
        return ip

#removing some accession numbers from a file
def remove_accession(copied_accession_file, accession):
        with open(copied_accession_file, 'r') as f:
                lines = f.readlines()
        with open(copied_accession_file, 'w') as f:
                for line in lines:
                        if line.strip('\n') != accession:
                                f.write(line)

#main function
def main(accession, accession_file):

        print('1. Clustal Omega 2. MAFFT 3. MUSCLE')
        #option = input()
        option = '2' #using only MAFFT for now
        write_file, out_file, temp_file, perc_idens = common_files()

        '''
        #if dca plot exists then exit function
        plot = 'dca_energy_plots/' + accession + '_dca_energies.png'
        if path.exists(plot):
                print('Already calculated')
                return
        '''
        #continue only if family is downloaded
        my_file = 'families/' + accession + '.fasta'
        if path.exists(my_file):

                #0th iteration 
                filename = accession
                print(filename, '%'*30)
                file = 'families/' + filename + '.fasta'
                copyfile(file, temp_file) #copy file to temp_file
                remove_dashes(temp_file, write_file) #remove dashes
                my_dict = {}

                test_seq, test_head = fasta_to_list(write_file)
                if len(test_seq) < 500: #exit if number of sequences is less than 500
                        return

                #cd-hit clustering
                try:
                        cdhit(write_file, out_file)
                except Exception as e:
                        print('Exception: ' + str(e))
                        return

                #getting some strings from ugly_strings.py file
                refined_alignment, plot, final_consensus, profile_hmm, hmm_emitted_sequences, combined_alignment = specific_files(filename)
                plot_file = '/media/Data/consensus/dca_energy_plots/' + filename + '_dca_energies.png'
                #if path.exists(plot_file) and os.stat(plot_file).st_size != 0: #return if dca energy plot already exists
                #    print('plot exists')
                #    return
                ts, th = fasta_to_list(out_file)
                if len(ts) < 500:
                    print('less than 500 sequences after cs-hit clusterting, shall not proceed')
                    return
                #plotting length distribution
                sequence_lengths = sequence_length_list(out_file)
                x = [i for i in range(len(sequence_lengths))]
                plt.scatter(x, sequence_lengths)
                plt.savefig(plot)
                plt.clf()
                plt.cla()
                plt.close()

                #0th alignment step
                mode = mode_of_list(sequence_lengths)[0] #mode of sequence lengths 
                try:
                        alignment(option, out_file, write_file)
                except Exception as e:
                        print('Exception: ' + str(e))
                        return

                iteration = 1
                #exit conditions
                #if number of sequences < 100
                #if length of alignment does not change in subsequent iterations
                #if length of alignment becomes to small i.e -15 the desired length (mode length)
                loa = 0
                condition = 'no_condition'
                #some variables which will be used for the break condition later
                x = 0.1*mode
                y = mode - x
                #try:
                while True:
                    print("Iteration Number: " + str(iteration) + '*'*30)
                    sequences, name_list = fasta_to_list(write_file)
                    number_of_sequences = len(sequences)
                    length_of_alignment = len(sequences[0])
                    print('Length of Alignment = ', length_of_alignment)
                    #here loa is the length of alignment from the previous iteration
                    #length_of_aligment is the length of alignment in the current iteration
                    print('Alignment length (previous iteration): ', loa, 'Alignment length (current iteration): ', length_of_alignment)
                    #saving break condition in a variable
                    if number_of_sequences < 500:
                        condition = 'condition_1'
                    elif length_of_alignment < y:
                        condition = 'condition_2'
                    elif loa == length_of_alignment:
                        condition = 'condition_3'
                    #checking break conditions
                    if number_of_sequences < 500 or length_of_alignment < y or loa == length_of_alignment:
                        f_tag = open('/media/Data/consensus/temp_files/break_tags.txt', 'a')
                        f_tag.write(filename + ' ' + condition) #write break condition along with filename in break_tags.txt 
                        copyfile(write_file, refined_alignment) #write final refined alignment to a file
                        f = open(final_consensus, 'w')
                        f.write('>consensus-from-refined-alignment' + '\n') #write final consensus sequence to a file
                        f.write(cs + '\n')
                        cwd = 'hmmbuild ' + profile_hmm + ' ' + refined_alignment #build profile hmm
                        os.system(cwd)
                        N = number_of_sequences
                        L = length_of_alignment
                        cwd = 'hmmemit -N ' + str(N) + ' -o ' + hmm_emitted_sequences + '-L' + str(L) + ' ' + profile_hmm #emit sequences from prpofile hmm
                        os.system(cwd)
                        os.chdir('/media/Data/consensus/hmm_emitted_sequences')
                        cwd = 'find | grep ' + filename + '_hmmsequences.fasta'#find emitted sequences file
                        ip = subprocess.check_output(cwd, shell=True)
                        ip = refine_filename(ip)
                        ip = 'hmm_emitted_sequences/' + ip
                        os.chdir('/media/Data/consensus')
                        realign(option, refined_alignment, ip, combined_alignment) #align emitted sequences with refined alignment
                        print('***********Final Consensus Sequence from refined alignment: ')
                        print(cs)
                        break
                    #if there is no break condition
                    pm = profile_matrix(sequences) #profile matrix
                    cs = consensus_sequence_nd(sequences, pm) #consensus sequence at current iteration
                    print(cs, len(cs), 'Consensus from refined alignment')
                    bad_sequence_numbers = find_bad_sequences(pm, sequences, name_list) #find bad sequences
                    sequences, name_list = remove_bad_sequences(sequences, name_list, bad_sequence_numbers) #remove bad sequences
                    list_to_fasta(sequences, name_list, temp_file) #convert new list without the bad sequences to fasta format
                    #remove dashes to make file ready for new alignment
                    remove_dashes(temp_file, out_file)
                    alignment(option, out_file, write_file) #realign
                    loa = copy.deepcopy(length_of_alignment) #copy length_of_alignment to loa
                    iteration += 1 #increment interation number

                #except Exception as e:
                #        print(str(e))
                #        time.sleep(5)
                #        return
                
                end = time.time() - start
                print('It took ' + str(end) + ' seconds to run the iterative alignment for: ' + filename)
                print('\n\n\n')
                #print('Time for next protein family/domain/motif\n')
                time.sleep(2)

if __name__ == '__main__':
        accession_file = 'temp_files/accession_list.txt'
        #copied_accession_file = 'temp_files/accession_list_copy.txt'
        #original_accession_file = 'temp_files/exceptions.txt'
        #copied_accession_file = 'temp_files/exceptions_copy.txt'
        #copyfile(original_accession_file, copied_accession_file)
        accession_list = []
        with open(accession_file, 'r') as f:
                for line in f:
                        accession_list.append(line.strip('\n')) #make a list of accession numbers from accession file

        #calling main function for all families in accession_list
        accession_list = ['PF00167', 'PF12079', 'PF05065', 'PF04398', 'PF00902']
        #accession_list = ['PF00902']
        for i in range(4, 7):
            for accession in accession_list: #for each family
                print('Iterative Alignment ', accession)
                time.sleep(2)
                t1 = threading.Thread(target = main, args = (accession, accession_file, )) #do iterative alignment
                t1.setDaemon(True)
                t1.start()
                t1.join()
                os.chdir('/media/Data/consensus/DCA') #change directory to call DCA script
                print('DCA calculation ', accession)
                time.sleep(2)
                t2 = threading.Thread(target = main_pydca, args = (accession, accession_file, )) #call dca script
                t2.setDaemon(True)
                t2.start()
                t2.join()
                os.chdir('/media/Data/consensus')
                print(accession, ' DONE!!!!')
                os.system('cp dca_energy_plots/' + accession + '_dca_energies.png ' + accession + '/' + accession + '_dca_energies_' + str(i) + '.png' )
                os.system('cp hmm_consensuses/' + accession + '_hmm_consensus.fasta ' + accession + '/'+ accession + '_hmm_consensus_' + str(i) + '.fasta')
                os.system('cp plot_stats/' + accession + '_plot_stats.txt ' + accession + '/' + accession + '_plot_stats_' + str(i) + '.txt')
                os.system('cp refined_alignments/' + accession + '_refined.fasta ' + accession + '/' + accession + '_refined_' + str(i) + '.fasta')
                os.system('cp refined_consensuses/' + accession + '_refined_consensus.fasta ' + accession + '/' + accession + '_refined_consensus_' + str(i) + '.fasta')
                os.system('cp all_consensus_sequences/' + accession + '_consensus.fasta ' + accession + '/' + accession + '_consensus_' + str(i) + '.fasta')
