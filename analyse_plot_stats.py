import sys
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from Bio import SeqIO
from os import path, stat

#a function to get all file names from the path of a directory
def get_plot_stat_files(path):

        filenames = [f for f in listdir(path) if isfile(join(path, f))]
        return filenames
#there is a separate stat file for each dca energy plot
#this function creates two dictionaries (for train(refined alignment) and test(hmm sequences) data) for each stat file
def get_plot_stats_from_file(filename):

        train_stats = {}
        test_stats = {}
        #checkpoint is a flag to distinguish if the read line is a part of train stats or test stats 
        #look at a stat file from plot_stats folder to understand this better
        checkpoint = '' 

        with open(filename, 'r') as f:
                for line in f:
                        if line.startswith('TRAIN'):
                                checkpoint = 'train'
                        elif line.startswith('TEST'):
                                checkpoint = 'test'
                        elif line.startswith('-'):
                                checkpoint = None
                        line = line.strip('\n')

                        if checkpoint == 'train' and not line.startswith('TRAIN'):
                                line = line.split(':')
                                train_stats[line[0]] = float(line[1])
                        elif checkpoint == 'test' and not line.startswith('TEST'):
                                line = line.split(':')
                                test_stats[line[0]] = float(line[1])

        return train_stats, test_stats

def analyse_stats(train_stats, test_stats):
        #train minus test
        mean_energy_difference, mode_energy_difference, consensus_energy_difference = 0, 0, 0
        train_mean = train_stats['Mean']
        mean_energy_difference = train_mean - test_stats['Mean']
        mode_energy_difference = train_stats['Mode'] - test_stats['Mode']
        consensus_energy_train = train_stats['Consensus energy']
        consensus_energy_difference = consensus_energy_train - test_stats['Consensus energy']

        a1, a2, b1, b2, b3, b4, b5, b6 = 0, 0, 0, 0, 0, 0, 0, 0

        if consensus_energy_difference < 0:
                a1 += 1
               # print('Consensus from refined alignment has lower energy than consensus from hmm emitted sequences')
        else:
                a2 += 1
                #print('Consensus from hmm emitted sequences have lower energy than consensus from refined alignment')

        two_deviations_above_mean_value_train = train_stats['Mean'] + 2*train_stats['SD']
        two_deviations_below_mean_value_train = train_stats['Mean'] - 2*train_stats['SD']
        if consensus_energy_train > train_stats['Max x']:
                b1 += 1
                #print('Consensus energy is out of the distribution (right side - high DCA energy)', '#######################')
        elif consensus_energy_train < train_stats['Min x']:
                b2 += 1
                #print('Consensus energy is out of the distribution (left side - low DCA energy)')
        elif consensus_energy_train > two_deviations_above_mean_value_train:
                b3 += 1
                #print('Consensus energy is two standard deviations above the mean but inside the distribution', '#########################')
        elif consensus_energy_train < two_deviations_below_mean_value_train:
                b4 += 1
                #print('Consensus energy is two standard deviations below the mean but inside the distribution')
        elif consensus_energy_train < two_deviations_above_mean_value_train and consensus_energy_train > train_mean:
                b5 += 1
                #print('Consensus energy is greater than mean but within two standard deviations')
        elif consensus_energy_train > two_deviations_below_mean_value_train and consensus_energy_train < train_mean:
                b6 +=1 
                #print('Consensus energy is below the mean but within two standard deviations')
        #total = a1 + a2 + b1 + b2 + b3 + b4 + b5 + b6
        return [a1, a2, b1, b2, b3, b4, b5, b6]

def scatter_plot(x, y, plot_name, x_label, y_label):
        '''
        if plot_name == 'refined_mode_vs_hmm_mode.png':
            fig, ax = plt.subplots()
            ax.scatter(x, y)
            line = mlines.Line2D([0, 1], [0, 1], color = 'red')
            transform = ax.transAxes
            line.set_transform(transform)
            ax.add_line(line)
        '''

        plot_name = 'cool_plots/' + plot_name
        plt.scatter(x, y)
        if plot_name == 'cool_plots/refined_mode_vs_hmm_mode.png':
            plt.plot(x, x, '-r', label = 'x = y')
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.axis('square')
        plt.savefig(plot_name, bbox_inches='tight')
        plt.clf()
        plt.cla()
        plt.close()

def find_alignment_length_refined(accession):

        refined_alignment_filename = 'refined_alignments/' + accession + '_refined.fasta'
        for record in SeqIO.parse(refined_alignment_filename, 'fasta'):
                loa = len(str(record.seq).upper())
                return loa

def find_number_of_sequences_refined(accession):

        refined_alignment_filename = 'refined_alignments/' + accession + '_refined.fasta'
        l = []
        for record in SeqIO.parse(refined_alignment_filename, 'fasta'):
            l.append(record.seq)
        return len(l)

def find_alignment_length(accession):

        fam_name = 'families/' + accession + '.fasta'
        for record in SeqIO.parse(fam_name, 'fasta'):
            loa = len(str(record.seq).upper())
            return loa

def find_number_of_sequences(accession):
        
        fam_name = 'families/' + accession + '.fasta'
        l = []
        for record in SeqIO.parse(fam_name, 'fasta'):
            l.append(record.seq)
        return len(l)

def main():

        path = 'plot_stats/'
        filenames = get_plot_stat_files(path)
        all_train_stats = []
        refined_consensus_energies = []
        refined_mode_energies = []
        refined_mean_energies = []
        hmm_mode_energies = []
        mode_minus_refined_consensus = []
        mode_minus_hmm_consensus = []
        normed_consensus_energies = []
        normed_mode_energies = []
        #variable for PF00075
        mode_energy_75, consensus_energy_75, mean_energy_75, hmm_mode_75, mode_minus_refined_cs_75, mode_minus_hmm_cs_75 = 0, 0, 0, 0, 0, 0
        l = [0, 0, 0, 0, 0, 0, 0, 0]
        hmm_cse, hmm_modes, refined_cse = [], [], []
        names_excp, nums_excp, loa_excp_refined, nos_excp_refined, loa_excp_original, nos_excp_original = [], [], [], [], [], []
        count = 0
       #for all plot_stats files
        for f in filenames:

                accession = f[0:7]
                loa = find_alignment_length_refined(accession)
                #varible values for PF00075
                if f == 'PF00075_plot_stats.txt':
                        mode_energy_75 = train_stats['Mode']
                        consensus_energy_75 = train_stats['Consensus energy']
                        mean_energy_75 = train_stats['Mean']
                        hmm_mode_75 = test_stats['Mode']
                        mode_minus_refined_cs_75 = abs(mode_energy_75 - consensus_energy_75)
                        mode_minus_hmm_cs_75 = abs(mode_energy_75 - test_stats['Consensus energy'])
                        normed_ce_75 = consensus_energy_75/loa
                        normed_me_75 = mode_energy_75/loa


                fname = path + f
                train_stats, test_stats = get_plot_stats_from_file(fname)
                res = analyse_stats(train_stats, test_stats)
                l = [sum(x) for x in zip(l, res)]
                all_train_stats.append(train_stats)
                refined_consensus_energies.append(train_stats['Consensus energy'])
                refined_mode_energies.append(train_stats['Mode'])
                refined_mean_energies.append(train_stats['Mean'])
                hmm_mode_energies.append(test_stats['Mode'])
                mode_minus_refined_consensus.append(abs(train_stats['Mode'] - train_stats['Consensus energy']))
                mode_minus_hmm_consensus.append(abs(train_stats['Mode'] - test_stats['Consensus energy']))
                #stats for exceptions
                if train_stats['Max x'] < train_stats['Consensus energy']:
                    count += 1
                    names_excp.append(accession)
                    nums_excp.append(count)
                    loa_excp_refined.append(find_alignment_length_refined(accession))
                    nos_excp_refined.append(find_number_of_sequences_refined(accession))
                    loa_excp_original.append(find_alignment_length(accession))
                    nos_excp_original.append(find_number_of_sequences(accession))

                #hmm_cse, hmm_modes, refined_cse = [], [], []
                if train_stats['Consensus energy'] > test_stats['Consensus energy']:
                    hmm_cse.append(test_stats['Consensus energy'])
                    hmm_modes.append(test_stats['Mode'])
                    refined_cse.append(train_stats['Consensus energy'])
                try:
                    normed_consensus_energies.append(train_stats['Consensus energy']/loa)
                    normed_mode_energies.append(train_stats['Mode']/loa)
                except:
                    continue

        #normed consensus energy
        scatter_plot(normed_consensus_energies, normed_mode_energies, 'normed_consensus_energy_vs_mode_energy.png',
                'Consensus Energy/Length of Alignment', 'Mode Energy/Length of Alignment')

        #consensus energy vs mode energy
        #plot_name = 'cool_plots/' + 'normed_consensus_energy_vs_mode_energy.png'
        scatter_plot(refined_consensus_energies, refined_mode_energies, 'consensus_energy_vs_mode_energy.png',
                'Consensus Energy', 'Mode Energy')

        #consensus energy vs mean energy
        scatter_plot(refined_consensus_energies, refined_mean_energies, 'consensus_energy_vs_mean_energy.png',
                'Consensus Energy', 'Mean Energy')

        #refined mode energy vs hmm mode enegy
        scatter_plot(refined_mode_energies, hmm_mode_energies, 'refined_mode_vs_hmm_mode.png',
                'Mode Energy (Refined Alignment)', 'Mode Energy (HMM Alignment)')

        #mode - refined consensus vs mode - hmm consensus
        scatter_plot(mode_minus_refined_consensus, mode_minus_hmm_consensus, 'mode_minus_consensuses.png',
                'Mode - Refined Consensus', 'Mode - HMM Consensus')

        #mode energy and consensus energy correlation
        #wow it's almost 1
        corr = pearsonr(refined_mode_energies, refined_consensus_energies)
        print('Correlation Coefficient (Consensus Energy and Mode Energy): ',corr[0])
        print('l: ',l)
        fig = plt.figure()
        ax = fig.add_axes([0, 0, 1, 1])
        labels = ['Refined MSA CS < HMM CS', 'Refined MSA CS > HMM CS']
        ax.bar(labels, l[0: 2])
        plt.ylabel('Number of sequences')
        plt.savefig('cool_plots/consensus_energies_1.png', bbox_inches = 'tight')
        plt.clf()
        plt.cla()
        plt.close()
        fig = plt.figure()
        labels = ['outside (right - high)', 'outside (left - low)', '> 2 SD inside', '< -2 SD inside', 'b/w mean and 2 SD', 'b/w mean and -2SD']
        ax = fig.add_axes([0, 0, 1, 1])
        ax.bar(labels, l[2: ])
        plt.xlabel('Position in DCA energy plot')
        plt.ylabel('Number of Sequences')
        plt.xticks(rotation = 90)
        plt.savefig('cool_plots/consensus_energies_2.png', bbox_inches = 'tight')
        plt.clf()
        plt.cla()
        plt.close()
        print(hmm_cse, hmm_modes, refined_cse)
        print(len(hmm_cse))
        plt.scatter(hmm_cse, hmm_modes, color = 'blue', label = 'HMM consensus vs HMM mode')
        plt.scatter(refined_cse, hmm_modes, color = 'red', label = 'Refined consensus vs HMM mode')
        plt.legend(loc = 'upper left')
        plt.title('Families with HMM consensus energy < Refined consensus energy')
        #plt.axis('square')
        plt.savefig('cool_plots/cs_vs_hmm_mode(hmm cs < refined cs).png', bbox_inches = 'tight')
        plt.clf()
        plt.cla()
        plt.close()
    
        plt.scatter(nums_excp, loa_excp_refined)
        plt.xlabel('Family')
        plt.ylabel('Length of Refined Alignment')
        plt.savefig('cool_plots/loa_refined_for_deviating_families.png', bbox_inches = 'tight')
        plt.clf()
        plt.cla()
        plt.close()

        plt.scatter(nums_excp, loa_excp_original)
        plt.xlabel('Family')
        plt.ylabel('Length of Original Alignment')
        plt.savefig('cool_plots/loa_original_for_deviating_families.png', bbox_inches = 'tight')
        plt.clf()
        plt.cla()
        plt.close()

        plt.scatter(nums_excp, nos_excp_refined)
        plt.xlabel('Family')
        plt.ylabel('Number of Sequences in Refined Alignment')
        plt.savefig('cool_plots/nos_refined_for_deviating_families.png', bbox_inches = 'tight')
        plt.clf()
        plt.cla()
        plt.close()

        plt.scatter(nums_excp, nos_excp_original)
        plt.xlabel('Family')
        plt.ylabel('Number of Sequences in Orignal Alignment')
        plt.savefig('cool_plots/nos_original_for_deviating_families.png', bbox_inches = 'tight')
        plt.clf()
        plt.cla()
        plt.close()

        print(names_excp, nums_excp, loa_excp_refined, loa_excp_original, nos_excp_refined, nos_excp_original)


if __name__ == '__main__':
        main()
