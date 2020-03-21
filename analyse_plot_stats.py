import sys
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from Bio import SeqIO


def get_plot_stat_files(path):

	filenames = [f for f in listdir(path) if isfile(join(path, f))]
	return filenames

def get_plot_stats_from_file(filename):

	train_stats = {}
	test_stats = {}
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

	if consensus_energy_difference < 0:
		print('Consensus from refined alignment has lower energy than consensus from hmm emitted sequences')
	else:
		print('Consensus from hmm emitted sequences have lower energy than consensus from refined alignment')

	two_deviations_above_mean_value_train = train_stats['Mean'] + 2*train_stats['SD']
	two_deviations_below_mean_value_train = train_stats['Mean'] - 2*train_stats['SD']
	if consensus_energy_train > train_stats['Max x']:
		print('Consensus energy is out of the distribution (right side - high DCA energy)', '#######################')
	elif consensus_energy_train < train_stats['Min x']:
		print('Consensus energy is out of the distribution (left side - low DCA energy)')
	elif consensus_energy_train > two_deviations_above_mean_value_train:
		print('Consensus energy is two standard deviations above the mean but inside the distribution', '#########################')
	elif consensus_energy_train < two_deviations_below_mean_value_train:
		print('Consensus energy is two standard deviations below the mean but inside the distribution')
	elif consensus_energy_train < two_deviations_above_mean_value_train and consensus_energy_train > train_mean:
		print('Consensus energy is greater than mean but within two standard deviations')
	elif consensus_energy_train > two_deviations_below_mean_value_train and consensus_energy_train < train_mean:
		print('Consensus energy is below the mean but within two standard deviations')

def scatter_plot(x, y, plot_name, x_label, y_label):

	plot_name = 'cool_plots/' + plot_name
	plt.scatter(x, y)
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.savefig(plot_name, bbox_inches='tight')
	plt.clf()
	plt.cla()
	plt.close()

def find_alignment_length(accession):

	refined_alignment_filename = 'refined_alignments/' + accession + '_refined.fasta'
	for record in SeqIO.parse(refined_alignment_filename, 'fasta'):
		loa = len(str(record.seq).upper())
		return loa

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
	mode_energy_75, consensus_energy_75, mean_energy_75, hmm_mode_75, mode_minus_refined_cs_75, mode_minus_hmm_cs_75 = 0, 0, 0, 0, 0, 0

	for f in filenames:
		#print(f, '************************')
		#sys.exit()
		accession = f[0:7]
		loa = find_alignment_length(accession)

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
		all_train_stats.append(train_stats)
		refined_consensus_energies.append(train_stats['Consensus energy'])
		refined_mode_energies.append(train_stats['Mode'])
		refined_mean_energies.append(train_stats['Mean'])
		hmm_mode_energies.append(test_stats['Mode'])
		mode_minus_refined_consensus.append(abs(train_stats['Mode'] - train_stats['Consensus energy']))
		mode_minus_hmm_consensus.append(abs(train_stats['Mode'] - test_stats['Consensus energy']))
		normed_consensus_energies.append(train_stats['Consensus energy']/loa)
		normed_mode_energies.append(train_stats['Mode']/loa)

	plot_name = 'cool_plots/' + 'consensus_energy_vs_mode_energy'
	plt.scatter(refined_consensus_energies, refined_mode_energies, color = 'blue')
	plt.scatter([consensus_energy_75], [mode_energy_75], color = 'red', label = 'PF00075')
	plt.legend(loc = 'upper left')
	plt.xlabel('Consensus Energy')
	plt.ylabel('Mode Energy')
	plt.savefig(plot_name, bbox_inches='tight')
	plt.clf()
	plt.cla()
	plt.close()

	plot_name = 'cool_plots/' + 'normed_consensus_energy_vs_mode_energy.png'
	plt.scatter(normed_consensus_energies, normed_mode_energies, color = 'blue')
	#plt.scatter([normed_ce_75], [normed_me_75], color = 'red', label = 'PF00075')
	#plt.legend(loc = 'upper left')
	plt.xlabel('Consensus Energy/Length of Alignment')
	plt.ylabel('Mode Energy/Length of Alignment')
	plt.savefig(plot_name, bbox_inches='tight')
	plt.clf()
	plt.cla()
	plt.close()


	plot_name = 'cool_plots/' + 'consensus_energy_vs_mean_energy'
	plt.scatter(refined_consensus_energies, refined_mean_energies, color = 'blue')
	plt.scatter([consensus_energy_75], [mean_energy_75], color = 'red', label = 'PF00075')
	plt.legend(loc = 'upper left')
	plt.xlabel('Consensus Energy')
	plt.ylabel('Mean Energy')
	plt.savefig(plot_name, bbox_inches='tight')
	plt.clf()
	plt.cla()
	plt.close()

	plot_name = 'cool_plots/' + 'refined_mode_vs_hmm_mode.png'
	plt.scatter(refined_mode_energies, hmm_mode_energies, color = 'blue')
	plt.scatter([mode_energy_75], [hmm_mode_75], color = 'red', label = 'PF00075')
	plt.legend(loc = 'upper left')
	plt.xlabel('Mode energies (refined alignment)')
	plt.ylabel('Mode energies (HMM alignment)')
	plt.savefig(plot_name, bbox_inches='tight')
	plt.clf()
	plt.cla()
	plt.close()

	plot_name = 'cool_plots/' + 'mode_minus_consensuses.png'
	plt.scatter(mode_minus_refined_consensus, mode_minus_hmm_consensus, color = 'blue')
	plt.scatter([mode_minus_refined_cs_75], [mode_minus_hmm_cs_75], color = 'red', label = 'PF00075')
	plt.legend(loc = 'upper left')
	plt.xlabel('(Mode Energy) - (Refined Consensus Energy)')
	plt.ylabel('(Mode Energy) - (HMM Consensus Energy)')
	plt.savefig(plot_name, bbox_inches='tight')
	plt.clf()
	plt.cla()
	plt.close()

	'''
	#consensus energy vs mode energy
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

	'''
	#mode energy and consensus energy correlation
	#wow it's almost 1
	corr = pearsonr(refined_mode_energies, refined_consensus_energies)
	print(corr[0])


if __name__ == '__main__':
	main()