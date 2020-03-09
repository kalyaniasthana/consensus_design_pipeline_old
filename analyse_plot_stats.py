import sys
from os import listdir
from os.path import isfile, join


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


def main():

	path = 'plot_stats/'
	filenames = get_plot_stat_files(path)
	for f in filenames:
		print(f, '************************')
		fname = path + f
		train_stats, test_stats = get_plot_stats_from_file(fname)
		print(analyse_stats(train_stats, test_stats))

if __name__ == '__main__':
	main()