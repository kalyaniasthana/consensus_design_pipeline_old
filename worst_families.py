from Bio import SeqIO
from analyse_plot_stats import *
from consensus import *
from tabulate import tabulate

def count_number_of_gaps(sequence):
    gaps = 0
    for char in sequence:
        if char == '-':
            gaps += 1
    return gaps

def main():

    path = 'worst_plots/'
    wf = get_plot_stat_files(path)
    #accessions_wf = [i[:7] for i in wf]
    wf = [i[:7] + '_plot_stats.txt' for i in wf]
    tuples_for_table = []

    for filename in wf:

        f = 'plot_stats/' + filename
        accession = filename[:7]

        train_stats, test_stats = get_plot_stats_from_file(f)
        ref_mean_energy = train_stats['Mean']
        ref_mode_energy = train_stats['Mode']
        ref_cs_energy = train_stats['Consensus energy']

        hmm_mean_energy = test_stats['Mean']
        hmm_mode_energy = test_stats['Mode']
        hmm_cs_energy = test_stats['Consensus energy']

        ref_alignment_file = 'refined_alignments/' + accession + '_refined.fasta'
        consensus_file = 'all_consensus_sequences/' + accession + '_consensus.fasta'
        
        ref_alignment_sequences, ignore = fasta_to_list(ref_alignment_file)
        no_of_seqs = len(ref_alignment_sequences)
        loa = len(ref_alignment_sequences[0])
    
        consensus_sequences, ignore = fasta_to_list(consensus_file)
        refined_cs = consensus_sequences[0]
        hmm_cs = consensus_sequences[1]

        refined_cs_gaps = count_number_of_gaps(refined_cs)
        hmm_cs_gaps = count_number_of_gaps(hmm_cs)
        
        mean_diff = ref_mean_energy - hmm_mean_energy
        mode_diff = ref_mode_energy - hmm_mode_energy
        cs_diff = ref_cs_energy - hmm_cs_energy

        tuples_for_table.append((accession, ref_mean_energy, ref_mode_energy, ref_cs_energy, hmm_mean_energy,
            hmm_mode_energy, hmm_cs_energy, mean_diff, mode_diff, cs_diff, no_of_seqs, loa, refined_cs_gaps, hmm_cs_gaps))

    headers = ['Family', 'r_mean_e', 'r_mode_en', 'r_cs_e', 'h_mean_e', 'h_mode_e', 'h_cs_e'
            , 'mean_diff', 'mode_diff', 'cs_diff', 'nos', 
            'loa', 'ngaps_r_cs', 'ngaps_h_cs']

    print(tabulate(tuples_for_table, headers, tablefmt = "grid"))
    print('Here: ')
    print('r_mean_e: mean energy of refined alignment')
    print('r_mode_e: mode energy of refined alignment')
    print('r_cs_e: energy of consensus from refined alignment')
    print('h_mean_e: mean energy of hmm sequences')
    print('h_mode_e: mode energy of hmm sequences')
    print('h_cs_e: energy of consensus from hmm sequences')
    print('mean_diff: r_mean_e - h_mean_e')
    print('mode_diff: r_mode_e - h_mode_e')
    print('cs_diff: r_cs_e - h_cs_e')
    print('nos: number of sequences in refined/hmm alignment')
    print('loa: length of refined/hmm alignment')
    print('ngaps_r_cs: number of gaps in consensus from refined alignment')
    print('ngaps_h_cs: number of gaps in consensus from hmm sequences')

if __name__ == '__main__':
    main()
