#each of these functions defines some strings, which are all paths to a file or a directory
#this is to make the larger scripts look cleaner

write_file, out_file, temp_file, perc_idens, refined_alignment, plot, final_consensus, profile_hmm, emitted_alignment, emitted_alignment_fasta, aligned_cs = '', '', '', '', '', '', '', '', '', '', ''

def common_files():
    write_file = 'temp_files/write.fasta'
    out_file = 'temp_files/output.fasta'
    temp_file = 'temp_files/temp.fasta'
    perc_idens = 'temp_files/percentage_identities_file'
    return write_file, out_file, temp_file, perc_idens

def specific_files(filename):
    refined_alignment = 'refined_alignments/' + filename + '_refined.fasta'
    plot = 'length_distributions/' + filename + '_length_distribution.png'
    final_consensus = 'refined_consensuses/' + filename + '_refined_consensus.fasta'
    profile_hmm = 'hmm_profiles/' + filename + '_profile.hmm'
    hmm_emitted_sequences = 'hmm_emitted_sequences/' + filename + '_hmmsequences.fasta'
    combined_alignment = 'combined_alignments/' + filename + '_combined.fasta'
    return refined_alignment, plot, final_consensus, profile_hmm, hmm_emitted_sequences, combined_alignment

def pydca_strings(filename):
    combined_file = '../combined_alignments/' + filename + '_combined.fasta'
    train_file = '../temp_files/train_file.fasta'
    test_file = '../temp_files/test_file.fasta'
    #only_refined = '../temp_files/only_refined_sequences.fasta'
    only_refined = '../refined_alignments/' + filename + '_refined.fasta'
    only_hmm = '../temp_files/only_hmm_emitted.fasta'
    dca_energy_plot = '../dca_energy_plots/' + filename + '_dca_energies.png'
    consensus_file = '../all_consensus_sequences/' + filename + '_consensus.fasta'
    combined_with_consensus = '../temp_files/combined_with_consensus.fasta'
    refined_consensus = '../refined_consensuses/' + filename + '_refined_consensus.fasta'
    hmm_consensus = '../hmm_consensuses/' + filename + '_hmm_consensus.fasta'
    return combined_file, train_file, test_file, only_refined, only_hmm, dca_energy_plot, refined_consensus, hmm_consensus, consensus_file, combined_with_consensus
