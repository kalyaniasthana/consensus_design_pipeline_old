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
	final_consensus = 'all_consensus_sequences/' + filename + '_consensus.fasta'
	profile_hmm = 'hmm_profiles/' + filename + '_profile.hmm'
	hmm_emitted_sequences = 'hmm_emitted_sequences/' + filename + '_hmmsequences.fasta'
	combined_alignment = 'combined_alignments/' + filename + '_combined.fasta'
	return refined_alignment, plot, final_consensus, profile_hmm, hmm_emitted_sequences, combined_alignment
