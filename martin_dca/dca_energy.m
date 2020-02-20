fasta_test = '/media/Data/consensus/temp_files/only_hmm_emitted.fasta';
fasta_train = '/media/Data/consensus/temp_files/train_file.fasta';
accession = 'PF02793';
[score_train, score_test] = calculate_dca_scores(fasta_train,fasta_test,accession)