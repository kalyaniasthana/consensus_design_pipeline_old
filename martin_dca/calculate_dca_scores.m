function [score_train, score_test] = calculate_dca_scores(fasta_train,fasta_test,accession)

% code calculates DCA energies for all training and test sequences
% (consensus or hmmemit seqs to be  presented as test sequences)
% before running, compile weightCalculator.c with mex inside Matlab

[ C, ~, Pi_pcred, N, q, ~, ~ ] = correls(fasta_train, 0.2, 0.5); % standard reweighting and pseudocount, can be changed
[eij,hi] = coupling_b( C, Pi_pcred, N, q );
score_train = score_fct( fasta_train , eij, hi, N, q);
score_test = score_fct( fasta_test , eij, hi, N, q);
cf = strcat('/media/Data/consensus/all_consensus_sequences/', accession);
consensus_file = strcat(cf, '_consensus.fasta');
score_consensus = score_fct( consensus_file, eij, hi, N, q);
fprintf('%d\n', score_consensus(1));
fprintf('%d', score_consensus(2));
figure;
histogram(-score_train, 'Normalization', 'prob', 'BinWidth', 20);
hold;
xlimtrain = get(gca, 'xlim');
file = fopen('/media/Data/consensus/temp_files/exceptions.txt', 'at');
histogram(-score_test, 'Normalization', 'prob', 'BinWidth', 20);
hold;
y1 = get(gca, 'ylim');
line([-score_consensus(1) -score_consensus(1)], y1, 'Color', 'g');
hold;
line([-score_consensus(2) -score_consensus(2)], y1, 'Color', 'r');
s = strcat('DCA energy for: ', accession);
xlabel(s);
legend('Refined MSA', 'HMM emitted MSA', 'Consensus from Refined MSA', 'Consensus from hmm emitted MSA');
plt = strcat('../dca_energy_plots/', accession);
plot_name = strcat(plt, '_dca_energies');
print(plot_name, '-dpng');
fprintf('\n\n%d\t', xlimtrain(2));
fprintf('%d\n\n', -score_consensus(1));
fprintf('%d\n', xlimtrain(2) > -score_consensus(1));
if xlimtrain(2) < -score_consensus(1)
	fprintf(file, '%s\n', accession);
fclose(file);
clear eij;
clear hi;
end