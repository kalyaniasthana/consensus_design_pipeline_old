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

%fprintf('%d\n', score_consensus(1));
%fprintf('%d\n', score_consensus(2));

figure;
histogram(-score_train, 'Normalization', 'prob', 'BinWidth', 20);
hold;
xlimtrain = get(gca, 'xlim');

fname = '/media/Data/consensus/plot_stats/';
fname_open = strcat(fname, accession);
fname = strcat(fname_open, '_plot_stats.txt');
file = fopen(fname, 'w');

histogram(-score_test, 'Normalization', 'prob', 'BinWidth', 20);
xlimtest = get(gca, 'xlim');
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

% train stats
mean_train = -mean(score_train);
mode_train = -mode(score_train);
sd_train = std(score_train);
consensus_train = -score_consensus(1);
min_x_train = -max(score_train);
max_x_train = -min(score_train);

fprintf(file, '%s\n', 'TRAIN STATS (Refined Alignment)');
fprintf(file, '%s\t%.6f\n', 'Mean: ', mean_train);
fprintf(file, '%s\t%.6f\n', 'Mode: ', mode_train);
fprintf(file, '%s\t%.6f\n', 'SD: ', sd_train);
fprintf(file, '%s\t%.6f\n', 'Consensus energy: ', consensus_train);
fprintf(file, '%s\t%.6f\n', 'Min x: ', min_x_train);
fprintf(file, '%s\t%.6f\n', 'Max x: ', max_x_train);

% test stats
mean_test = -mean(score_test);
mode_test = -mode(score_test);
sd_test = std(score_test);
consensus_test = -score_consensus(2);
min_x_test = -max(score_test);
max_x_test = -min(score_test);

fprintf(file, '%s\n', '-------------------------------------------------------');
fprintf(file, '%s\n', 'TEST STATS (HMM emitted sequences)');
fprintf(file, '%s\t%.6f\n', 'Mean: ', mean_test);
fprintf(file, '%s\t%.6f\n', 'Mode: ', mode_test);
fprintf(file, '%s\t%.6f\n', 'SD: ', sd_test);
fprintf(file, '%s\t%.6f\n', 'Consensus energy: ', consensus_test);
fprintf(file, '%s\t%.6f\n', 'Min x: ', min_x_test);
fprintf(file, '%s\t%.6f\n', 'Max x: ', max_x_test);

fclose(file);

clear eij;
clear hi;
clear plot_name;
clear plt;
clear s;
clear file;
clear xlimtrain;
clear score_consensus;
clear consensus_file;
clear cf;
clear score_train;
clear score_test;
clear C;
clear Pi_pcred;
clear N;
clear q;
clear all;
close all;
end