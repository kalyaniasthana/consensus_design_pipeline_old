function [ C, Pi, Pi_pcred, N, q, align, faa ] = correls(nomefile, theta, pseudocount_weight)

[Pij_true,Pi_true,N,q,align,faa] = read_alignment(nomefile, theta);
[Pij,Pi] = with_pc(Pij_true,Pi_true,pseudocount_weight,N,q);

C = Compute_C(Pij,Pi,N,q);  % covariance matrix

[~,Pi_pcred] = with_pc(Pij_true,Pi_true,0.001,N,q);

end

function [Pij_true,Pi_true,N,q,align,faa] = read_alignment(nomefile, theta)

[N,M,align,faa] = read_alignment_fasta(nomefile); 

% reweighting (only if M<20000 since pdist becomes large)

W = ones(1,M);

if( theta > 0.0 )
   if ( M<200 )
       W = (1./(1+sum(squareform(pdist(align,'hamm')<theta))));
   else
       align_vec = zeros(1,N*M);
       for i=1:M
           align_vec( 1+(i-1)*N : i*N ) = align( i, : );
       end
       tic
       W = weightCalculator(align_vec, theta, M, N);
       toc
       clear align_vec
   end
end
Meff=sum(W);

q = max(max(align));

[N,M,Meff,q]

% frequencies of amino acid occurrence

Pij_true = zeros(N,N,q,q);
Pi_true = zeros(N,q);

for j=1:M
    for i=1:N
        Pi_true(i,align(j,i)) = Pi_true(i,align(j,i)) + W(j);
    end
end
Pi_true = Pi_true/Meff;

for l=1:M
    for i=1:N-1
        for j=i+1:N
            Pij_true(i,j,align(l,i),align(l,j)) = Pij_true(i,j,align(l,i),align(l,j)) + W(l);
            Pij_true(j,i,align(l,j),align(l,i)) = Pij_true(i,j,align(l,i),align(l,j));
        end
    end
end
Pij_true = Pij_true/Meff;

scra = eye(q,q);
for i=1:N
    for alpha=1:q
        for beta=1:q
            Pij_true(i,i,alpha,beta) = Pi_true(i,alpha) * scra(alpha,beta);
        end
    end
end
end

function [N,M,Z,X] = read_alignment_fasta(nomefile)
X = fastaread(nomefile);
N = size(X(1).Sequence,2);
M = size(X,1);
Z = zeros(M,N);
for i=1:M
    for j=1:N
        Z(i,j)=letter2number(X(i).Sequence(j));
    end
end
end

function x=letter2number(a)
switch(a)

    case '-'
         x=1;
    case 'A'    
        x=2;    
    case 'C'    
        x=3;
    case 'D'
        x=4;
    case 'E'  
        x=5;
    case 'F'
        x=6;
    case 'G'  
        x=7;
    case 'H'
        x=8;
    case 'I'  
        x=9;
    case 'K'
        x=10;
    case 'L'  
        x=11;
    case 'M'
        x=12;
    case 'N'  
        x=13;
    case 'P'
        x=14;
    case 'Q'
        x=15;
    case 'R'
        x=16;
    case 'S'  
        x=17;
    case 'T'
        x=18;
    case 'V'
        x=19;
    case 'W'
        x=20;
    case 'Y'
        x=21;
    otherwise
        x=1;
end
end

function [Pij,Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight,N,q)

Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*ones(N,N,q,q);
Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*ones(N,q);

scra = eye(q);

for i=1:N
    for alpha = 1:q
        for beta = 1:q
           Pij(i,i,alpha,beta) =  (1.-pseudocount_weight)*Pij_true(i,i,alpha,beta) + pseudocount_weight/q*scra(alpha,beta);
        end
    end
end 

end

function C = Compute_C(Pij,Pi,N,q)

C=zeros(N*(q-1),N*(q-1));

for i=1:N
    for j=1:N
        for alpha=1:q-1
            for beta=1:q-1
                 C(mapkey(i,alpha,q),mapkey(j,beta,q)) = Pij(i,j,alpha,beta) - Pi(i,alpha)*Pi(j,beta);
            end
        end
    end
end
end

function A=mapkey(i,alpha,q)
A = (q-1)*(i-1)+alpha;
end
