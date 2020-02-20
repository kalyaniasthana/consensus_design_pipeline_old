function [ score ] = score_fct( nomefile , eij, hi, N, q)

[N1,M,Z] = read_alignment_fasta(nomefile);

%[N1,M]

score = zeros(M,1);

for m = 1:M
   score(m) = - energy(m, Z, eij, hi, N, q);  
end

% score = score - mean(score);

end

function [N1,M,Z] = read_alignment_fasta(nomefile)
X = fastaread(nomefile);
N1 = size(X(1).Sequence,2);
M = size(X,1);
Z = zeros(M,N1);
for i=1:M
    for j=1:N1
        Z(i,j)=letter2number(X(i).Sequence(j));
    end
end
end

function x=letter2number(a)
switch(a)

    % full AA alphabet

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

function e = energy(m, Z, eij, hi, N, q) 

e = 0;

for i = 1:N
        e = e - hi( i , Z(m,i) );
        for j = 1:N
            %if ( i ~= j )
                e = e - eij(i,j,Z(m,i),Z(m,j))/2;
            %end
        end
end
            
end