function [eij,hi] = coupling_b( C, Pi, N, q )

invC = inv(C);
%file = fopen('couplings.txt','w');
eij = zeros(N,N,q,q);
for i=1:N
   for j=1:N
      for a = 1:(q-1)
         for b = 1:(q-1)
            eij(i,j,a,b) = 0;                   
            eij(i,j,a,b) = - invC(mapkey(i,a,q),mapkey(j,b,q)) ;   % comment for profile model
            %fprintf(file, '%d, %d, %d, %d, %d\n', i, j, a, b, eij(i,j,a,b));
          end
      end
   end
end
%fclose(file);
%file = fopen('fields.txt', 'w');
hi = zeros(N,q);
for i = 1:N
    for a = 1:(q-1)
        hi(i,a) = log( Pi(i,a)/Pi(i,q) );          % standard mean field
        for j = 1:N
            for b = 1:(q-1)
              hi(i,a) = hi(i,a) - eij(i,j,a,b)*Pi(j,b);
            end
        %fprintf(file, '%d, %d, %d\n', i, a, hi(i,a));
        end
    end
end
%fclose(file);
end

function A=mapkey(i,alpha,q)
A = (q-1)*(i-1)+alpha;
end
