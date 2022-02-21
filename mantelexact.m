function [manstat, p1, index, np, r] = mantelexact(mata, matb)
%MANTELEXACT performs an EXACT permuation test of the Mantel
%             statistic for two n * n matrices 'mata' and 'matb'
% mata and matb are a pair of n*n matrices
% manstat is the Mantel statistic for the random permutation
% p1 is a one-tailed p-value
% index is the number permutations with a statistic >= manstat
% np is the number of permutations evaulated
% r is the correlation between the off-diag elements
% NOTES:
%  1. The main diagonal is zeroed out and ignored

tic;
n=size(mata, 1);
ic=1;
mata = mata.*~eye(n);
matb = matb.*~eye(n);
manstat=sum(sum(mata .* matb));
x = zeros(n.*(n-1),1);
y = zeros(n.*(n-1),1);
ict = 0;
for i = 1:n
    for j = 1:n
        if i == j
            continue
        end
        ict = ict + 1;
        x(ict) = mata(i,j);
        y(ict) = matb(i,j);
    end
end
r = corrcoef(x,y);

index = 1;              % Number of perms. more extreme than Rho (one-tail)
np=1;                   % Total perms
ipivot=n-1;
a = 1:n;
c = ones(1,n);

while ipivot > 0        % SYSTEMATIC generation of all permutations

  if a(ipivot) > a(ipivot+1)     % If object on left > right
    ipivot = ipivot-1;           % then reduce the marker by 1.
    continue
  else
    mark=a(ipivot);              % Otherwise, free up the allocation
    for j = ipivot:n             % for the marker and everything right
      c(a(j))=0;
    end
    for j = mark+1:n             % Set the object in the 'ipivot'
      if c(j) == 0               % spot to its next largest avail.
        a(ipivot)=j;             % value
        c(j)=1;
        break
      end
    end                          % Fill-up everything to the right
    mark2=ipivot;                % in order
    for j = 1:n
      if c(j) == 0
        mark2=mark2+1;
        a(mark2)=j;
        c(j)=1;
      end
    end
                                 % Increment permutation counter
    np=np+1;                     % and evaluate
    manval = sum(sum(mata.*matb(a,a)));
    if manval >= manstat 
      index = index + 1;
    end
    ipivot=n-1;
  end
end
p1 = index ./ np;                % One-tailed p-value
toc
