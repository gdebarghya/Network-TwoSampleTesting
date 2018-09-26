function [test,pval] = NormalityTest(A,B,sig)

% returns acceptance/rejection for proposed Frobenius norm based test
% assumes A, B are defined as a cell array of length m, each cell being an
% adjacency matrix (stored as a sparse matrix)

m1 = floor(0.5*min(length(A),length(B)));
n = size(A{1},1);

SA1 = sparse(n,n); SA2 = sparse(n,n); 
SB1 = sparse(n,n); SB2 = sparse(n,n);

for i = 1:m1
    SA1 = SA1 + A{i}; SA2 = SA2 + A{i+m1};
    SB1 = SB1 + B{i}; SB2 = SB2 + B{i+m1};
end
nummat = triu(SA1-SB1,1).*triu(SA2-SB2,1);
denmat = triu(SA1+SB1,1).*triu(SA2+SB2,1);
testStat = sum(nummat(:))/sqrt(sum(denmat(:)));  

%p-value and acceptance/rejection
pval = 2*normcdf(-abs(testStat));
test = (pval<=sig);

% alternative way to compute test
% thresh = norminv(1-sig/2,0,1);
% test = abs(testStat)>thresh;