function [test,pval] = NormalityTest(A,B,sig)

% Returns acceptance/rejection for the proposed normality based test (Asymp-Normal)
% Note: 
% All graphs are assumed to unweighted, undirected, and defined on a common vertex set. 
% Sample size in each population must be at least 2 
%
% Input:
% A: cell array containing networks in 1st population; each cell is a sparse adjacency matrix
% B: cell array containing networks in 2nd population, each cell being a sparse adjacency matrix
% sig: significance level for acceptance of null hypothesis
%
% Output:
% test: 1 if null is rejected, 0 otherwise
% pval: p-value for the test

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
