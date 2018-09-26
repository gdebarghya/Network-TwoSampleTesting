function [test,pval] = GraphChi2Test(A,B,sig)

% Returns acceptance/rejection for a chi2-type test (AsympChi2)
% Note: All graphs are assumed to unweighted, undirected, and defined on a common vertex set
%
% Input:
% A: cell array containing networks in 1st population; each cell is a sparse adjacency matrix
% B: cell array containing networks in 2nd population, each cell being a sparse adjacency matrix
% sig: significance level for acceptance of null hypothesis
%
% Output:
% test: 1 if null is rejected, 0 otherwise
% pval: p-value for the test

mA = length(A); mB = length(B);
n = size(A{1},1); 

% Vectorisation of adjacencies
vecA = zeros(mA,n^2);
for i = 1:mA
    C = triu(A{i},1);
    vecA(i,:) = C(:).';
end
vecB = zeros(mB,n^2);
for i = 1:mB
    C = triu(B{i},1);
    vecB(i,:) = C(:).';
end
ind = (sum(vecA+vecB,1)~=0);
vecA = vecA(:,ind); vecB = vecB(:,ind); 
% removes everything below and on diagonal, and also those edges that are 
% not observed in any graph in the population

meanDiff = mean(vecA,1) - mean(vecB,1);
varAll = var(vecA)/mA + var(vecB)/mB;

if sum((meanDiff~=0).*(varAll==0))>0   % ?/0 case
    pval = 0;
    test = 1;
    return
end

ind = (varAll~=0);  
meanDiff = meanDiff(ind); varAll = varAll(ind); % to avoid 0/0 or ?/0
% assume 0 < Pij, Qij < 1 and hence zero variance is not theoretically possible 

%testStat = (mA*mB/(mA+mB))*sum((meanDiff.^2)./varAll) % test statistic with pooled variance
testStat = sum((meanDiff.^2)./varAll);  

%p-value and acceptance/rejection
pval = 1 - chi2cdf(testStat,n*(n-1)/2);
test = (pval<=sig);

% alternative way to compute test
% thresh = chi2inv(1-sig,n*(n-1)/2)
% test = testStat>thresh
