function [test,pval] = GraphChi2Test(A,B,sig)

% returns acceptance/rejection for a chi2-type test
% assumes A, B are defined as a cell array of length m, each cell being an
% adjacency matrix (stored as a sparse matrix)

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