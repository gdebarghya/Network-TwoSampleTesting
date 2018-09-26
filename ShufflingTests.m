function [tFro,tOp,pvalFro,pvalOp] = ShufflingTests(A,B,sig,bs)

% permutation based bootstrapped variants of tests in arxiv 1707.00833
% assumes A, B are defined as a cell array of length m, each cell being an
% adjacency matrix (stored as a sparse matrix)

m = min(length(A),length(B));
C = cat(1,A,B); 

testStat = computeStat(A(1:m),B(1:m));

%Bootstrapping via randomly permuting all labels 
bsStat = zeros(bs,2);
for b = 1:bs
    ind = randperm(length(C));
    bsStat(b,:) = computeStat(C(ind(1:m)),C(ind((m+1):(2*m)))); 
end
bsStat = sort(bsStat,'descend');

%p-values and acceptance/rejection
pvalFro = (sum(bsStat(:,1)>=testStat(1))+0.5)/bs; % +0.5 for continuity correction
pvalOp = (sum(bsStat(:,2)>=testStat(2))+0.5)/bs;
tFro = (pvalFro<=sig);
tOp = (pvalOp<=sig);

%alternatively
% thresh = bsStat(floor(sig*bs),:)
% tFro = testStat(1)>thresh(1);
% tOp = testStat(2)>thresh(2);

%%%

function stats = computeStat(A1,B1)
m1 = floor(length(A1)/2);
n = size(A1{1},1);

SA1 = sparse(n,n); SA2 = sparse(n,n); SB1 = sparse(n,n); SB2 = sparse(n,n);
for i = 1:m1
    SA1 = SA1 + A1{i}; SA2 = SA2 + A1{i+m1};
    SB1 = SB1 + B1{i}; SB2 = SB2 + B1{i+m1};
end
nummat = triu(SA1-SB1,1).*triu(SA2-SB2,1);
denmat = triu(SA1+SB1,1).*triu(SA2+SB2,1);
stats(1) = sum(nummat(:))/sqrt(sum(denmat(:)));  % frobenius norm stat
stats(2) = svds(SA1+SA2-SB1-SB2,1)/sqrt(norm(SA1+SA2+SB1+SB2,1)); % operator norm stat


    
    
    