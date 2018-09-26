function [tASE,tAdj,pvalASE,pvalAdj] = LowRankTests(A,B,r,sig,bs)


testStat = computeStat(A{1},B{1},r);

% bootstrap
bsStat = zeros(bs,2); % we do twice, once from A and once from B

%%% generate samples from E[A]
[u,s,v] = svds(double(A{1}),r); s(isnan(s))=0; u(isnan(u))=0; v(isnan(v))=0;
EA = u*s*v'; EA = EA - diag(diag(EA)); EA = min(1,max(0,EA));
model = struct('name','IER','n',size(EA,1),'P',EA); 
for b = 1:bs
    C = genSparseGraph(2,model);
    bsStat(b,:) = computeStat(C{1},C{2},r);
end
px = (sum(bsStat>repmat(testStat,bs,1),1)+0.5)/bs;

%%% generate samples from E[B]
[u,s,v] = svds(double(B{1}),r); s(isnan(s))=0; u(isnan(u))=0; v(isnan(v))=0;
EB = u*s*v'; EB = EB - diag(diag(EB)); EB = min(1,max(0,EB));
model.P = EB;
for b = 1:bs
    C = genSparseGraph(2,model);
    bsStat(b,:) = computeStat(C{1},C{2},r);
end
py = (sum(bsStat>repmat(testStat,bs,1),1)+0.5)/bs;

%p-values and acceptance/rejection
pvalASE = max(px(1),py(1));
tASE = (pvalASE<=sig);
pvalAdj = max(px(2),py(2));
tAdj = (pvalAdj<=sig);


function stats = computeStat(A1,B1,r)

[u1,s1,v1] = svds(double(A1),r); s1(isnan(s1))=0; u1(isnan(u1))=0; v1(isnan(v1))=0;
[u2,s2,v2] = svds(double(B1),r); s2(isnan(s2))=0; u2(isnan(u2))=0; v2(isnan(v2))=0; 

% ASE statistic
% min_W||X-YW||_F solved using orthogonal Procrustes problem, whose solution 
% is Wopt = UV' where X'Y = USV' (see wikipedia)
X = u1*sqrt(s1); Y = u2*sqrt(s2);
[u,~,v] = svd(X'*Y); W = u*v'; 
stats(1) = norm(X-Y*W,'fro');

% EPA statistic (difference of estimated population adjacencies)
stats(2) = norm(u1*s1*v1'-u2*s2*v2','fro');

