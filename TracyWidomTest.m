function [test,pval] = TracyWidomTest(A,B,r,sig)

% A,B = two cell arrays, each containing a large sparse matrix (both of same size)
% r = number of communities or can be a vector with all community labels
% sig = significance level

n = size(A{1},1);

if length(r)==1
    idx = spectralClustering(A,B,r);
else
    % if we precompute communities and pass it as r. In this case, we
    % assume that the community id-s start from 1 to a maximum value
    idx = r;                    
    r = max(idx);
end
    

% computing the scaled difference matrix
C = A{1}-B{1}; 
for i = 1:r
    for j = 1:r
        if i==j
            temp = A{1}(idx==i,idx==j);
            Pij = sum(temp(:))/(size(temp,1)*size(temp,1)-1); % takes into account the diagonal is zero 
            temp = B{1}(idx==i,idx==j);
            Qij = sum(temp(:))/(size(temp,1)*size(temp,1)-1);
        else
            temp = A{1}(idx==i,idx==j);
            Pij = mean(temp(:));
            temp = B{1}(idx==i,idx==j);
            Qij = mean(temp(:));
        end
        denom = sqrt((n-1)*(Pij*(1-Pij)+Qij*(1-Qij)));
        if denom==0             %this should not occur unless graphs are too sparse
            denom = 1e-5;
        end
        C(idx==i,idx==j) = C(idx==i,idx==j)/denom;
    end
end
C(isnan(C)) = 0;

% compute test statistic and p-value
testStat = n^(2/3)*(svds(C,1)-2);
pval = min(1,2*computeTWpval(testStat));
test = (pval<=sig);  
% we double p-value for Bonferroni correction. 
% We know lambda_1 and -lambda_n are TW.
% Spectral norm = max(lambda_1, -lambda_n), and hence, we need Bonferroni


function p = computeTWpval(x)
load TW_beta1_CDF.mat TW_arg TW_CDF
ind = max(1,sum(TW_arg<=x));
%disp(int2str([ind length(TW_CDF)]))
p = 1 - TW_CDF(ind);

function idx = spectralClustering(A,B,r)
% spectral clustering to find common block structure
C = (A{1}+B{1})/2;
d = sum(C,2); d(d==0) = 1; d = 1./sqrt(d);
C = C.*(d*d'); 
[~,~,vec] = svds(double(C),r);  % we use dominant singular vectors to tackle both homophilic and heterophilic cases
normv = sqrt(sum(vec.^2,2)); normv(normv==0) = 1; 
vec = vec./repmat(normv,1,r);
idx = kmeans(vec,r);
