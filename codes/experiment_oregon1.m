clearvars
clc
load('oregon_networks.mat')
numgrps = 9;
numnode = size(oregonNetworks{1,1},1);

pval = zeros(numgrps); 

%% experiment oregon_1a
for i = 1:numgrps
    for j = 1:numgrps
        [~,pval(i,j)] = NormalityTest(oregonNetworks(i,:),oregonNetworks(j,:),0.05);
    end
end
nlpv = -log(pval)

%% experiment oregon_1b
trials = 100;
p_all = 0.2:0.025:0.4;
k = round(0.01*numnode)*ones(size(p_all));
[pow_p,pval_p] = detectPlanting(oregonNetworks,k,p_all,trials,numnode,numgrps);

mean(-log(pval_p),3)

%% experiment oregon_1c
k_all = round((0.01:0.001:0.02)*numnode);
p = 20./k_all;
[pow_k,pval_k] = detectPlanting(oregonNetworks,k_all,p,trials,numnode,numgrps);

mean(-log(pval_k),3)

save(strcat('results/expt_oregon1_trial',int2str(trials),'.mat'))

%%
function [pow,pval] = detectPlanting(oregonNetworks,k_all,p_all,trials,numnode,numgrps)
pow = zeros(length(k_all),numgrps,trials);
pval = zeros(length(k_all),numgrps,trials);
tic
for t = 1:trials
    if mod(t,100)==0
        disp(int2str(t)), toc
    end
    
    for j = 1:length(k_all)
        k = k_all(j);
        p = p_all(j);
        
        % choose a subset of vertices
        ind = randperm(numnode,k);
        
        % Planting a sub-graph
        %planting = sparse(ones(k)-eye(k));  % planted clique
        planting = sprand(k,k,p);
        planting = triu(planting,1); planting = planting + planting';
        
        for i = 1:numgrps
            temp = oregonNetworks(i,:);
            for l = 1:length(temp)
                temp{l}(ind,ind) = planting;
            end
            [pow(j,i,t),pval(j,i,t)] = NormalityTest(oregonNetworks(i,:),temp,0.05);
        end
    end
end
end