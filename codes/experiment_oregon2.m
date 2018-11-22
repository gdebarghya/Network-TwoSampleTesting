clearvars
clc
load('oregon_networks.mat')
numgrps = 18;

numnode = size(oregonNetworks{1,1},1);
oregonNetworks = oregonNetworks';
oregonNetworks = oregonNetworks(:);

%% experiment oregon2 (randomly perturb few connections)
trials = 100;
e_all = 0:25:300; 
pow_e = zeros(length(e_all),numgrps,trials);
pval_e = zeros(length(e_all),numgrps,trials);

tic
parfor t = 1:trials
    if mod(t,10)==0
        disp(int2str(t))
    end
    [pow_e(:,:,t),pval_e(:,:,t)] = detectPerturb(oregonNetworks,commID,e_all,numnode,numgrps);
end
toc
pow_all = mean(pow_e,3); 
pval_all = mean(pval_e,3);
[mean(pow_all,2) mean(pval_all,2)]

save(strcat('results/expt_oregon2_trial',int2str(trials),'.mat'))

%%
function [pow,pval] = detectPerturb(oregonNetworks,commID,e_all,numnode,numgrps)
pow = zeros(length(e_all),numgrps);
pval = zeros(length(e_all),numgrps);

for j = 1:length(e_all)
    ne = e_all(j);
    
    % Perturbation matrix
    ind = randi(numnode,ne,2);
    pert = sparse(ind(:,1),ind(:,2),1,numnode,numnode);
    pert = double(pert|pert');
    pert = pert - diag(diag(pert));
    
    
    for i = 1:numgrps
        temp = oregonNetworks(i);
        temp{1} = mod(temp{1}+pert,2);
        [pow(j,i),pval(j,i)] = TracyWidomTest(oregonNetworks(i),temp,commID,0.05);
    end
end
end