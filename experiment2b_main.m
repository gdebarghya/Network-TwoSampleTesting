function experiment2b_main (trial)

% trial = number of independent runs (we use trial = 1000)
n_all = 100:100:1000; % number of vertices
rho_all = [0.25 0.5 1 2 4]; % controls sparsity of graphs
bs = 200; % number of bootstrap samples generated
sig = 0.05; % significance level
r = 2; % r = rank parameter used by ASE and Population adjacency 
% tests, and controls number of blocks for Tracy-Widom test 


p = 0.1; q = 0.05; % block model edge probabilities under null when rho=1
eps = 0.04; % controls separation between two models under alternative 
model1 = struct('name','2SBM','n',0,'k',2,'p',0,'q',0);
model2 = struct('name','2SBM','n',0,'k',2,'p',0,'q',0);

%% Null hypothesis (both models same)

power0_ASEBoot = zeros(length(n_all),length(rho_all));
power0_AdjBoot = zeros(length(n_all),length(rho_all));
power0_TWTest = zeros(length(n_all),length(rho_all));

tic
for in = 1:length(n_all)
    model1.n = n_all(in); % n varies
    for ir = 1:length(rho_all)
        model1.p = rho_all(ir)*p; model1.q = rho_all(ir)*q; % rho varies

        disp(int2str([n_all(in) rho_all(ir)]))
        pow = runTests(r,trial,model1,model1,sig,bs);
        power0_ASEBoot(in,ir) = pow(1);
        power0_AdjBoot(in,ir) = pow(2);
        power0_TWTest(in,ir) = pow(3);
        toc
    end
    
    save(strcat('results/temp.mat'))
end

clear pow ir in
[power0_ASEBoot power0_AdjBoot power0_TWTest]
save(strcat('results/expt2b_trial',int2str(trial),'.mat'))


%% Alternative hypothesis (model2.p = model1.p + eps) 

power1_ASEBoot = zeros(length(n_all),length(rho_all));
power1_AdjBoot = zeros(length(n_all),length(rho_all));
power1_TWTest = zeros(length(n_all),length(rho_all));

for in = 1:length(n_all)
    model1.n = n_all(in); model2.n = n_all(in);
    for ir = 1:length(rho_all)
        model1.p = rho_all(ir)*p; model1.q = rho_all(ir)*q;
        model2.p = rho_all(ir)*(p+eps); model2.q = rho_all(ir)*q;

        disp(int2str([n_all(in) rho_all(ir)]))
        pow = runTests(r,trial,model1,model2,sig,bs);
        power1_ASEBoot(in,ir) = pow(1);
        power1_AdjBoot(in,ir) = pow(2);
        power1_TWTest(in,ir) = pow(3);
        toc
    end
    
    save(strcat('results/temp.mat'))
end

clear pow ir in
[power1_ASEBoot power1_AdjBoot power1_TWTest]
save(strcat('results/expt2b_trial',int2str(trial),'.mat'))

elapsedTime = toc;

function pow = runTests(r,trial,model1,model2,sig,bs)
% this function runs all the tests. We put it as a separate function to
% pararellize the loop over trials

ASEBoot = zeros(trial,1);
AdjBoot = zeros(trial,1);
TWTest = zeros(trial,1);

parfor t = 1:trial
    A = genSparseGraph(1,model1); % sample size, m=1
    B = genSparseGraph(1,model2);
        
    warning('off','MATLAB:svds:MultNotCorrectLargest');
    % suppresses warning that svds may not return correct singular values.
    % This occurs often if rank parameter r is not set correctly
    
    [ASEBoot(t),AdjBoot(t),~,~] = LowRankTests(A,B,r,sig,bs); 
    % bootstrapped tests assuming low rank structure 
    % includes bootstrap test of Tang et al (arxiv:1403.7249) using ASE and
    % our suggested improvement based on estimating population adjacency
    % r = rank of pulation adjacency (used by both tests)
     
    [TWTest(t),~] = TracyWidomTest(A,B,r,sig);    
    % test proposed in paper for m=1 (based on Tracy-Widom law)
    % r = number of communities used for SBM approximation

end

pow = mean([ASEBoot AdjBoot TWTest],1);
