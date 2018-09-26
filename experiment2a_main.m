function experiment2a_main (trial)

% trial = number of independent runs (we use trial = 1000)
n_all = 100:100:1000; % number of vertices
bs = 200; % number of bootstrap samples generated
sig = 0.05; % significance level
r_all = [2 4]; % r = rank parameter used by ASE and Population adjacency 
% tests, and controls number of blocks for Tracy-Widom test 


p = 0.1; q = 0.05; % block model edge probabilities under null
eps = 0.04; % controls separation between two models under alternative 
model1 = struct('name','2SBM','n',0,'k',2,'p',p,'q',q);
model2 = struct('name','2SBM','n',0,'k',2,'p',p+eps,'q',q);

%% Null hypothesis (both models same)

power0_ASEBoot = zeros(length(n_all),length(r_all));
power0_AdjBoot = zeros(length(n_all),length(r_all));
power0_TWTest = zeros(length(n_all),length(r_all));

tic
for in = 1:length(n_all)
    model1.n = n_all(in); % n varies
    for ir = 1:length(r_all)
        disp(int2str([n_all(in) r_all(ir)]))
        pow = runTests(r_all(ir),trial,model1,model1,sig,bs);
        power0_ASEBoot(in,ir) = pow(1);
        power0_AdjBoot(in,ir) = pow(2);
        power0_TWTest(in,ir) = pow(3);
        toc
    end
end

clear pow ir in
[power0_ASEBoot power0_AdjBoot power0_TWTest]
save(strcat('results/expt2a_trial',int2str(trial),'.mat'))


%% Alternative hypothesis (model2.p = model1.p + eps) 

power1_ASEBoot = zeros(length(n_all),length(r_all));
power1_AdjBoot = zeros(length(n_all),length(r_all));
power1_TWTest = zeros(length(n_all),length(r_all));

for in = 1:length(n_all)
    model1.n = n_all(in); model2.n = n_all(in);
    for ir = 1:length(r_all)
        disp(int2str([n_all(in) r_all(ir)]))
        pow = runTests(r_all(ir),trial,model1,model2,sig,bs);
        power1_ASEBoot(in,ir) = pow(1);
        power1_AdjBoot(in,ir) = pow(2);
        power1_TWTest(in,ir) = pow(3);
        toc
    end
end

clear pow ir in
[power1_ASEBoot power1_AdjBoot power1_TWTest]
save(strcat('results/expt2a_trial',int2str(trial),'.mat'))

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
    % our suggested improvement based on estimating population adjacency (EPA)
    % r = rank of population adjacency (used by both tests)
     
    [TWTest(t),~] = TracyWidomTest(A,B,r,sig);    
    % test proposed in paper for m=1 (based on Tracy-Widom law)
    % r = number of communities used for SBM approximation

end

pow = mean([ASEBoot AdjBoot TWTest],1);
