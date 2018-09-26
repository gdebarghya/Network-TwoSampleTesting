function experiment1a_main (trial)
% Runs the experiments for Figure 1 in the paper
% Input:
% trial = number of independent runs (we use trial = 1000)
%
% Output (saved in mat file):
% Test power (percentage of rejections) for Boot-Frobenius, Boot-Spectral, Asymp-Normal, Asymp-Chi2

n_all = 100:100:1000; % number of vertices
m_all = [2 4]; % sample size
bs = 200; % number of bootstrap samples generated
sig = 0.05; % significance level

p = 0.1; q = 0.05; % block model edge probabilities under null
eps = 0.04; % controls separation between two models under alternative 
model1 = struct('name','2SBM','n',0,'k',2,'p',p,'q',q);
model2 = struct('name','2SBM','n',0,'k',2,'p',p+eps,'q',q);

%% Null hypothesis (both models same)

power0_FroShuff = zeros(length(n_all),length(m_all));
power0_OpShuff = zeros(length(n_all),length(m_all));
power0_NorTest = zeros(length(n_all),length(m_all));
power0_ChiTest = zeros(length(n_all),length(m_all));

tic
for in = 1:length(n_all)
    model1.n = n_all(in); % n varies
    for im = 1:length(m_all)
        disp(int2str([n_all(in) m_all(im)]))
        pow = runTests(m_all(im),trial,model1,model1,sig,bs);
        power0_FroShuff(in,im) = pow(1);
        power0_OpShuff(in,im) = pow(2);
        power0_NorTest(in,im) = pow(3);
        power0_ChiTest(in,im) = pow(4);
        toc
    end
    save(strcat('results/temp.mat'))
end

clear pow im in
[power0_FroShuff power0_OpShuff power0_NorTest power0_ChiTest]
save(strcat('results/expt1a_trial',int2str(trial),'.mat'))


%% Alternative hypothesis (model2.p = model1.p + eps) 

power1_FroShuff = zeros(length(n_all),length(m_all));
power1_OpShuff = zeros(length(n_all),length(m_all));
power1_NorTest = zeros(length(n_all),length(m_all));
power1_ChiTest = zeros(length(n_all),length(m_all));

for in = 1:length(n_all)
    model1.n = n_all(in); model2.n = n_all(in);
    for im = 1:length(m_all)
        disp(int2str([n_all(in) m_all(im)]))
        pow = runTests(m_all(im),trial,model1,model2,sig,bs);
        power1_FroShuff(in,im) = pow(1);
        power1_OpShuff(in,im) = pow(2);
        power1_NorTest(in,im) = pow(3);
        power1_ChiTest(in,im) = pow(4);
        toc
    end
    save(strcat('results/temp.mat'))
end

clear pow im in
[power1_FroShuff power1_OpShuff power1_NorTest power1_ChiTest]
save(strcat('results/expt1a_trial',int2str(trial),'.mat'))

elapsedTime = toc;

function pow = runTests(m,trial,model1,model2,sig,bs)
% this function runs all the tests. We put it as a separate function to
% pararellize the loop over trials

FroShuff = zeros(trial,1);
OpShuff = zeros(trial,1);
NorTest = zeros(trial,1);
ChiTest = zeros(trial,1);

if m==1
 return
end

parfor t = 1:trial
    A = genSparseGraph(m,model1);
    B = genSparseGraph(m,model2);
        
    [FroShuff(t),OpShuff(t),~,~] = ShufflingTests(A,B,sig,bs); 
    % permutation based bootstrapped variants of tests in Ghoshdastidar et
    % al (arxiv:1707.00833)
        
    [NorTest(t),~] = NormalityTest(A,B,sig);    
    % test proposed in paper for m>=2 based on asymptotic normality
    
    [ChiTest(t),~] = GraphChi2Test(A,B,sig);    
    % chi2-type test similar to Ginestet et al (arxiv:1407.5525)

end

pow = mean([FroShuff OpShuff NorTest ChiTest],1);
