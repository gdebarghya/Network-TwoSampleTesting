function experiment1b_main (trial)

% trial = number of independent runs (we use trial = 1000)
n = 500; % number of vertices
m_all = 2:2:10; % sample size
rho_all = [0.25 0.5 1 2 4]; % controls sparsity of graphs
bs = 200; % number of bootstrap samples generated
sig = 0.05; % significance level

p = 0.1; q = 0.05; % block model edge probabilities under null when rho=1
eps = 0.04; % controls separation between two models under alternative 
model1 = struct('name','2SBM','n',n,'k',2,'p',0,'q',0);
model2 = struct('name','2SBM','n',n,'k',2,'p',0,'q',0);

%% Null hypothesis (both models same)

power0_FroShuff = zeros(length(rho_all),length(m_all));
power0_OpShuff = zeros(length(rho_all),length(m_all));
power0_NorTest = zeros(length(rho_all),length(m_all));
power0_ChiTest = zeros(length(rho_all),length(m_all));

tic
for ir = 1:length(rho_all)
    model1.p = rho_all(ir)*p; model1.q = rho_all(ir)*q; % rho varies
    for im = 1:length(m_all)
        disp(int2str([rho_all(ir) m_all(im)]))
        pow = runTests(m_all(im),trial,model1,model1,sig,bs); 
        power0_FroShuff(ir,im) = pow(1);
        power0_OpShuff(ir,im) = pow(2);
        power0_NorTest(ir,im) = pow(3);
        power0_ChiTest(ir,im) = pow(4);
        toc
    end
    save(strcat('results/temp.mat'))
end

clear pow im ir
[power0_FroShuff power0_OpShuff power0_NorTest power0_ChiTest]
save(strcat('results/expt1b_trial',int2str(trial),'.mat'))


%% Alternative hypothesis (model2.p = model1.p + rho*eps) 

power1_FroShuff = zeros(length(rho_all),length(m_all));
power1_OpShuff = zeros(length(rho_all),length(m_all));
power1_NorTest = zeros(length(rho_all),length(m_all));
power1_ChiTest = zeros(length(rho_all),length(m_all));

for ir = 1:length(rho_all)
    model1.p = rho_all(ir)*p; model1.q = rho_all(ir)*q;
    model2.p = rho_all(ir)*(p+eps); model2.q = rho_all(ir)*q;
    for im = 1:length(m_all)
        disp(int2str([rho_all(ir) m_all(im)]))
        pow = runTests(m_all(im),trial,model1,model2,sig,bs);
        power1_FroShuff(ir,im) = pow(1);
        power1_OpShuff(ir,im) = pow(2);
        power1_NorTest(ir,im) = pow(3);
        power1_ChiTest(ir,im) = pow(4);
        toc
    end
    save(strcat('results/temp.mat'))
end

clear pow im ir
[power1_FroShuff power1_OpShuff power1_NorTest power1_ChiTest]
save(strcat('results/expt1b_trial',int2str(trial),'.mat'))

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
