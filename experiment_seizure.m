clc
clearvars

%% Read data and set parameters
datalab = csvread('seizure.csv',1,1);
numinst = size(datalab,1);
numlabel = max(datalab(:,end));
numattr = size(datalab,2)-1;

trials = 1000;
numpart = 4;
numgrps = numlabel*numpart/2;

pval_m1 = zeros(numgrps,numgrps,trials);
pval_m2 = zeros(numgrps,numgrps,trials);
pow_m1 = zeros(numgrps,numgrps,trials);
pow_m2 = zeros(numgrps,numgrps,trials);

tic
for t = 1:trials
    if mod(t,100)==0
        disp(int2str(t)), toc
    end
    data = datalab(randperm(numinst),:);
    label = data(:,end);
    data = data(:,1:(end-1));
        
    %% Partition data
    groups = cell(numlabel,numpart);
    for i = 1:numlabel
        ind = find(label==i);
        blocksize = ceil(length(ind)/numpart);
        c = 1;
        for j = 1:numpart
            d = min(c+blocksize-1,length(ind));
            groups{i,j} = ind(c:d);
            c = d+1;
        end
    end
    
    %% Compute correlation and networks
    numedge = round(0.1*numattr^2);
    allgraph = cell(size(groups));
    for i = 1:numlabel
        for j = 1:numpart
            allgraph{i,j} = sparse(numattr,numattr);
            seldata = data(groups{i,j},:);
            c = corrcoef(seldata);
            d = sort(c(:),'descend');
            thresh = d(numedge+numattr);
            c = sparse(c>=thresh); c = c - diag(diag(c));
            allgraph{i,j} = c;
        end
    end
    
    allgraph1 = reshape(allgraph',2,numlabel*numpart/2).';
    allgraph2 = allgraph1(:,1);
    
    %% m = 2
    for i = 1:numgrps
        for j = i:numgrps
            [pow,pval] = NormalityTest(allgraph1(i,:),allgraph1(j,:),0.05);
            pval_m2(i,j,t) = pval;
            pval_m2(j,i,t) = pval;
            pow_m2(i,j,t) = pow;
            pow_m2(j,i,t) = pow;
        end
    end
    
    %% test for m=1
    for i = 1:numgrps
        for j = i:numgrps
            [pow,pval] = TracyWidomTest(allgraph2(i),allgraph2(j),10,0.05);
            pval_m1(i,j,t) = pval;
            pval_m1(j,i,t) = pval;
            pow_m1(i,j,t) = pow;
            pow_m1(j,i,t) = pow;
        end
    end
end

mpow1 = mean(pow_m1,3)
mpow2 = mean(pow_m2,3)
mlnp1 = mean(-log(pval_m1),3)
mlnp2 = mean(-log(pval_m2),3)

save(strcat('results/expt_seizure_trials',int2str(trials),'.mat'))
