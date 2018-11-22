clc
clearvars
fileNames = dir('oregon');

allIndices = cell(9,2);

c1 = 0; c2 = 0; totaledge = 0;
for l = 1:length(fileNames)
    if contains(fileNames(l).name,'oregon1')
        c1 = c1 + 1; i = c1; j = 1;
    elseif contains(fileNames(l).name,'oregon2')
        c2 = c2 + 1; i = c2; j = 2;
    else
        continue
    end
    
    fileID = fopen(strcat(fileNames(l).folder,'/',fileNames(l).name));
    allIndices{i,j} = textscan(fileID,'%d %d','HeaderLines',5); 
    totaledge = totaledge + length(allIndices{i,j}{1});
    % note: the first pair on line 5 is actually redundant
    fclose(fileID);
end

% determine number of nodes and relabel all nodes
allnodes = zeros(2*totaledge,1);
c = 0;
for i = 1:size(allIndices,1)
    for j = 1:size(allIndices,2)
        l = length(allIndices{i,j}{1});
        allnodes((c+1):(c+l)) = allIndices{i,j}{1};
        c = c + l;
        allnodes((c+1):(c+l)) = allIndices{i,j}{2};
        c = c + l;
    end
end

[relabnodes,labnames] = grp2idx(allnodes);
numnode = length(labnames); % total distinct nodes

% relabelling all nodes
c = 0;
for i = 1:size(allIndices,1)
    for j = 1:size(allIndices,2)
        l = length(allIndices{i,j}{1});
        allIndices{i,j}{1}(:) = relabnodes((c+1):(c+l));
        c = c + l;
        allIndices{i,j}{2}(:) = relabnodes((c+1):(c+l));
        c = c + l;
    end
end

% create sparse networks. one can also combine this step with previous
% deg = zeros(numnode,1);
oregonNetworks = cell(9,2);
combinedgraph = sparse([],[],[],numnode,numnode,2*totaledge);

for i = 1:size(allIndices,1)
    for j = 1:size(allIndices,2)
        temp = allIndices{i,j};
        oregonNetworks{i,j} = sparse(double([temp{1}; temp{2}]),...
                          double([temp{2}; temp{1}]),1,numnode,numnode);
        combinedgraph = combinedgraph|oregonNetworks{i,j};
%        deg = deg + sum(oregonNetworks{i,j},2);
    end
end    

% community detection from all networks combined
[I,J,~] = find(triu(combinedgraph));

fileID = fopen('temp-graph.txt','w');
fprintf(fileID,'# Undirected graph \n# Union of all oregon graphs\n');
fprintf(fileID,'# Nodes: %d Edges: %d \n# FromNodeId \tToNodeId',numnode,length(I));
for l = 1:length(I)
    fprintf(fileID,'\n%-5d\t%d',I(l),J(l));
end
fclose(fileID);

% run bigclam from command line
!snap-master/examples/bigclam/bigclam -i:temp-graph.txt -c:50
numcomm = 50;  % have not tried to pass numcomm to command line

% read labels
comms = cell(numcomm,1);
fileID = fopen('cmtyvv.txt','r');
textLine = fgets(fileID); 
l = 1;
while ischar(textLine)
    comms{l} = sscanf(textLine,'%d ');
    comms{l} = unique(comms{l});
    l = l+1;
    textLine = fgets(fileID); 
end
fclose(fileID);

commDeg = zeros(1,numcomm);
commID = zeros(numnode,1);
for i = 1:numnode
    for l = 1:numcomm
        commDeg(l) = sum(comms{l}==i)*sum(combinedgraph(i,comms{l}));
        % if i in community, then return its degree within that community
    end 
    [mx,commID(i)] = max(commDeg);
    if sum(commDeg==mx)>1
        ind = find(commDeg==mx);
        j = randperm(length(ind),1);
        commID(i) = j;
    end
end

% remove additionally created files and save data
!rm cmtyvv.txt
!rm temp-graph.txt
!rm graph.gexf
save('oregon_networks.mat','oregonNetworks','combinedgraph','commID');