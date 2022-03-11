function [clusterNamesUp,clusterNamesDown,net8angle_Up,net8angle_Down,numNets, YeoNetNames] = NAME_CLUSTERS_UP_DOWN(centroids)

%Provide names for clusters based on angular distance to binary Yeo
%System Vectors
%returns vector where 1 indicates a "+" state and 0 indicates a "-" state
%centroids = kClusterCentroids;

[nparc,numClusters] = size(centroids);
if nparc == 333
    load('cuttingEdgeFC/gorgon_networklabels.mat'); network7labels = gordon_numbers(1:nparc,2);
    numNets=12;
    YeoNetNames={'default mode', 'somatomotor', 'visual', 'fronto-parietal', 'auditory', 'cingulo-parietal', 'retrosplenial-temporal', 'ventral-attn', 'dorsal-attn', 'salience', 'cingulo-opercular', 'none'}
end
if nparc == 392
    tmp=readtable('data/CC400_Yeo7_Map.csv');
    network7labels=tmp.Yeo7Label;
    numNets=9;
    YeoNetNames = {'VIS', 'SOM', 'ATTN-DORS', 'ATTN-VEN', 'LIMB', 'FPN', 'DMN', 'SUB', 'CEREB'}
elseif nparc == 86
    tmp=load('data/fs86_to_yeo_map.csv')';
    numNets=9;
    network7labels=tmp;
    YeoNetNames = {'VIS', 'SOM', 'ATTN-DORS', 'ATTN-VEN', 'LIMB', 'FPN', 'DMN', 'SUB', 'CEREB'}

elseif nparc == 268
    load('networklabels.mat'); network7labels = networks_net(1:nparc);
    numNets=8;
    YeoNetNames = {'MED FRONT', 'FPN', 'DMN', 'SUB', 'MOTOR', 'VIS I', 'VIS II', 'VIS III'};
end

network8labels=network7labels;
binaryNetVectors = ones(nparc,numNets) .* repmat((1:numNets),[nparc 1]); 
binaryNetVectors = double(binaryNetVectors == network8labels);


% calculate cosine of angle between binary state vector and centroids

centroids_up = centroids .* (centroids > 0);
centroids_down = -1 * centroids .* (centroids < 0);     % make negative activity positive and get rid of positive activity

net8angle_Up = zeros(numClusters,numNets);
net8angle_Down = zeros(numClusters,numNets);

for K = 1:numClusters
    for B = 1:numNets
        net8angle_Up(K,B) = dot(centroids_up(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
        net8angle_Down(K,B) = dot(centroids_down(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
    end
end

% get index of minimum and assign names

clusterNamesUp = cell(numClusters,1);
clusterNamesDown = cell(numClusters,1);
for K = 1:numClusters       % for up and down separately, calculate closest network
    Up_ind = find(net8angle_Up(K,:) == max(net8angle_Up(K,:)));
    Down_ind = find(net8angle_Down(K,:) == max(net8angle_Down(K,:)));    %cos(0) = 1 so need max not min (like for E.D.)
    clusterNamesUp{K} = [YeoNetNames{Up_ind},'+'];
    clusterNamesDown{K} = [YeoNetNames{Down_ind},'-'];
end