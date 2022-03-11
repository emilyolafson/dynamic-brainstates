function [net7angle, numNets, YeoNetNames] = NAME_CLUSTERS_ANGLE(centroids)
%Provide names for clusters based on angular distance to binary Yeo
%System Vectors
%returns vector where 1 indicates a "+" state and 0 indicates a "-" state
%centroids = kClusterCentroids;
%Contains Yeo parcellations for fs86, shen268, gordon333, cc400

[nparc,numClusters] = size(centroids);

if nparc == 333
    load('cuttingEdgeFC/gorgon_networklabels.mat');
    network7labels = gordon_numbers(1:nparc,2);
    numNets=12;
    YeoNetNames={'default mode', 'somatomotor', 'visual', 'fronto-parietal', 'auditory', 'cingulo-parietal', 'retrosplenial-temporal', 'ventral-attn', 'dorsal-attn', 'salience', 'cingulo-opercular', 'none'}
end
if nparc == 392
    disp('hello')
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
    load('networklabels.mat');
    network7labels = networks_net(1:nparc);
    numNets=8;
    YeoNetNames = {'MED FRONT', 'FPN', 'DMN', 'SUB', 'MOTOR', 'VIS I', 'VIS II', 'VIS III'};
end
% make a matrix where each column corresponds to a labeled Yeo system in Lausanne parcellation
% the columns are binary vectors indicated whether a region belongs to corresponding Yeo system

binaryNetVectors = ones(nparc,numNets) .* repmat((1:numNets),[nparc 1]); 
binaryNetVectors = double(binaryNetVectors == network7labels);
% then duplicate this matrix, multiply by -1 and horizontally concatenate to
% provide separate names for when systems are low amplitude

binaryNetVectors = [binaryNetVectors, -1*binaryNetVectors];

% calculate cosine of angle between binary state vector and centroids

net7angle = zeros(numClusters,numNets*2);

for K = 1:numClusters
    for B = 1:(numNets*2)
        net7angle(K,B) = dot(centroids(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
    end
end

% get index of minimum and assign names

clusterNamesInit = cell(numClusters,1);
plusminus = true(numClusters,1);
for K = 1:numClusters
    ind = find(net7angle(K,:) == max(net7angle(K,:)));    %cos(0) = 1 so need max not min (like for E.D.)
    %clusterNamesInit{K} = YeoNetNames{ind};
    if ind > numNets
        plusminus(K) = false;
    end
end

%clusterNames = cellstr(clusterNamesInit);

%sort by name then plus-minus
%
%[clusterNamesSort,I] = sort(clusterNamesInit);
%[~,I2] = sort(plusminus(I));
%clusterNamesSort = clusterNamesSort(I2);
%reorderClusters = I(I2);
%}
