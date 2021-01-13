%%
% Start here with your analysis: 
% load concatenated TS with [T, nparc] size.
% Specify numclusters
% 50 partitions are compared for mutual information
% the partition which scores the most mutual information with all 49 other
% cluseters is selected, plotted, and saved for further analysis

%%
clear all; close all;clc
basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);
addpath(genpath('code'))
%% set inputs

savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
distanceMethod = 'correlation'; % distance metric for clustering, we used correlation
numClusters = 6; % set number of clusters -- this requires a separate process of your choosing to select the best number of clusters. See https://arxiv.org/abs/1007.1075 for a useful discussion. 
nreps = 50;	% how many times to repeat clustering. will choose lowest error solution

%% load BOLD data

% replace TS with your BOLD data formatted as a T-by-nparc matrix
% where T is the number of time points and nparc is the number of ROIs
load LSDgsr_cat.mat TS_gsr % HCP U100 sample in lausanne250 parcellation

TR = 220;

[T,nparc] = size(TS_gsr);

%% generate 50 partitions that will be compared pair-wise for mutual information

parts = NaN(19580,50);

for i=1:50
    [parts(:,i),~,sumd] = kmeans(TS_gsr,numClusters,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',1000);
end

%% calculate adjusted mutual information for every pair of partitions

ami_results = NaN(50,50);

for i=1:50
    for j=1:50
        ami_results(i,j) = ami(parts(:,i),parts(:,j));
    end
end

% assess 
[m,ind] = max(sum(ami_results,1)); %ind corresponds to the partition which has the highest mutual information with all other partitions
partition = parts(:,ind); % take partition that has most agreement with all other for further analysis

% plot
f = figure;

imagesc(ami_results); title('Adjusted Mutal Information between Partitions'); colorbar;
axis square; set(gca,'FontSize',8); 
f.PaperUnits = 'inches';
f.PaperSize = [4 2];
f.PaperPosition = [0 0 4 2];
saveas(f,fullfile(savedir,['AMI_k',num2str(numClusters),'.pdf']));


%% compute centroids and plot
centroids = GET_CENTROIDS(TS_gsr,partition,numClusters);
% name clusters based on alignment with Yeo resting state networks
clusterNames = NAME_CLUSTERS_ANGLE(centroids);  % need to add a prior Yeo partition labels for your parcellation
[clusterNamesUp,clusterNamesDown] = NAME_CLUSTERS_UP_DOWN(centroids);  % need to add a prior Yeo partition labels for your parcellation

f = figure;
subplot(1,2,1); imagesc(centroids); title('Centroids'); xticks(1:numClusters); xticklabels(clusterNames);
colormap('plasma'); axis square; colorbar; set(gca,'FontSize',8); COLOR_TICK_LABELS(true,false,numClusters);
subplot(1,2,2); imagesc(corr(centroids)); title('Centroid Similarity'); colorbar; caxis([-1 1]); 
colormap('plasma'); axis square; set(gca,'FontSize',8); xticks(1:numClusters); yticks(1:numClusters); 
xticklabels(clusterNames); yticklabels(clusterNames); xtickangle(90);
COLOR_TICK_LABELS(true,true,numClusters);
f.PaperUnits = 'inches';
f.PaperSize = [4 2];
f.PaperPosition = [0 0 4 2];
saveas(f,fullfile(savedir,['Centroids_k',num2str(numClusters),'.pdf']));

%% save 

save(fullfile(savedir,['Partition',num2str(numClusters),'.mat']), 'parts', 'ami_results', 'ind', 'partition', 'clusterNames', 'centroids')
