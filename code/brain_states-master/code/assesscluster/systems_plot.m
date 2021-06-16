addpaths;
masterdir = fullfile(basedir,'results',name_root);

%%
load(fullfile(masterdir,'clusterAssignments',['k',num2str(numClusters),name_root,'.mat']));
overallPartition = clusterAssignments.(['k',num2str(numClusters)]).partition;
centroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
overallNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
savedir = fullfile(masterdir,'analyses','centroids');
mkdir(savedir);

%%
nparc=268
six=centroids(:,6);
one=centroids(:,1);
rest=centroids(:,2:5);

new=[one,six,rest];
centroids=new
[nparc,numClusters] = size(centroids);

[~,~,~,net8angle] = NAME_CLUSTERS_ANGLE(centroids);

%%
[up,dn,net8angle_Up,net8angle_Down] = NAME_CLUSTERS_UP_DOWN(centroids)
YeoNetNames = {'MED FRONT', 'FPN', 'DMN', 'SUB', 'MOTOR', 'VIS I', 'VIS II', 'VIS III'};
numNets = numel(YeoNetNames);

%% make radial plots
clusterColors = GET_CLUSTER_COLORS(numClusters);
numNets=8
clusterColors = hex2rgb(clusterColors);
netAngle = linspace(0,2*pi,numNets+1);
thetaNames = YeoNetNames; thetaNames{9} = '';
f=figure;
for K = 1:numClusters
    ax = subplot(1,numClusters,K,polaraxes); hold on
    polarplot(netAngle,[net8angle_Up(K,:) net8angle_Up(K,1)],'k');
    polarplot(netAngle,[net8angle_Down(K,:) net8angle_Down(K,1)],'r');
    thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
    rticks([0.4 0.8]); rticklabels({'0.4','0.8'});
%     for L = 1:numNets
%         ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
%         YeoColors(L,:), ax.ThetaTickLabel{L});
%     end
    set(ax,'FontSize',6);
    %title(overallNames{K},'Color',clusterColors(K,:),'FontSize',8);
end


saveas(f,fullfile(savedir,['SystemsRadial_k',num2str(numClusters),name_root,'.pdf']));

% for source data file
save(fullfile(savedir,['Fig2b__YeoSystemAlignment_k',num2str(numClusters),'.mat']),'netAngle','net7angle_Up','net7angle_Down');
