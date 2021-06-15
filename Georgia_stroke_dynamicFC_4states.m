
baselineFM=load('baselineFM.mat');
baselineFM=baselineFM.basline;
save('baseline_scores.mat', 'baselineFM')

numscans1_11=[4 3 3 3 3 3 3 3 3 3 3;4 3 3 3 3 3 3 4 3 3 3;3 3 3 3 3 3 3 3 3 3 3;4 4 3 3 3 3 3 3 3 3 3 ; 3 3 3 3 3 0 3 3 3 3 3];
numscans12_23=[3 3 3 3 3 2 2 2 2 2 2 2; 3 3 3 4 3 2 2 2 2 2 2 2 ;3 3 3 3 3 2 2 2 0 2 2 2;0 3 3 3 3 2 2 2 0 2 2 2;0 3 4 3 3 2 2 2 0 2 2 2];
nscans = [numscans1_11, numscans12_23];
nscan = ones(5,24)+4;
nscans = [nscans, nscan];

%% 
allsubs=load('data/ts_shenGSR_may29.mat'); 
allsubs=allsubs.newaverg;
allsubs_fmri=reshape(allsubs', [1 235]); %reshape goes column-wise; all S1 listed, then S2,..

% Determine the optimal number of clusters using variance explained
allsubs_fmri=allsubs_fmri(find(~cellfun(@isempty,allsubs_fmri)))
allsubs_fmri=allsubs_fmri';
np= [];
count=[];
leng2=[];

for i=1:229
    disp(i)
    a=allsubs_fmri{i};
    b=cell2mat(a);
    leng_original(i)=size(b,1);
    if leng_original(i)<200
        b=b(1:leng_original(i),:);
        count=[count; i];
    else %set scans to all have the same length (10 mins)
        b=b(1:200,:);
    end
    leng2(i)=size(b,1);
    np=[np;b];
end

np=zscore(np,[],1); % z-score along columns

partition=ones(59147,4)
c=1
l=1
for i=1:47
    nsess=5;
    if i==6
        nsess=4
    end
    if i==12
        nsess=3
    end
    if i==20
        nsess=2
    end
    for j=1:nsess
        clusters2(i,j)=size(partition(c:c+leng2(l)-1),2)
        c=c+leng2(i);
        l=l+1;
    end
end
%%
distanceMethod= 'correlation'
nreps = 50;
totSum=[];

% should take about 5 mins to run for each cluster #
clear cluster_output2
clear sumd
for i=4
    [cluster_output2(:,i),~,sumd]=kmeans(np,i,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',1000);
    totSum(i)=sum(sumd);
end

% Elbow plots
for i=1:12
  disp(num2str(i))
  [~,~,sumd]=kmeans(np,i); %sumd= within-cluster sums of point-to-centroid distances
  totSum(i)=sum(sumd); % Inertia
  avgDist(i)=mean(sumd); % Distortion
end


% plot cluster quality
close all;
plot(avgDist, 'or')
hold on;
plot(avgDist, '-k')
xticks(1:12)
xticklabels({"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"})
xlabel('Number of clusters')
ylabel({'Average distance between', 'cluster centroid & member point'})
set(gca, 'FontSize', 15)



%% Visualize cluster centroids
% compute centroids and plot
best_number_of_clusters =4 %define
centroids = GET_CENTROIDS(np,cluster_output2(:,best_number_of_clusters),best_number_of_clusters);


% radial plots

[~,~,~,net8angle] = NAME_CLUSTERS_ANGLE(centroids);
[up,dn,net8angle_Up,net8angle_Down] = NAME_CLUSTERS_UP_DOWN(centroids)
YeoNetNames = {'MED FRONT', 'FPN', 'DMN', 'SUB', 'MOTOR', 'VIS I', 'VIS II', 'VIS III'};
numNets = numel(YeoNetNames);

[nparc,numClusters] = size(centroids);
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



%% once the number of clusters is determined, find the partition that best represents your clustering (i.e. is most similar to every other randomly initialized clustering)
numClusters = best_number_of_clusters;
parts = NaN(size(np,1),50);
% maybe takes a while?
for i=1:50
    [parts(:,i),~,sumd]=kmeans(np,numClusters,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',1000);
end

%% calculate adjusted mutual information for every pair of partitions
ami_results = NaN(50,50);
for i=1:50
    for j=1:50
        ami_results(i,j) = ami(parts(:,i),parts(:,j));
    end
end
ind=1
ind=best_number_of_clusters;
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

%% count number of each cluster per scan
partition=cluster_output2(:,4)
clear A

z=1;
l=1;
clear stroke_count
clear control_count
clear *count_state*
clear *FO*
for i=1:47
    nsess=5;
    if i==6
        nsess=4;
        j=5
        stroke_count_state1{i,j}=[0]
        stroke_count_state2{i,j}=[0]
        stroke_count_state3{i,j}=[0]
        stroke_count_state4{i,j}=[0]

        stroke_FO1{i,j}=[0]
        stroke_FO2{i,j}=[0]
        stroke_FO3{i,j}=[0]
        stroke_FO4{i,j}=[0]
    end
    if i==12
        nsess=3;
        for j=4:5
            stroke_count_state1{i,j}=[0]
            stroke_count_state2{i,j}=[0]
            stroke_count_state3{i,j}=[0]
            stroke_count_state4{i,j}=[0]
            stroke_FO1{i,j}=[0]
            stroke_FO2{i,j}=[0]
            stroke_FO3{i,j}=[0]
            stroke_FO4{i,j}=[0]
        end
    end

    if i==20
        nsess=2;
        for j=3:5
        stroke_count_state1{i,j}=[0]
        stroke_count_state2{i,j}=[0]
        stroke_count_state3{i,j}=[0]
        stroke_count_state4{i,j}=[0]
        stroke_FO1{i,j}=[0]
        stroke_FO2{i,j}=[0]
        stroke_FO3{i,j}=[0]
        stroke_FO4{i,j}=[0]
        end
    end
    for j=1:nsess
        clusters{i,j}=partition(z:z+leng2(l)-1);
        if i<24
            %stroke
            counts=hist(clusters{i,j}, 1:4);
            stroke_count_state1{i,j}=counts(1);
            stroke_count_state2{i,j}=counts(2);
            stroke_count_state3{i,j}=counts(3);
            stroke_count_state4{i,j}=counts(4);
            stroke_FO1{i,j}=stroke_count_state1{i,j}/clusters2(i,j);
            stroke_FO2{i,j}=stroke_count_state2{i,j}/clusters2(i,j);
            stroke_FO3{i,j}=stroke_count_state3{i,j}/clusters2(i,j);
            stroke_FO4{i,j}=stroke_count_state4{i,j}/clusters2(i,j);
            dwell=zeros(size(clusters{i,j},1),1)
            c=1;
            for k=1:4
                c=1
                dwell=zeros(size(clusters{i,j},1),1)

                timeseries=clusters{i,j};
                appear=timeseries==k
                appear=[appear;NaN]
                for p=1:size(appear,1)-1
                    if appear(p)==0 
                        continue;
                    elseif appear(p)==1 && appear(p+1)==0
                        dwell(c)=dwell(c)+1;
                        c=c+1; % if that state is over, start new count.
                        continue;
                    elseif appear(p)==1 && appear(p+1)==1
                        dwell(c)=dwell(c)+1;
                        continue;
                    elseif appear(p+1)==NaN
                    	exit;
                    end
                end
                dwell_all_stroke{i,j,k}=dwell(1:c);
                dwell_avg_stroke(i,j,k)=mean(dwell(1:c));

            end
            %FO
        else
            %control
            counts=hist(clusters{i,j}, 1:4);
            control_count_state1{i-23,j}=counts(1);
            control_count_state2{i-23,j}=counts(2);
            control_count_state3{i-23,j}=counts(3);
            control_count_state4{i-23,j}=counts(4);
            control_FO1{i-23,j}=control_count_state1{i-23,j}/clusters2(i,j);
            control_FO2{i-23,j}=control_count_state2{i-23,j}/clusters2(i,j);
            control_FO3{i-23,j}=control_count_state3{i-23,j}/clusters2(i,j);
            control_FO4{i-23,j}=control_count_state4{i-23,j}/clusters2(i,j);
            dwell=zeros(size(clusters{i,j},1),1)
            c=1;
            for k=1:4
                c=1
                dwell=zeros(size(clusters{i,j},1),1)

                timeseries=clusters{i,j};
                appear=timeseries==k
                appear=[appear;NaN]
                for p=1:size(appear,1)-1
                    if appear(p)==0 
                        continue;
                    elseif appear(p)==1 && appear(p+1)==0
                        dwell(c)=dwell(c)+1;
                        c=c+1; % if that state is over, start new count.
                        continue;
                    elseif appear(p)==1 && appear(p+1)==1
                        dwell(c)=dwell(c)+1;
                        continue;
                    elseif appear(p+1)==NaN
                    	exit;
                    end
                end
                dwell_all_control{i-23,j,k}=dwell(1:c);
                dwell_avg_control(i-23,j,k)=mean(dwell(1:c));
            end
        end
        z=z+leng2(l);
        l=l+1;
    end
end

% transition probabilities.
%GET_TRANS_PROBS(partition,subjInd,numClusters)
subjID=[];
for i=1:size(leng2,2)
    subjID = [subjID; ones(leng2(i), 1)*i];
end
[transProbs,transitionProbabilityMats,numTransitions] = GET_TRANS_PROBS(partition,subjID, 4);  

[~,nopersist_transitionProbabilityMats] = GET_TRANS_PROBS_NO_PERSIST(partition, subjID);  

nsessions=[1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 6 6 6 6 7 7 7 7 7 8 8 8 8 8 9 9 9 9 9 10 10 10 10 10 11 11 11 11 11 12 12 12 13 13 13 13 13 14 14 14 14 14 15 15 15 15 15 16 16 16 16 16 17 17 17 17 17 18 18 18 18 18 19 19 19 19 19 20 20 21 21 21 21 21 22 22 22 22 22 23 23 23 23 23];
vec=[]
for i=24:47
    vec=[vec;ones(5, 1)*i];
end
nsessions=[nsessions,vec'];

for i=1:47
    vec=nsessions==i;
    tmp=transitionProbabilityMats(vec,:,:);
    for j=1:size(tmp,1) %nsessions
        for p=1:4
            tran_prob{i,j}=squeeze(tmp(j,:,:));
        end
    end
end

for i=1:47
    vec=nsessions==i;
    tmp=nopersist_transitionProbabilityMats(vec,:,:);
    for j=1:size(tmp,1) %nsessions
        for p=1:4
            tran_prob_nopersist{i,j}=squeeze(tmp(j,:,:));
        end
    end
end

tran_prob_stroke=tran_prob(1:23,:)
tran_prob_control=tran_prob(24:47,:)

tran_prob_stroke_nopersist=tran_prob_nopersist(1:23,:)
tran_prob_control_nopersist=tran_prob_nopersist(24:47,:)

% get appearance rates
clear *appearance*
for i=1:47
    nsess=5;
    if i==20
        nsess=2
        for j=3:5
        stroke_appearance1{i,j}=[0]
        stroke_appearance2{i,j}=[0]
        stroke_appearance3{i,j}=[0]
        stroke_appearance4{i,j}=[0]
        end
    end
    if i==12
        nsess=3
        for j=4:5
        stroke_appearance1{i,j}=[0]
        stroke_appearance2{i,j}=[0]
        stroke_appearance3{i,j}=[0]
        stroke_appearance4{i,j}=[0]
        end
    end 
    if i==6
        nsess=4
        for j=5
        stroke_appearance1{i,j}=[0]
        stroke_appearance2{i,j}=[0]
        stroke_appearance3{i,j}=[0]
        stroke_appearance4{i,j}=[0]
        end
    end
    for j=1:nsess
        if i<=23
            vector=clusters{i,j};
            appear1=0;
            appear2=0;
            appear3=0;
            appear4=0;
            index = 1
            while index <= clusters2(i,j)
                disp(['cluster of TR: ', num2str(index)])
                c=1;
                if vector(index)==1
                    disp('New appearance: state 1')
                    appear1=appear1+1;
                    try 
                        while vector(index+c)==1
                            disp('Still in state 1')
                            c=c+1;
                        end
                    end
                    disp(['Leaving state 1'])
                    index=index+c;
                    continue;
                elseif vector(index)==2
                    disp('New appearance: state 2')
                    appear2=appear2+1;
                    try
                        while vector(index+c)==2
                            disp('Still in state 2')
                            c=c+1;
                        end
                    end
                    disp(['Leaving state 2'])
                    index=index+c;
                    continue;
                elseif vector(index)==3
                    disp('New appearance: state 3')
                    appear3=appear3+1;
                    try 
                        while vector(index+c)==3
                            disp('Still in state 3')
                            c=c+1;
                        end
                    end
                    index=index+c;
                    disp(['Leaving state 3'])
                    continue;
                elseif vector(index)==4
                    disp('New appearance: state 4')
                    appear4=appear4+1;
                    try
                        while vector(index+c)==4
                            disp('Still in state 4')
                            c=c+1;
                        end
                    end
                    disp(['Leaving state 4'])
                    index=index+c;
                    continue;
                end
            end
            stroke_appearance1{i,j}=appear1/(clusters2(i,j)*3/60);
            stroke_appearance2{i,j}=appear2/(clusters2(i,j)*3/60);
            stroke_appearance3{i,j}=appear3/(clusters2(i,j)*3/60);
            stroke_appearance4{i,j}=appear4/(clusters2(i,j)*3/60);
        else
            vector=clusters{i,j};
            appear1=0;
            appear2=0;
            appear3=0;
            appear4=0;
            index = 1
            while index <= clusters2(i,j)
                disp(['cluster of TR: ', num2str(index)])
                c=1;
                if vector(index)==1
                    disp('New appearance: state 1')
                    appear1=appear1+1;
                    try 
                        while vector(index+c)==1
                            disp('Still in state 1')
                            c=c+1;
                        end
                    end
                    disp(['Leaving state 1'])
                    index=index+c;
                    continue;
                elseif vector(index)==2
                    disp('New appearance: state 2')
                    appear2=appear2+1;
                    try
                        while vector(index+c)==2
                            disp('Still in state 2')
                            c=c+1;
                        end
                    end
                    disp(['Leaving state 2'])
                    index=index+c;
                    continue;
                elseif vector(index)==3
                    disp('New appearance: state 3')
                    appear3=appear3+1;
                    try 
                        while vector(index+c)==3
                            disp('Still in state 3')
                            c=c+1;
                        end
                    end
                    index=index+c;
                    disp(['Leaving state 3'])
                    continue;
                elseif vector(index)==4
                    disp('New appearance: state 4')
                    appear4=appear4+1;
                    try
                        while vector(index+c)==4
                            disp('Still in state 4')
                            c=c+1;
                        end
                    end
                    disp(['Leaving state 4'])
                    index=index+c;
                    continue;
                end
            end
            control_appearance1{i-23,j}=appear1/(clusters2(i,j)*3/60);
            control_appearance2{i-23,j}=appear2/(clusters2(i,j)*3/60);
            control_appearance3{i-23,j}=appear3/(clusters2(i,j)*3/60);
            control_appearance4{i-23,j}=appear4/(clusters2(i,j)*3/60);
        end
    end
end


% Set data to 0 when subjects missing fmri

%stroke
stroke_count_state1=cell2mat(stroke_count_state1)
stroke_count_state2=cell2mat(stroke_count_state2)
stroke_count_state3=cell2mat(stroke_count_state3)
stroke_count_state4=cell2mat(stroke_count_state4)
stroke_FO1=cell2mat(stroke_FO1)
stroke_FO2=cell2mat(stroke_FO2)
stroke_FO3=cell2mat(stroke_FO3)
stroke_FO4=cell2mat(stroke_FO4)
stroke_appearance1=cell2mat(stroke_appearance1)
stroke_appearance2=cell2mat(stroke_appearance2)
stroke_appearance3=cell2mat(stroke_appearance3)
stroke_appearance4=cell2mat(stroke_appearance4)

stroke_FO1(20,3)=NaN
stroke_FO2(20,3)=NaN
stroke_FO3(20,3)=NaN
stroke_FO4(20,3)=NaN

stroke_FO1(20,4)=NaN
stroke_FO2(20,4)=NaN
stroke_FO3(20,4)=NaN
stroke_FO4(20,4)=NaN

stroke_FO1(20,5)=NaN
stroke_FO2(20,5)=NaN
stroke_FO3(20,5)=NaN
stroke_FO4(20,5)=NaN

stroke_FO1(12,4)=NaN
stroke_FO2(12,4)=NaN
stroke_FO3(12,4)=NaN
stroke_FO4(12,4)=NaN

stroke_FO1(12,5)=NaN
stroke_FO2(12,5)=NaN
stroke_FO3(12,5)=NaN
stroke_FO4(12,5)=NaN

stroke_FO1(6,5)=NaN
stroke_FO2(6,5)=NaN
stroke_FO3(6,5)=NaN
stroke_FO4(6,5)=NaN


stroke_appearance1(20,3)=NaN
stroke_appearance2(20,3)=NaN
stroke_appearance3(20,3)=NaN
stroke_appearance4(20,3)=NaN

stroke_appearance1(20,4)=NaN
stroke_appearance2(20,4)=NaN
stroke_appearance3(20,4)=NaN
stroke_appearance4(20,4)=NaN

stroke_appearance1(20,5)=NaN
stroke_appearance2(20,5)=NaN
stroke_appearance3(20,5)=NaN
stroke_appearance4(20,5)=NaN

stroke_appearance1(12,4)=NaN
stroke_appearance2(12,4)=NaN
stroke_appearance3(12,4)=NaN
stroke_appearance4(12,4)=NaN

stroke_appearance1(12,5)=NaN
stroke_appearance2(12,5)=NaN
stroke_appearance3(12,5)=NaN
stroke_appearance4(12,5)=NaN

stroke_appearance1(6,5)=NaN
stroke_appearance2(6,5)=NaN
stroke_appearance3(6,5)=NaN
stroke_appearance4(6,5)=NaN



dwell_avg_stroke(20,5,:)=NaN
dwell_avg_stroke(20,4,:)=NaN
dwell_avg_stroke(20,3,:)=NaN
dwell_avg_stroke(12,5,:)=NaN
dwell_avg_stroke(12,4,:)=NaN
dwell_avg_stroke(12,5,:)=NaN
dwell_avg_stroke(6,5,:)=NaN

%controls
control_count_state1=cell2mat(control_count_state1)
control_count_state2=cell2mat(control_count_state2)
control_count_state3=cell2mat(control_count_state3)
control_count_state4=cell2mat(control_count_state4)
control_FO1=cell2mat(control_FO1)
control_FO2=cell2mat(control_FO2)
control_FO3=cell2mat(control_FO3)
control_FO4=cell2mat(control_FO4)
control_appearance1=cell2mat(control_appearance1)
control_appearance2=cell2mat(control_appearance2)
control_appearance3=cell2mat(control_appearance3)
control_appearance4=cell2mat(control_appearance4)

%% Analysis - across all sessions, difference between stroke & control transition probabilities

% perm test
nperms = 100000;
%pvals_twotail = PERM_TEST(tran_prob_stroke,tran_prob_control,nperms);
    % INPUT:
    % A and B: NxKxK stack of KxK transition matrices for N subjects and two conditions, A and B
    % nperms: number of permutations, pretty fast so 10000 is good to get a precise p-value
    
 
% all sessions combined
allstroke=reshape(tran_prob_stroke,1,[])
tmp=cell2mat(allstroke(i));
allstroke_trans=[]
for i=1:size(allstroke,2)
    tmp=cell2mat(allstroke(i));
    if isempty(tmp)
        continue
    end
     allstroke_trans=cat(3, allstroke_trans,tmp);
end
allstroke_trans=permute(allstroke_trans, [3 1 2])
%format control
allctl=reshape(tran_prob_control,1,[])
tmp=cell2mat(allctl(i));
allcontrol_trans=[]
for i=1:size(allctl,2)
    tmp=cell2mat(allctl(i));
    if isempty(tmp)
        continue
    end
     allcontrol_trans=cat(3, allcontrol_trans,tmp);
end
allcontrol_trans=permute(allcontrol_trans, [3 1 2])

pvals_twotail = PERM_TEST(allstroke_trans,allcontrol_trans,nperms);

s1stroke_trans=[]
s1control_trans=[]
clear p1
clear stats
close all;
clear allstroke
clear allstroke_trans
clear allcontrol_trans

% across all 5 sessions.
for a=1:4
    for b=1:4
        allstroke=reshape(tran_prob_stroke,1,[])
        for i=1:size(allstroke,2)
            disp(i)
            tmp=cell2mat(allstroke(i));
            if isempty(tmp)
                continue
            end
            allstroke_trans(i)=tmp(a,b);
        end
        allcontrol=reshape(tran_prob_control,1,[])
        for i=1:size(allcontrol,2)
            disp(i)
            tmp=cell2mat(allcontrol(i));
            allcontrol_trans(i)=tmp(a,b);
        end
        [h, p1(a,b),~, stat]=ttest2(allstroke_trans, allcontrol_trans)
        stats(a,b)=stat.tstat;
    end
end
[h, ~, ~, p_adj] = fdr_bh(pvals_twotail, 0.05,'pdep')

%plot results across all sessions
close all;
imagesc(stats)
colorbar;
title('T-statistics: Stroke - Control')
signif=pvals_twotail<0.05;
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.02, '*', 'FontSize', 30, 'Color', 'white')
end

set(gca, 'FontSize', 15)
xticks(1:4)
xticklabels({"1", "2","3","4"})
yticks(1:4)
yticklabels({"1", "2","3","4"})
caxis([-3,3])
set(gcf, 'Position', [0 0 600 500])
xlabel('Next state')
ylabel('Current state')
colormap(brewermap([],'RdBu'))



%% Transition probabilities: separate sessions
% across only session 1
allstroke=tran_prob_stroke(:,1)
allstroke_trans=[]
for i=1:size(allstroke,1)
    tmp=cell2mat(allstroke(i));
    if isempty(tmp)
        continue
    end
     allstroke_trans=cat(3, allstroke_trans,tmp)
end
allstroke_trans=permute(allstroke_trans, [3 1 2])
%format control
allctl=tran_prob_control(:,1)
tmp=cell2mat(allctl(i));
allcontrol_trans=[]
for i=1:size(allctl,1)
    tmp=cell2mat(allctl(i));
    if isempty(tmp)
        continue
    end
     allcontrol_trans=cat(3, allcontrol_trans,tmp);
end
allcontrol_trans=permute(allcontrol_trans, [3 1 2])
pvals_twotail_S1 = PERM_TEST(allstroke_trans,allcontrol_trans,nperms);
% only across session 2
allstroke=tran_prob_stroke(:,2)
allstroke_trans=[]
for i=1:size(allstroke,1)
    tmp=cell2mat(allstroke(i));
    if isempty(tmp)
        continue
    end
     allstroke_trans=cat(3, allstroke_trans,tmp)
end
allstroke_trans=permute(allstroke_trans, [3 1 2])
%format control
allctl=tran_prob_control(:,2)
tmp=cell2mat(allctl(i));
allcontrol_trans=[]
for i=1:size(allctl,1)
    tmp=cell2mat(allctl(i));
    if isempty(tmp)
        continue
    end
     allcontrol_trans=cat(3, allcontrol_trans,tmp);
end
allcontrol_trans=permute(allcontrol_trans, [3 1 2])
pvals_twotail_S2 = PERM_TEST(allstroke_trans,allcontrol_trans,nperms);
% session 3
allstroke=tran_prob_stroke(:,3)
allstroke_trans=[]
for i=1:size(allstroke,1)
    tmp=cell2mat(allstroke(i));
    if isempty(tmp)
        continue
    end
     allstroke_trans=cat(3, allstroke_trans,tmp)
end
allstroke_trans=permute(allstroke_trans, [3 1 2])
%format control
allctl=tran_prob_control(:,3)
tmp=cell2mat(allctl(i));
allcontrol_trans=[]
for i=1:size(allctl,1)
    tmp=cell2mat(allctl(i));
    if isempty(tmp)
        continue
    end
     allcontrol_trans=cat(3, allcontrol_trans,tmp);
end
allcontrol_trans=permute(allcontrol_trans, [3 1 2])
pvals_twotail_S3 = PERM_TEST(allstroke_trans,allcontrol_trans,nperms);
% session 4
allstroke=tran_prob_stroke(:,4)
allstroke_trans=[]
for i=1:size(allstroke,1)
    tmp=cell2mat(allstroke(i));
    if isempty(tmp)
        continue
    end
     allstroke_trans=cat(3, allstroke_trans,tmp)
end
allstroke_trans=permute(allstroke_trans, [3 1 2])
%format control
allctl=tran_prob_control(:,4)
tmp=cell2mat(allctl(i));
allcontrol_trans=[]
for i=1:size(allctl,1)
    tmp=cell2mat(allctl(i));
    if isempty(tmp)
        continue
    end
     allcontrol_trans=cat(3, allcontrol_trans,tmp);
end
allcontrol_trans=permute(allcontrol_trans, [3 1 2])
pvals_twotail_S4 = PERM_TEST(allstroke_trans,allcontrol_trans,nperms);
% session 5
allstroke=tran_prob_stroke(:,5)
allstroke_trans=[]
for i=1:size(allstroke,1)
    tmp=cell2mat(allstroke(i));
    if isempty(tmp)
        continue
    end
     allstroke_trans=cat(3, allstroke_trans,tmp)
end
allstroke_trans=permute(allstroke_trans, [3 1 2])
%format control
allctl=tran_prob_control(:,5)
tmp=cell2mat(allctl(i));
allcontrol_trans=[]
for i=1:size(allctl,1)
    tmp=cell2mat(allctl(i));
    if isempty(tmp)
        continue
    end
     allcontrol_trans=cat(3, allcontrol_trans,tmp);
end
allcontrol_trans=permute(allcontrol_trans, [3 1 2])
pvals_twotail_S5 = PERM_TEST(allstroke_trans,allcontrol_trans,nperms);


% plot results of permutation pvalue testing
% (need to run code below first to get t-stats)
tiledlayout(5,1,'padding', 'none')
% session 1
nexttile;imagesc(stats1); colorbar;title('Session 1'); 
signif=pvals_twotail_S1<0.05;
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.02, '*', 'FontSize', 30, 'Color', 'white')
end
caxis([-3, 3])
yticks(1:4)
ylabel('Current state')
xlabel('Next state')

% session 2
nexttile;imagesc(stats2); colorbar;title('Session 2'); 
signif=pvals_twotail_S2<0.05;
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.02, '*', 'FontSize', 30, 'Color', 'white')
end
caxis([-3, 3])
yticks(1:4)
ylabel('Current state')
xlabel('Next state')

% session 3
nexttile;imagesc(stats3); colorbar;title('Session 3'); 
signif=pvals_twotail_S3<0.05;
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.02, '*', 'FontSize', 30, 'Color', 'white')
end
caxis([-3, 3])
yticks(1:4)
ylabel('Current state')
xlabel('Next state')

% session 4
nexttile;imagesc(stats4); colorbar;title('Session 4'); 
signif=pvals_twotail_S4<0.05;
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.02, '*', 'FontSize', 30, 'Color', 'white')
end
caxis([-3, 3])
yticks(1:4)
ylabel('Current state')
xlabel('Next state')

% session 5
nexttile;imagesc(stats5); colorbar;title('Session 5'); 
signif=pvals_twotail_S5<0.05;
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.02, '*', 'FontSize', 30, 'Color', 'white')
end
caxis([-3, 3])
yticks(1:4)
ylabel('Current state')
xlabel('Next state')
colormap(brewermap([],'RdBu'))

pall=[pvals_twotail_S1,pvals_twotail_S2,pvals_twotail_S3,pvals_twotail_S4,pvals_twotail_S5]
[h, ~, ~, p_adj] = fdr_bh(pvals_twotail_S1, 0.05,'pdep')


clear p1
clear stats
close all;

%%  get t-stats for transition probs for separate sessions

%test significance between transition probabilities for each session
%separately
s1stroke_trans=[]
s1control_trans=[]
% session 1
for a=1:4
    for b=1:4
        s1stroke=tran_prob_stroke(:,1)
        for i=1:size(s1stroke,1)
            tmp=cell2mat(s1stroke(i));
            s1stroke_trans(i)=tmp(a,b);
        end
        s1control=tran_prob_control(:,1)
        for i=1:size(s1control,1)
            tmp=cell2mat(s1control(i));
            s1control_trans(i)=tmp(a,b);
        end
        [h, p1(a,b),~, stat]=ttest2(s1stroke_trans, s1control_trans)
        stats1(a,b)=stat.tstat;

    end
end
[h, ~, ~, p_adj] = fdr_bh(p1, 0.05,'pdep')

%correlate w recovery/impairment
s1stroke=tran_prob_stroke(:,1)
alltransprob=[]
for a=1:4
    for b=1:4
        for i=1:size(s1stroke,1)
            tmp=cell2mat(s1stroke(i));
            alltransprob(i)=tmp(a,b);
        end
        [rho,pvals(a,b)]=corr(fm_5-fm_1,alltransprob', 'rows', 'complete')
    end
end

[rho, p]=corr(fm_1,s1stroke_trans', 'rows', 'complete')

% session 2
s1stroke_trans=[]
s1control_trans=[]
clear p2
for a=1:4
    for b=1:4
        s1stroke=tran_prob_stroke(:,2)
        for i=1:size(s1stroke,1)
            tmp=cell2mat(s1stroke(i));
            s1stroke_trans(i)=tmp(a,b);
        end
         s1control=tran_prob_control(:,2)
        for i=1:size(s1control,1)
            tmp=cell2mat(s1control(i));
            s1control_trans(i)=tmp(a,b);
        end
        [h, p2(a,b),~, stat]=ttest2(s1stroke_trans, s1control_trans)
         stats2(a,b)=stat.tstat;

    end
end
[h, ~, ~, p_adj] = fdr_bh(p2, 0.05,'pdep')

% session 3
s1stroke_trans=[]
s1control_trans=[]
clear p3
for a=1:4
    for b=1:4
        s1stroke=tran_prob_stroke(:,3)
        for i=1:size(s1stroke,1)
            if i==20
                continue
            end
            tmp=cell2mat(s1stroke(i));
            s1stroke_trans(i)=tmp(a,b);
        end
         s1control=tran_prob_control(:,3)
        for i=1:size(s1control,1)
            tmp=cell2mat(s1control(i));
            s1control_trans(i)=tmp(a,b);
        end
        [h, p3(a,b),~, stat]=ttest2(s1stroke_trans, s1control_trans)
         stats3(a,b)=stat.tstat;

    end
end
[h, ~, ~, p_adj] = fdr_bh(p3, 0.05,'pdep')

% session 4
s1stroke_trans=[]
s1control_trans=[]
clear p4
for a=1:4
    for b=1:4
        s1stroke=tran_prob_stroke(:,4)
        for i=1:size(s1stroke,1)
            if i==20
                continue;
            elseif i==12
                continue;
            end
            tmp=cell2mat(s1stroke(i));
            s1stroke_trans(i)=tmp(a,b);
        end
         s1control=tran_prob_control(:,4)
        for i=1:size(s1control,1)
            tmp=cell2mat(s1control(i));
            s1control_trans(i)=tmp(a,b);
        end
        [h, p4(a,b),~, stat]=ttest2(s1stroke_trans, s1control_trans)
         stats4(a,b)=stat.tstat;

    end
end
[h, ~, ~, p_adj] = fdr_bh(p4, 0.05,'pdep')

% session 5
s1stroke_trans=[]
s1control_trans=[]
clear p5
for a=1:4
    for b=1:4
        s1stroke=tran_prob_stroke(:,5)
        for i=1:size(s1stroke,1)
             if i==20
                continue;
            elseif i==12
                continue;
             elseif i==6
                 continue;
            end
            tmp=cell2mat(s1stroke(i));
            s1stroke_trans(i)=tmp(a,b);
        end
        s1control=tran_prob_control(:,5)
        for i=1:size(s1control,1)
            tmp=cell2mat(s1control(i));
            s1control_trans(i)=tmp(a,b);
        end
        [h, p5(a,b),~, stat]=ttest2(s1stroke_trans, s1control_trans)
        stats5(a,b)=stat.tstat;
    end
end
[h, ~, ~, p_adj] = fdr_bh(p5, 0.05,'pdep')

% visualize results
tiledlayout(1,5)
nexttile;imagesc(p1); caxis([0 0.05]);nexttile;imagesc(p2);caxis([0 0.05]);nexttile;imagesc(p3);caxis([0 0.05]);nexttile;imagesc(p4);caxis([0 0.05]);nexttile;imagesc(p5);caxis([0 0.05]);
colormap 'hot'

%stroke 1 vs stroke 5
s1stroke_trans=[]
for a=1:4
    for b=1:4
        s1stroke=tran_prob_stroke(:,1)
        for i=1:size(s1stroke,1)
             if i==20
                continue;
            elseif i==12
                continue;
             elseif i==6
                 continue;
            end
            tmp=cell2mat(s1stroke(i));
            s1stroke_trans(i)=tmp(a,b);
        end
        s5stroke=tran_prob_stroke(:,5)
        for i=1:size(s5stroke,1)
             if i==20
                continue;
            elseif i==12
                continue;
             elseif i==6
                 continue;
            end
            tmp=cell2mat(s5stroke(i));
            s5stroke_trans(i)=tmp(a,b);
        end
        [h, p5(a,b),~, ~]=ttest(s1stroke_trans, s5stroke_trans)
    end
end
imagesc(p5)

[h, ~, ~, p_adj] = fdr_bh(p5, 0.05,'pdep')

imagesc(p_adj)


%% Analysis - Fractional occupancy and dwell time between stroke & control
% FO differences across all sessions, stroke and controls;
figure('Position', [0 0 900 300])
idx=[zeros(115, 1); ones(120,1); ...
    ones(115, 1)*2; ones(120,1)*3;...
    ones(115, 1)*4; ones(120,1)*5;...
    ones(115, 1)*6; ones(120,1)*7
    ]

violin=violinplot([reshape(stroke_FO1, 1,[]),reshape(control_FO1, 1,[]),...
    reshape(stroke_FO2, 1,[]),reshape(control_FO2, 1,[]),...
    reshape(stroke_FO3, 1,[]),reshape(control_FO3, 1,[]),...
    reshape(stroke_FO4, 1,[]),reshape(control_FO4, 1,[])
    ], idx)
ylabel('Fractional Occupancy')
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]
violin(3).ViolinColor=[0.4 0.4 0.8]
violin(4).ViolinColor=[1 0.4 0.4]
violin(5).ViolinColor=[0.4 0.4 0.8]
violin(6).ViolinColor=[1 0.4 0.4]
violin(7).ViolinColor=[0.4 0.4 0.8]
violin(8).ViolinColor=[1 0.4 0.4]

xticks(1:2:12)
xticklabels({'State 1','State 2','State 3', 'State 4'})
set(gca, 'FontSize', 15)

[h, p(1), ci, stats]=ttest2(reshape(stroke_FO1, 1,[]),reshape(control_FO1, 1,[]))
[h, p(2), ci, stats]=ttest2(reshape(stroke_FO2, 1,[]),reshape(control_FO2, 1,[]))
[h, p(3), ci, stats]=ttest2(reshape(stroke_FO3, 1,[]),reshape(control_FO3, 1,[]))
[h, p(4), ci, stats]=ttest2(reshape(stroke_FO4, 1,[]),reshape(control_FO4, 1,[]))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')


%% Dwell time differences across all sessions, stroke and controls;
figure('Position', [0 0 900 300])
idx=[zeros(115, 1); ones(120,1); ...
    ones(115, 1)*2; ones(120,1)*3;...
    ones(115, 1)*4; ones(120,1)*5;...
    ones(115, 1)*6; ones(120,1)*7
    ]

violin=violinplot([reshape(dwell_avg_stroke(:,:,1), 1,[]),reshape(dwell_avg_control(:,:,1), 1,[]),...
    reshape(dwell_avg_stroke(:,:,2), 1,[]),reshape(dwell_avg_control(:,:,2), 1,[]),...
    reshape(dwell_avg_stroke(:,:,3), 1,[]),reshape(dwell_avg_control(:,:,3), 1,[]),...
    reshape(dwell_avg_stroke(:,:,4), 1,[]),reshape(dwell_avg_control(:,:,4), 1,[])
    ], idx)

violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]
violin(3).ViolinColor=[0.4 0.4 0.8]
violin(4).ViolinColor=[1 0.4 0.4]
violin(5).ViolinColor=[0.4 0.4 0.8]
violin(6).ViolinColor=[1 0.4 0.4]
violin(7).ViolinColor=[0.4 0.4 0.8]
violin(8).ViolinColor=[1 0.4 0.4]
xticks(1:2:12)
xticklabels({'State 1','State 2','State 3', 'State 4'})
set(gca, 'FontSize', 15)
ylabel('Dwell time (TRs)')

[h, p(3), ci, stats]=ttest2(reshape(dwell_avg_stroke(:,:,3), 1,[]),reshape(dwell_avg_control(:,:,3), 1,[]))
[h, p(4), ci, stats]=ttest2(reshape(dwell_avg_stroke(:,:,4), 1,[]),reshape(dwell_avg_control(:,:,4), 1,[]))
[h, p(1), ci, stats]=ttest2(reshape(dwell_avg_stroke(:,:,1), 1,[]),reshape(dwell_avg_control(:,:,1), 1,[]))
[h, p(2), ci, stats]=ttest2(reshape(dwell_avg_stroke(:,:,2), 1,[]),reshape(dwell_avg_control(:,:,2), 1,[]))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')


%% Appearance rate
%  differences across all sessions, stroke and controls;
figure('Position', [0 0 900 300])
idx=[zeros(115, 1); ones(120,1); ...
    ones(115, 1)*2; ones(120,1)*3;...
    ones(115, 1)*4; ones(120,1)*5;...
    ones(115, 1)*6; ones(120,1)*7
    ]

violin=violinplot([reshape(stroke_appearance1, 1,[]),reshape(control_appearance1, 1,[]),...
    reshape(stroke_appearance2, 1,[]),reshape(control_appearance2, 1,[]),...
    reshape(stroke_appearance3, 1,[]),reshape(control_appearance3, 1,[]),...
    reshape(stroke_appearance4, 1,[]),reshape(control_appearance4, 1,[])], idx)

test=[reshape(stroke_appearance1, 1,[]),reshape(control_appearance1, 1,[])]
ylabel('Appearance Rate')
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]
violin(3).ViolinColor=[0.4 0.4 0.8]
violin(4).ViolinColor=[1 0.4 0.4]
violin(5).ViolinColor=[0.4 0.4 0.8]
violin(6).ViolinColor=[1 0.4 0.4]
violin(7).ViolinColor=[0.4 0.4 0.8]
violin(8).ViolinColor=[1 0.4 0.4]

xticks(1:2:12)
xticklabels({'State 1','State 2','State 3', 'State 4'})
set(gca, 'FontSize', 15)

[h, p(1), ci, stats]=ttest2(reshape(stroke_appearance1, 1,[]),reshape(control_appearance1, 1,[]))
[h, p(2), ci, stats]=ttest2(reshape(stroke_appearance2, 1,[]),reshape(control_appearance2, 1,[]))
[h, p(3), ci, stats]=ttest2(reshape(stroke_appearance3, 1,[]),reshape(control_appearance3, 1,[]))
[h, p(4), ci, stats]=ttest2(reshape(stroke_appearance4, 1,[]),reshape(control_appearance4, 1,[]))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

%% Severe vs moderate stroke subjects
idx_severe=baselineFM > 50
idx_mod=baselineFM <= 50
sev=sum(idx_severe)
mod=sum(idx_mod)


%% Fractional occupancy - moderate to severe
% across all sessions
severe_FO1=stroke_FO1(idx_severe,:)
mod_FO1=stroke_FO1(idx_mod,:)
severe_FO2=stroke_FO2(idx_severe,:)
mod_FO2=stroke_FO2(idx_mod,:)
severe_FO3=stroke_FO3(idx_severe,:)
mod_FO3=stroke_FO3(idx_mod,:)
severe_FO4=stroke_FO4(idx_severe,:)
mod_FO4=stroke_FO4(idx_mod,:)

ids=[zeros(1,50), ones(1,65)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([reshape(severe_FO1(:,:),[],1);reshape(mod_FO1(:,:), [],1)],ids)
nexttile;
violinplot([reshape(severe_FO2(:,:),[],1);reshape(mod_FO2(:,:), [],1)],ids)
nexttile;
violinplot([reshape(severe_FO3(:,:),[],1);reshape(mod_FO3(:,:), [],1)],ids)
nexttile;
violinplot([reshape(severe_FO4(:,:),[],1);reshape(mod_FO4(:,:), [],1)],ids)
set(gcf, 'Position', [0 0 1000 1000])

[h, p(1), ~, ~]=ttest2(reshape(severe_FO1(:,:),[],1), reshape(mod_FO1(:,:), [],1))
[h, p(2), ~, ~]=ttest2(reshape(severe_FO2(:,:),[],1), reshape(mod_FO2(:,:), [],1))
[h, p(3), ~, ~]=ttest2(reshape(severe_FO3(:,:),[],1), reshape(mod_FO3(:,:), [],1))
[h, p(4), ~, ~]=ttest2(reshape(severe_FO4(:,:),[],1), reshape(mod_FO4(:,:), [],1))



%% Dwell time - moderate vs severe
% across single sessions
severe_FO1=dwell_avg_stroke(idx_severe,:,1)
mod_FO1=dwell_avg_stroke(idx_mod,:,1)
severe_FO2=dwell_avg_stroke(idx_severe,:,2)
mod_FO2=dwell_avg_stroke(idx_mod,:,2)
severe_FO3=dwell_avg_stroke(idx_severe,:,3)
mod_FO3=dwell_avg_stroke(idx_mod,:,3)
severe_FO4=dwell_avg_stroke(idx_severe,:,4)
mod_FO4=dwell_avg_stroke(idx_mod,:,4)

%across all sessions
ids=[zeros(1,50), ones(1,65)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([reshape(severe_FO1(:,:),[],1);reshape(mod_FO1(:,:), [],1)],ids)
nexttile;
violinplot([reshape(severe_FO2(:,:),[],1);reshape(mod_FO2(:,:), [],1)],ids)
nexttile;
violinplot([reshape(severe_FO3(:,:),[],1);reshape(mod_FO3(:,:), [],1)],ids)
nexttile;
violinplot([reshape(severe_FO4(:,:),[],1);reshape(mod_FO4(:,:), [],1)],ids)
set(gcf, 'Position', [0 0 1000 1000])

[h, p(1), ~, ~]=ttest2(reshape(severe_FO1(:,:),[],1), reshape(mod_FO1(:,:), [],1))
[h, p(2), ~, ~]=ttest2(reshape(severe_FO2(:,:),[],1), reshape(mod_FO2(:,:), [],1))
[h, p(3), ~, ~]=ttest2(reshape(severe_FO3(:,:),[],1), reshape(mod_FO3(:,:), [],1))
[h, p(4), ~, ~]=ttest2(reshape(severe_FO4(:,:),[],1), reshape(mod_FO4(:,:), [],1))

%% Appearance rate - moderate vs severe
% across single sessions
severe_FO1=stroke_appearance1(idx_severe,:)
mod_FO1=stroke_appearance1(idx_mod,:)
severe_FO2=stroke_appearance2(idx_severe,:)
mod_FO2=stroke_appearance2(idx_mod,:)
severe_FO3=stroke_appearance3(idx_severe,:)
mod_FO3=stroke_appearance3(idx_mod,:)
severe_FO4=stroke_appearance4(idx_severe,:)
mod_FO4=stroke_appearance4(idx_mod,:)

ids=[zeros(1,sev), ones(1,mod)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([severe_FO1(:,1);mod_FO1(:,1)],ids)
nexttile;
violinplot([severe_FO2(:,1);mod_FO2(:,1)],ids)
nexttile;
violinplot([severe_FO3(:,1);mod_FO3(:,1)],ids)
nexttile;
violinplot([severe_FO4(:,1);mod_FO4(:,1)],ids)
set(gcf, 'Position', [0 0 1000 1000])

[h, p, ~, ~]=ttest2(severe_FO1(:,1), mod_FO1(:,1))
[h, p, ~, ~]=ttest2(severe_FO2(:,1), mod_FO2(:,1))
[h, p, ~, ~]=ttest2(severe_FO3(:,1), mod_FO3(:,1))
[h, p, ~, ~]=ttest2(severe_FO4(:,1), mod_FO4(:,1))

% across all sessions
ids=[zeros(1,50), ones(1,65)]
tiledlayout(1,4,'padding', 'none')
nexttile;
viol=violinplot([reshape(severe_FO1(:,:),[],1);reshape(mod_FO1(:,:), [],1)],ids)
viol(1).ViolinColor=[0.4 0.4 0.8]
viol(2).ViolinColor=[1 0.4 0.4]
ylim([1, 5])
title('State 1')
xticklabels({"Severe", "Moderate"})
nexttile;

viol=violinplot([reshape(severe_FO2(:,:),[],1);reshape(mod_FO2(:,:), [],1)],ids)
viol(1).ViolinColor=[0.4 0.4 0.8]
viol(2).ViolinColor=[1 0.4 0.4]
ylim([1, 5])
title('State 2')
xticklabels({"Severe", "Moderate"})

nexttile;

viol=violinplot([reshape(severe_FO3(:,:),[],1);reshape(mod_FO3(:,:), [],1)],ids)
viol(1).ViolinColor=[0.4 0.4 0.8]
viol(2).ViolinColor=[1 0.4 0.4]
ylim([1, 5])
title('State 3')
xticklabels({"Severe", "Moderate"})

nexttile;

viol=violinplot([reshape(severe_FO4(:,:),[],1);reshape(mod_FO4(:,:), [],1)],ids)
viol(1).ViolinColor=[0.4 0.4 0.8]
viol(2).ViolinColor=[1 0.4 0.4]
ylim([1, 5])
title('State 4')
xticklabels({"Severe", "Moderate"})

set(gcf, 'Position', [0 0 1000 1000])

[h, p(1), ~, ~]=ttest2(reshape(severe_FO1(:,:),[],1), reshape(mod_FO1(:,:), [],1))
[h, p(2), ~, ~]=ttest2(reshape(severe_FO2(:,:),[],1), reshape(mod_FO2(:,:), [],1))
[h, p(3), ~, ~]=ttest2(reshape(severe_FO3(:,:),[],1), reshape(mod_FO3(:,:), [],1))
[h, p(4), ~, ~]=ttest2(reshape(severe_FO4(:,:),[],1), reshape(mod_FO4(:,:), [],1))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

%% Fractional occupancy - session-specific comparison
ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_FO1(:,1);control_FO1(:,1)],ids)
nexttile;
violinplot([stroke_FO2(:,1);control_FO2(:,1)],ids)
nexttile;
violinplot([stroke_FO3(:,1);control_FO3(:,1)],ids)
nexttile;
violinplot([stroke_FO4(:,1);control_FO4(:,1)],ids)

[h,p(1),~,~]=ttest2(stroke_FO1(:,1),control_FO1(:,1))
[h,p(2),~,~]=ttest2(stroke_FO2(:,1),control_FO2(:,1))
[h,p(3),~,~]=ttest2(stroke_FO3(:,1),control_FO3(:,1))
[h,p(4),~,~]=ttest2(stroke_FO4(:,1),control_FO4(:,1))

% plot session w sig diffs
close all;
figure()
ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violin=violinplot([stroke_FO1(:,2);control_FO1(:,2)],ids)
title('State 1')
xticklabels({"Stroke", "Control"})
set(gca, 'FontSize', 15)
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]
nexttile;
violin=violinplot([stroke_FO2(:,2);control_FO2(:,2)],ids)
title('State 2')
xticklabels({"Stroke", "Control"})
set(gca, 'FontSize', 15)
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]

nexttile;
violin=violinplot([stroke_FO3(:,2);control_FO3(:,2)],ids)
title('State 3')
xticklabels({"Stroke", "Control"})
set(gca, 'FontSize', 15)
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]

nexttile;
violin=violinplot([stroke_FO4(:,2);control_FO4(:,2)],ids)
title('State 4')
xticklabels({"Stroke", "Control"})
sgtitle('Fractial occupancy - 2 weeks post-stroke')
set(gca, 'FontSize', 15)
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]

[h,p(5),~,~]=ttest2(stroke_FO1(:,2),control_FO1(:,2))
[h,p(6),~,~]=ttest2(stroke_FO2(:,2),control_FO2(:,2))
[h,p(7),~,~]=ttest2(stroke_FO3(:,2),control_FO3(:,2))
[h,p(8),~,~]=ttest2(stroke_FO4(:,2),control_FO4(:,2))


figure()
ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_FO1(:,3);control_FO1(:,3)],ids)
nexttile;
violinplot([stroke_FO2(:,3);control_FO2(:,3)],ids)
nexttile;
violinplot([stroke_FO3(:,3);control_FO3(:,3)],ids)
nexttile;
violinplot([stroke_FO4(:,3);control_FO4(:,3)],ids)
[h,p(9),~,~]=ttest2(stroke_FO1(:,3),control_FO1(:,3))
[h,p(10),~,~]=ttest2(stroke_FO2(:,3),control_FO2(:,3))
[h,p(11),~,~]=ttest2(stroke_FO3(:,3),control_FO3(:,3))
[h,p(12),~,~]=ttest2(stroke_FO4(:,3),control_FO4(:,3))



ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_FO1(:,4);control_FO1(:,4)],ids)
nexttile;
violinplot([stroke_FO2(:,4);control_FO2(:,4)],ids)
nexttile;
violinplot([stroke_FO3(:,4);control_FO3(:,4)],ids)
nexttile;
violinplot([stroke_FO4(:,4);control_FO4(:,4)],ids)

[h,p(13),~,~]=ttest2(stroke_FO1(:,4),control_FO1(:,4))
[h,p(14),~,~]=ttest2(stroke_FO2(:,4),control_FO2(:,4))
[h,p(15),~,~]=ttest2(stroke_FO3(:,4),control_FO3(:,4))
[h,p(16),~,~]=ttest2(stroke_FO4(:,4),control_FO4(:,4))


ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_FO1(:,5);control_FO1(:,5)],ids)
nexttile;
violinplot([stroke_FO2(:,5);control_FO2(:,5)],ids)
nexttile;
violinplot([stroke_FO3(:,5);control_FO3(:,5)],ids)
nexttile;
violinplot([stroke_FO4(:,5);control_FO4(:,5)],ids)

[h,p(17),~,~]=ttest2(stroke_FO1(:,5),control_FO1(:,5))
[h,p(18),~,~]=ttest2(stroke_FO2(:,5),control_FO2(:,5))
[h,p(19),~,~]=ttest2(stroke_FO3(:,5),control_FO3(:,5))
[h,p(20),~,~]=ttest2(stroke_FO4(:,5),control_FO4(:,5))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')


%% Dwell time - specific sessions
% stroke vs control

idx=[ones(23,1);zeros(24,1)]
% session 1 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,1,i);dwell_avg_control(:,1,i)], idx)
    [h, p(i), ci, stats]=ttest2(dwell_avg_stroke(:,1,i),dwell_avg_control(:,1,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end

% session 2 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,2,i);dwell_avg_control(:,2,i)], idx)
    [h, p(i+4), ci, stats]=ttest2(dwell_avg_stroke(:,2,i),dwell_avg_control(:,2,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end

% session 3 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,3,i);dwell_avg_control(:,3,i)], idx)
    [h, p(i+8), ci, stats]=ttest2(dwell_avg_stroke(:,3,i),dwell_avg_control(:,3,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end

% session 4 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,4,i);dwell_avg_control(:,4,i)], idx)
    [h, p(i+12), ci, stats]=ttest2(dwell_avg_stroke(:,4,i),dwell_avg_control(:,4,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end

% session 5 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,5,i);dwell_avg_control(:,5,i)], idx)
    [h, p(i+16), ci, stats]=ttest2(dwell_avg_stroke(:,5,i),dwell_avg_control(:,5,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

%% Trajectory of FO in state2
idx=[ones(23,1); ones(23,1)*2;ones(23,1)*3;ones(23,1)*4;ones(23,1)*5];

viol=violinplot([stroke_FO2(:,1);stroke_FO2(:,2);stroke_FO2(:,3);stroke_FO2(:,4);stroke_FO2(:,5)], idx)
title('State 2')
ylim([0, 0.6])
viol(1).ViolinColor=[0.5 0.5 0.8]
viol(2).ViolinColor=[0.35 0.35 0.8]
viol(3).ViolinColor=[0.30 0.30 1]
viol(4).ViolinColor=[0.25 0.25 1]
ylabel('Fractional Occupancy')

viol=violinplot([stroke_FO3(:,1);stroke_FO3(:,2);stroke_FO3(:,3);stroke_FO3(:,4);stroke_FO3(:,5)], idx)
title('State 3')
ylim([0, 0.6])
viol(1).ViolinColor=[0.5 0.5 0.8]
viol(2).ViolinColor=[0.35 0.35 0.8]
viol(3).ViolinColor=[0.30 0.30 1]
viol(4).ViolinColor=[0.25 0.25 1]
plot([stroke_FO3(:,1),stroke_FO3(:,2),stroke_FO3(:,3),stroke_FO3(:,4),stroke_FO3(:,5)]', 'LineWidth', 3)
alldata=mean([stroke_FO3(:,1),stroke_FO3(:,2),stroke_FO3(:,3),stroke_FO3(:,4),stroke_FO3(:,5)]',2, 'omitnan')
hold on;
plot(alldata, '-r', 'LineWidth', 4)
ylabel('Fractional Occupancy')

%% relationship to impairment and recovery

% load motor scores
fm_dir=strcat('/Users/emilyolafson/GIT/stroke-graph-matching/data/');
fuglmeyer=readtable(strcat(fm_dir, 'fuglmeyer_allpts.csv'));
fm_1=fuglmeyer.Var2;
fm_2=fuglmeyer.Var3;
fm_3=fuglmeyer.Var4;
fm_4=fuglmeyer.Var5;
fm_5=fuglmeyer.Var6;

fm_1(23)=NaN;
fm_1(22)=NaN;
fm_3(20)=NaN;
fm_4(12)=NaN;
fm_4(20)=NaN;
fm_5(20)=NaN;
fm_5(6)=NaN;

%% Fractional Occupancy- related to impairment at same scan

[rho, p(1)]=corr(stroke_FO1(:,1), fm_1, 'rows', 'complete')
[rho, p(2)]=corr(stroke_FO2(:,1), fm_1, 'rows', 'complete')
[rho, p(3)]=corr(stroke_FO3(:,1), fm_1, 'rows', 'complete')
[rho, p(4)]=corr(stroke_FO4(:,1), fm_1, 'rows', 'complete')

[rho, p(5)]=corr(stroke_FO1(:,2), fm_2, 'rows', 'complete')
[rho, p(6)]=corr(stroke_FO2(:,2), fm_2, 'rows', 'complete')
[rho, p(7)]=corr(stroke_FO3(:,2), fm_2, 'rows', 'complete')
[rho, p(8)]=corr(stroke_FO4(:,2), fm_2, 'rows', 'complete')

[rho, p(9)]=corr(stroke_FO1(:,3), fm_3, 'rows', 'complete')
[rho, p(10)]=corr(stroke_FO2(:,3), fm_3, 'rows', 'complete')
[rho, p(11)]=corr(stroke_FO3(:,3), fm_3, 'rows', 'complete')
[rho, p(12)]=corr(stroke_FO4(:,3), fm_3, 'rows', 'complete')

[rho, p(13)]=corr(stroke_FO1(:,4), fm_4, 'rows', 'complete')
[rho, p(14)]=corr(stroke_FO2(:,4), fm_4, 'rows', 'complete')
[rho, p(15)]=corr(stroke_FO3(:,4), fm_4, 'rows', 'complete')
[rho, p(16)]=corr(stroke_FO4(:,4), fm_4, 'rows', 'complete')

[rho, p(17)]=corr(stroke_FO1(:,5), fm_5, 'rows', 'complete')
[rho, p(18)]=corr(stroke_FO2(:,5), fm_5, 'rows', 'complete')
[rho, p(19)]=corr(stroke_FO3(:,5), fm_5, 'rows', 'complete')
[rho, p(20)]=corr(stroke_FO4(:,5), fm_5, 'rows', 'complete')

% diff in parameters related to change in fm score?
fm12=fm_2-fm_1;
fm23=fm_3-fm_2
fm34=fm_4-fm_3;
fm45=fm_5-fm_4;

fo_12=stroke_FO1(:,2)-stroke_FO1(:,1);
fo_23=stroke_FO1(:,3)-stroke_FO1(:,2);
fo_34=stroke_FO1(:,4)-stroke_FO1(:,3);
fo_45=stroke_FO1(:,5)-stroke_FO1(:,4);

[rho,p(1)]=corr(fm12, fo_12, 'rows', 'complete')
[rho,p(2)]=corr(fm23, fo_23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, fo_34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, fo_45, 'rows', 'complete')

%%  Dwell time - related to impairment at same scan
[rho, p(1)]=corr(dwell_avg_stroke(:,1,1), fm_1, 'rows', 'complete')
[rho, p(2)]=corr(dwell_avg_stroke(:,1,2), fm_1, 'rows', 'complete')
[rho, p(3)]=corr(dwell_avg_stroke(:,1,3), fm_1, 'rows', 'complete')
[rho, p(4)]=corr(dwell_avg_stroke(:,1,4), fm_1, 'rows', 'complete')

[rho, p(5)]=corr(dwell_avg_stroke(:,2,1), fm_2, 'rows', 'complete')
[rho, p(6)]=corr(dwell_avg_stroke(:,2,2), fm_2, 'rows', 'complete')
[rho, p(7)]=corr(dwell_avg_stroke(:,2,3), fm_2, 'rows', 'complete')
[rho, p(8)]=corr(dwell_avg_stroke(:,2,4), fm_2, 'rows', 'complete')

[rho, p(9)]=corr(dwell_avg_stroke(:,3,1), fm_3, 'rows', 'complete')
[rho, p(10)]=corr(dwell_avg_stroke(:,3,2), fm_3, 'rows', 'complete')
[rho, p(11)]=corr(dwell_avg_stroke(:,3,3), fm_3, 'rows', 'complete')
[rho, p(12)]=corr(dwell_avg_stroke(:,3,4), fm_3, 'rows', 'complete')

[rho, p(13)]=corr(dwell_avg_stroke(:,4,1), fm_4, 'rows', 'complete')
[rho, p(14)]=corr(dwell_avg_stroke(:,4,2), fm_4, 'rows', 'complete')
[rho, p(15)]=corr(dwell_avg_stroke(:,4,3), fm_4, 'rows', 'complete')
[rho, p(16)]=corr(dwell_avg_stroke(:,4,4), fm_4, 'rows', 'complete')

[rho, p(17)]=corr(dwell_avg_stroke(:,5,1), fm_5, 'rows', 'complete')
[rho, p(18)]=corr(dwell_avg_stroke(:,5,2), fm_5, 'rows', 'complete')
[rho, p(19)]=corr(dwell_avg_stroke(:,5,3), fm_5, 'rows', 'complete')
[rho, p(20)]=corr(dwell_avg_stroke(:,5,4), fm_5, 'rows', 'complete')

%state 1
dt12=dwell_avg_stroke(:,2,1)-dwell_avg_stroke(:,1,1)
dt23=dwell_avg_stroke(:,3,1)-dwell_avg_stroke(:,2,1)
dt34=dwell_avg_stroke(:,4,1)-dwell_avg_stroke(:,3,1)
dt45=dwell_avg_stroke(:,5,1)-dwell_avg_stroke(:,4,1)

clear p
[rho,p(1)]=corr(fm12, dt12, 'rows', 'complete')
[rho,p(2)]=corr(fm23, dt23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, dt34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, dt45, 'rows', 'complete')

% state 2
dt12=dwell_avg_stroke(:,2,2)-dwell_avg_stroke(:,1,2)
dt23=dwell_avg_stroke(:,3,2)-dwell_avg_stroke(:,2,2)
dt34=dwell_avg_stroke(:,4,2)-dwell_avg_stroke(:,3,2)
dt45=dwell_avg_stroke(:,5,2)-dwell_avg_stroke(:,4,2)

clear p
[rho,p(1)]=corr(fm12, dt12, 'rows', 'complete')
[rho,p(2)]=corr(fm23, dt23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, dt34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, dt45, 'rows', 'complete')

% state 3
dt12=dwell_avg_stroke(:,2,3)-dwell_avg_stroke(:,1,3)
dt23=dwell_avg_stroke(:,3,3)-dwell_avg_stroke(:,2,3)
dt34=dwell_avg_stroke(:,4,3)-dwell_avg_stroke(:,3,3)
dt45=dwell_avg_stroke(:,5,3)-dwell_avg_stroke(:,4,3)

clear p
[rho,p(1)]=corr(fm12, dt12, 'rows', 'complete', 'Type', 'Spearman')
[rho,p(2)]=corr(fm23, dt23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, dt34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, dt45, 'rows', 'complete')

close all
figure()
scatter(fm12, dt12, 'ok', 'filled')
xlabel('Change in Fugl-Meyer scores 1-2 weeks post-stroke')
ylabel('Change in State 3 dwell time 1-2 weeks post-stroke')
text(30, -0.4, 'Corr = 0.53, p = 0.013 (unc)', 'FontSize', 14)
set(gca, 'FontSize', 14)
lsline
title('State 3 change in dwell time vs. change in FM')
[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

% state 4
dt12=dwell_avg_stroke(:,2,4)-dwell_avg_stroke(:,1,4)
dt23=dwell_avg_stroke(:,3,4)-dwell_avg_stroke(:,2,4)
dt34=dwell_avg_stroke(:,4,4)-dwell_avg_stroke(:,3,4)
dt45=dwell_avg_stroke(:,5,4)-dwell_avg_stroke(:,4,4)

clear p
[rho,p(1)]=corr(fm12, dt12, 'rows', 'complete')
[rho,p(2)]=corr(fm23, dt23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, dt34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, dt45, 'rows', 'complete')

%% Appearance rate
clear p
[rho, p(1)]=corr(stroke_appearance1(:,1), fm_1, 'rows', 'complete')
[rho, p(2)]=corr(stroke_appearance2(:,1), fm_1, 'rows', 'complete')
[rho, p(3)]=corr(stroke_appearance3(:,1), fm_1, 'rows', 'complete')
[rho, p(4)]=corr(stroke_appearance4(:,1), fm_1, 'rows', 'complete')

[rho, p(5)]=corr(stroke_appearance1(:,2), fm_2, 'rows', 'complete')
[rho, p(6)]=corr(stroke_appearance2(:,2), fm_2, 'rows', 'complete')
[rho, p(7)]=corr(stroke_appearance3(:,2), fm_2, 'rows', 'complete')
[rho, p(8)]=corr(stroke_appearance4(:,2), fm_2, 'rows', 'complete')

[rho, p(9)]=corr(stroke_appearance1(:,3), fm_3, 'rows', 'complete')
[rho, p(10)]=corr(stroke_appearance2(:,3), fm_3, 'rows', 'complete')
[rho, p(11)]=corr(stroke_appearance3(:,3), fm_3, 'rows', 'complete')
[rho, p(12)]=corr(stroke_appearance4(:,3), fm_3, 'rows', 'complete')

[rho, p(13)]=corr(stroke_appearance1(:,4), fm_4, 'rows', 'complete')
[rho, p(14)]=corr(stroke_appearance2(:,4), fm_4, 'rows', 'complete')
[rho, p(15)]=corr(stroke_appearance3(:,4), fm_4, 'rows', 'complete')
[rho, p(16)]=corr(stroke_appearance4(:,4), fm_4, 'rows', 'complete')

[rho, p(17)]=corr(stroke_appearance1(:,5), fm_5, 'rows', 'complete')
[rho, p(18)]=corr(stroke_appearance2(:,5), fm_5, 'rows', 'complete')
[rho, p(19)]=corr(stroke_appearance3(:,5), fm_5, 'rows', 'complete')
[rho, p(20)]=corr(stroke_appearance4(:,5), fm_5, 'rows', 'complete')

% diff in parameters related to change in fm score?
fm12=fm_2-fm_1;
fm23=fm_3-fm_2
fm34=fm_4-fm_3;
fm45=fm_5-fm_4;

clear p
fo_12=stroke_appearance1(:,2)-stroke_appearance1(:,1);
fo_23=stroke_appearance1(:,3)-stroke_appearance1(:,2);
fo_34=stroke_appearance1(:,4)-stroke_appearance1(:,3);
fo_45=stroke_appearance1(:,5)-stroke_appearance1(:,4);

[rho,p(1)]=corr(fm12, fo_12, 'rows', 'complete')
[rho,p(2)]=corr(fm23, fo_23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, fo_34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, fo_45, 'rows', 'complete')

clear p
fo_12=stroke_appearance2(:,2)-stroke_appearance2(:,1);
fo_23=stroke_appearance2(:,3)-stroke_appearance2(:,2);
fo_34=stroke_appearance2(:,4)-stroke_appearance2(:,3);
fo_45=stroke_appearance2(:,5)-stroke_appearance2(:,4);

[rho,p(1)]=corr(fm12, fo_12, 'rows', 'complete')
[rho,p(2)]=corr(fm23, fo_23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, fo_34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, fo_45, 'rows', 'complete')


clear p
fo_12=stroke_appearance3(:,2)-stroke_appearance3(:,1);
fo_23=stroke_appearance3(:,3)-stroke_appearance3(:,2);
fo_34=stroke_appearance3(:,4)-stroke_appearance3(:,3);
fo_45=stroke_appearance3(:,5)-stroke_appearance3(:,4);

[rho,p(1)]=corr(fm12, fo_12, 'rows', 'complete')
[rho,p(2)]=corr(fm23, fo_23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, fo_34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, fo_45, 'rows', 'complete')




% difference relative to controls correlated with recovery/impairment?
avgctl=mean(control_FO1(:,:),2)
meanctl=mean(avgctl)
stdctl=std(avgctl,1)
zscoretst1=(stroke_FO1(:,1)-meanctl*ones(23,1))./stdctl

avgctl=mean(control_FO2(:,:),2)
meanctl=mean(avgctl)
stdctl=std(avgctl,1)
zscoretst2=(stroke_FO2(:,1)-meanctl*ones(23,1))./stdctl

avgctl=mean(control_FO3(:,:),2)
meanctl=mean(avgctl)
stdctl=std(avgctl,1)
zscoretst3=(stroke_FO3(:,1)-meanctl*ones(23,1))./stdctl

avgctl=mean(control_FO4(:,:),2)
meanctl=mean(avgctl)
stdctl=std(avgctl,1)
zscoretst4=(stroke_FO4(:,1)-meanctl*ones(23,1))./stdctl

[rho,p]=corr(zscoretst1, fm_5, 'rows', 'complete')
[rho,p]=corr(zscoretst2, fm_5, 'rows', 'complete')
[rho,p]=corr(zscoretst3, fm_5, 'rows', 'complete')
[rho,p]=corr(zscoretst4, fm_5, 'rows', 'complete')

% dwell time

% difference relative to controls correlated with recovery/impairment?
state=1
avgctl=squeeze(mean(dwell_avg_control(:,:,:),2))
meanctl=mean(avgctl(:,state))
stdctl=std(avgctl(:,state),1)
zscoretst1=(dwell_avg_stroke(:,1,state)-meanctl*ones(23,1))./stdctl

state=2
avgctl=squeeze(mean(dwell_avg_control(:,:,:),2))
meanctl=mean(avgctl(:,state))
stdctl=std(avgctl(:,state),1)
zscoretst2=(dwell_avg_stroke(:,1,state)-meanctl*ones(23,1))./stdctl

state=3
avgctl=squeeze(mean(dwell_avg_control(:,:,:),2))
meanctl=mean(avgctl(:,state))
stdctl=std(avgctl(:,state),1)
zscoretst3=(dwell_avg_stroke(:,1,state)-meanctl*ones(23,1))./stdctl

state=4
avgctl=squeeze(mean(dwell_avg_control(:,:,:),2))
meanctl=mean(avgctl(:,state))
stdctl=std(avgctl(:,state),1)
zscoretst4=(dwell_avg_stroke(:,1,state)-meanctl*ones(23,1))./stdctl

[rho,p]=corr(zscoretst1, fm_5, 'rows', 'complete')
[rho,p]=corr(zscoretst2, fm_5, 'rows', 'complete')
[rho,p]=corr(zscoretst3, fm_5, 'rows', 'complete')
[rho,p]=corr(zscoretst4, fm_5, 'rows', 'complete')





