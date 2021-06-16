
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
for i=4:5
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


%% once the number of clusters is determined, find the partition that best represents your clustering (i.e. is most similar to every other randomly initialized clustering)
numClusters = 4;
parts = NaN(size(np,1),50);
sumd=[]
% maybe takes a while?
for i=1:50
    disp(i)
    [parts(:,i),~,sumd]=kmeans(np,numClusters,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',1000);
    disp(sumd)
end
save('/Users/emilyolafson/GIT/dynamic-brainstates/results/cluster_output_k4k5.mat', 'cluster_output2')
save('/Users/emilyolafson/GIT/dynamic-brainstates/results/sum_squared_distk4k5.mat', 'totSum')
save('/Users/emilyolafson/GIT/dynamic-brainstates/results/partitions_k4_50reps.mat', 'parts')

%% calculate adjusted mutual information for every pair of partitions
ami_results = NaN(50,50);
for i=1:50
    for j=1:50
        ami_results(i,j) = ami(parts(:,i),parts(:,j));
    end
end
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

%% Visualize cluster centroids
% compute centroids and plot
best_number_of_clusters =4 %define
centroids = GET_CENTROIDS(np,partition,best_number_of_clusters);
atlasblobs_list=load('gummibrain/atlasblobs_saved.mat');
atlasblobs_list=atlasblobs_list.atlasblobs_list;

whichatlas={'shen268'};

% gummibrain code (need to install from kjamison repo) https://github.com/kjamison/atlasblobs
results_dir='/Users/emilyolafson/GIT/dynamic-brainstates/results/'
for i=1:best_number_of_clusters
    close all;
    %cc400_data needs to be a 1x392 vector
    %shen368 needs to be 1x268
    data = centroids(:,i);
    data_min=-0.8;
    data_max=0.8;
    img=display_atlas_blobs(data,atlasblobs_list,...
        'atlasname',whichatlas,...
        'render',true,...
        'backgroundimage',false,...
        'crop',true,...
        'colormap',cmap,...
        'alpha', data);
    %rescale(mean_fi(2,:).^2)
    %'alpha',rescale(pinv(rescale(abs(cc400_fi')))).^2
    %'alpha',cc400_t_facealpha
    %'alpha', rescale(abs(mean_fs86_size_featimp))
    figure('Position', [0 0 2000 2000]);
    imshow(img);
    c=colorbar('SouthOutside', 'fontsize', 20);
    c.Label.String='Colorbar Label';
    set(gca,'colormap',cmap);
    caxis([data_min data_max]);
    %title('TT', 'fontsize',18);
   % annotation('textbox',[.835 .38 .1 .2],'String','RH','EdgeColor','none', 'fontsize',20,'color','white')
   % annotation('textbox',[.12 .38 .1 .2],'String','LH','EdgeColor','none', 'fontsize',20,'color','white')
   % annotation('textbox',[.45 .9 .1 .04],'String','Lateral','EdgeColor','none', 'fontsize',20,'color','white', 'horizontalalignment','center','backgroundcolor','black')
   % annotation('textbox',[.45 .21 .1 .04],'String','Medial','EdgeColor','none', 'fontsize',20,'color','white', 'horizontalalignment','center', 'backgroundcolor','black')
    set(gcf, 'Position', [0 0 2000 2000]);
    saveas(gcf, strcat(results_dir, '/state', num2str(i), '.png'))
    pause(2)
end

  
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
f=figure('Position', [0 0 5000 5000]);
for K = 1:numClusters
    ax = subplot(1,numClusters,K,polaraxes); hold on
    polarplot(netAngle,[net8angle_Up(K,:) net8angle_Up(K,1)],'k', 'LineWidth', 3);
    polarplot(netAngle,[net8angle_Down(K,:) net8angle_Down(K,1)],'r', 'LineWidth', 3);
    thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
    rticks([0.4 0.8]); rticklabels({'0.4','0.8'});rlim([0 0.6])
%     for L = 1:numNets
%         ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
%         YeoColors(L,:), ax.ThetaTickLabel{L});
%     end
    set(ax,'FontSize',15);
    %title(overallNames{K},'Color',clusterColors(K,:),'FontSize',8);
end


%% count number of each cluster per scan
%partition=cluster_output2(:,4)
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

