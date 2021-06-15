
baselineFM=load('baselineFM.mat');
baselineFM=baselineFM.basline;
save('baseline_scores.mat', 'baselineFM')

numscans1_11=[4 3 3 3 3 3 3 3 3 3 3;4 3 3 3 3 3 3 4 3 3 3;3 3 3 3 3 3 3 3 3 3 3;4 4 3 3 3 3 3 3 3 3 3 ; 3 3 3 3 3 0 3 3 3 3 3];
numscans12_23=[3 3 3 3 3 2 2 2 2 2 2 2; 3 3 3 4 3 2 2 2 2 2 2 2 ;3 3 3 3 3 2 2 2 0 2 2 2;0 3 3 3 3 2 2 2 0 2 2 2;0 3 4 3 3 2 2 2 0 2 2 2];
nscans = [numscans1_11, numscans12_23];
nscan = ones(5,24)+4;
nscans = [nscans, nscan];

%%
allsubs2=load('ts_GSR_shen268_allsub_allsess.mat', 'allsubs')
allsubs2=allsubs2.allsubs;
allsubs_fmri=reshape(allsubs2', [1 235]); %reshape goes column-wise; all S1 listed, then S2,..

% Determine the optimal number of clusters using variance explained
allsubs_fmri=allsubs_fmri(find(~cellfun(@isempty,allsubs_fmri)))
allsubs_fmri=allsubs_fmri';

np= [];
for i=1:229
    a=allsubs_fmri{i};
    b=cell2mat(a);
    leng(i)=size(b,1);
    %b=b(4:173,:);
    np=[np;b];
end
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
        clusters(i,j)=size(partition(c:c+leng(l)-1),2)
        c=c+leng(i);
        l=l+1;
    end
end

%load in functional mri scans from all 5 sessions
allsubs=load('data/ts_shenGSR_may29.mat'); 
allsubs=allsubs.newaverg;
allsubs_fmri=reshape(allsubs', [1 235]); %reshape goes column-wise; all S1 listed, then S2,..

% Determine the optimal number of clusters using variance explained
allsubs_fmri=allsubs_fmri(find(~cellfun(@isempty,allsubs_fmri)))
allsubs_fmri=allsubs_fmri';
np= [];
for i=1:229
    a=allsubs_fmri{i};
    b=cell2mat(a);
    leng2(i)=size(b,1);
    %b=b(4:173,:);
    np=[np;b];
end

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
distanceMethod= 'sqeuclidean'
nreps = 50;
totSum=[];

% should take about 5 mins to run for each cluster #
clear cluster_output
clear sumd
for i=6
    [cluster_output(:,i),~,sumd]=kmeans(np,i,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',1000);
    totSum(i)=sum(sumd);
end

% Elbow plots
for i=1:5
  disp(num2str(i))
  [~,~,sumd]=kmeans(np,i); %sumd= within-cluster sums of point-to-centroid distances
  totSum(i)=sum(sumd); % Inertia
  avgDist(i)=mean(sumd); % Distortion
end

best_number_of_clusters = 6 %define

%% Visualize cluster centroids
% compute centroids and plot
centroids = GET_CENTROIDS(np,cluster_output(:,best_number_of_clusters),best_number_of_clusters);
save('/Volumes/MyPassport/centroids.mat', 'centroids')
% name clusters based on alignment with Yeo resting state networks
clusterNames = [1:6];

f = figure;
imagesc(centroids); title('Centroids'); xticks(1:best_number_of_clusters); xticklabels(clusterNames);
colormap('parula'); axis square; colorbar; set(gca,'FontSize',12); COLOR_TICK_LABELS(true,false,best_number_of_clusters);

xticklabels(clusterNames); yticklabels(clusterNames); xtickangle(90);
COLOR_TICK_LABELS(true,true,best_number_of_clusters);
f.PaperUnits = 'inches';
f.PaperSize = [10 10];
f.PaperPosition = [0 0 4 2];

%% visualize clusters on brain
visualise_kmeans(centroids)

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
partition=cluster_output(:,4)
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
[transitionProbability,transitionProbabilityMats,numTransitions] = GET_TRANS_PROBS(partition,subjID, 4);  

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

tran_prob_stroke=tran_prob(1:23,:)
tran_prob_control=tran_prob(24:47,:)

%test significance between transition probabilities
s1stroke_trans=[]
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
        [h, p1(a,b),~, ~]=ttest2(s1stroke_trans, s1control_trans)
    end
end
[h, ~, ~, p_adj] = fdr_bh(p1, 0.05,'pdep')


% session 2
s1stroke_trans=[]
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
        [h, p2(a,b),~, ~]=ttest2(s1stroke_trans, s1control_trans)
        c=c+1;
    end
end
[h, ~, ~, p_adj] = fdr_bh(p2, 0.05,'pdep')

% session 3
s1stroke_trans=[]
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
        [h, p3(a,b),~, ~]=ttest2(s1stroke_trans, s1control_trans)
        
    end
end
[h, ~, ~, p_adj] = fdr_bh(p3, 0.05,'pdep')

% session 4
s1stroke_trans=[]
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
        [h, p4(a,b),~, ~]=ttest2(s1stroke_trans, s1control_trans)
    end
end
[h, ~, ~, p_adj] = fdr_bh(p4, 0.05,'pdep')

% session 5
s1stroke_trans=[]
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
        [h, p5(a,b),~, ~]=ttest2(s1stroke_trans, s1control_trans)
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


[h, ~, ~, p_adj] = fdr_bh(p5, 0.05,'pdep')

imagesc(p_adj)

%stroke
stroke_count_state1=cell2mat(stroke_count_state1)
stroke_count_state2=cell2mat(stroke_count_state2)
stroke_count_state3=cell2mat(stroke_count_state3)
stroke_count_state4=cell2mat(stroke_count_state4)
stroke_FO1=cell2mat(stroke_FO1)
stroke_FO2=cell2mat(stroke_FO2)
stroke_FO3=cell2mat(stroke_FO3)
stroke_FO4=cell2mat(stroke_FO4)

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

ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_FO1(:,2);control_FO1(:,1)],ids)
nexttile;
violinplot([stroke_FO2(:,2);control_FO2(:,1)],ids)
nexttile;
violinplot([stroke_FO3(:,2);control_FO3(:,1)],ids)
nexttile;
violinplot([stroke_FO4(:,2);control_FO4(:,1)],ids)

ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_FO1(:,3);control_FO1(:,1)],ids)
nexttile;
violinplot([stroke_FO2(:,3);control_FO2(:,1)],ids)
nexttile;
violinplot([stroke_FO3(:,3);control_FO3(:,1)],ids)
nexttile;
violinplot([stroke_FO4(:,3);control_FO4(:,1)],ids)

ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_FO1(:,4);control_FO1(:,1)],ids)
nexttile;
violinplot([stroke_FO2(:,4);control_FO2(:,1)],ids)
nexttile;
violinplot([stroke_FO3(:,4);control_FO3(:,1)],ids)
nexttile;
violinplot([stroke_FO4(:,4);control_FO4(:,1)],ids)



%% across all stroke and controls;
idx=[zeros(115, 1); ones(120,1)]
tiledlayout(2,2,'padding','none')
nexttile;
violinplot([reshape(stroke_FO1, 1,[]),reshape(control_FO1, 1,[])], idx)
nexttile;
violinplot([reshape(stroke_FO2, 1,[]),reshape(control_FO2, 1,[])], idx)
nexttile;
violinplot([reshape(stroke_FO3, 1,[]),reshape(control_FO3, 1,[])], idx)
nexttile;
violinplot([reshape(stroke_FO4, 1,[]),reshape(control_FO4, 1,[])], idx)


%%

controlavg1=mean(control_FO1,2)
controlavg2=mean(control_FO2,2)
controlavg3=mean(control_FO3,2)
controlavg4=mean(control_FO4,2)

tiledlayout(4,1,'padding','none')
nexttile;
idx=[ones(23,1); ones(23,1)*2;ones(23,1)*3;ones(23,1)*4;ones(23,1)*5;ones(24,1)*6];
violinplot([stroke_FO1(:,1);stroke_FO1(:,2);stroke_FO1(:,3);stroke_FO1(:,4);stroke_FO1(:,5); controlavg1], idx)
title('State 1')
ylim([0, 0.6])
nexttile;
violinplot([stroke_FO2(:,1);stroke_FO2(:,2);stroke_FO2(:,3);stroke_FO2(:,4);stroke_FO2(:,5); controlavg2], idx)
title('State 2')
ylim([0, 0.6])
nexttile;
violinplot([stroke_FO3(:,1);stroke_FO3(:,2);stroke_FO3(:,3);stroke_FO3(:,4);stroke_FO3(:,5); controlavg3], idx)
title('State 3')
ylim([0, 0.6])
nexttile;
violinplot([stroke_FO4(:,1);stroke_FO4(:,2);stroke_FO4(:,3);stroke_FO4(:,4);stroke_FO4(:,5); controlavg4], idx)
title('State 4')
ylim([0, 0.6])

% stroke at session 1 vs control
% stroke at session 5 vs control
[~, p_one(1), ~, stats]=ttest2(stroke_FO1(:,1), control_FO1(:,1))
[~,  p_five(1), ~, stats]=ttest2(stroke_FO1(:,5), control_FO1(:,5))

[~,  p_one(2), ~, stats]=ttest2(stroke_FO2(:,1), control_FO2(:,1))
[~,  p_five(2), ~, stats]=ttest2(stroke_FO2(:,5), control_FO2(:,5))

[~, p_one(3), ~, stats]=ttest2(stroke_FO3(:,1), control_FO3(:,1))
[~,  p_five(3), ~, stats]=ttest2(stroke_FO3(:,5), control_FO3(:,5))

[~,  p_one(4), ~, stats]=ttest2(stroke_FO4(:,1), control_FO4(:,1))
[~, p_five(4), ~, stats]=ttest2(stroke_FO4(:,5), control_FO4(:,5))

% stroke at session 1 vs stroke at session 5
% stroke at session 5 vs control
[~, p_one(1), ~, stats]=ttest2(stroke_FO1(:,1), control_FO1(:,1))
[~,  p_five(1), ~, stats]=ttest2(stroke_FO1(:,5), control_FO1(:,5))

[~,  p_one(2), ~, stats]=ttest2(stroke_FO2(:,1), control_FO2(:,5))
[~,  p_five(2), ~, stats]=ttest2(stroke_FO2(:,5), control_FO2(:,5))

[~, p_one(3), ~, stats]=ttest2(stroke_FO3(:,1), control_FO3(:,1))
[~,  p_five(3), ~, stats]=ttest2(stroke_FO3(:,5), control_FO3(:,5))

[~,  p_one(4), ~, stats]=ttest2(stroke_FO4(:,1), control_FO4(:,1))
[~, p_five(4), ~, stats]=ttest2(stroke_FO4(:,5), control_FO4(:,5))


%% Stroke session 1 vs. stroke session 5 vs. Control avg state
tiledlayout(4,1,'padding','none')
nexttile;
idx=[ones(23,1);ones(23,1)*5;ones(24,1)*6];
boxplot([stroke_FO1(:,1);stroke_FO1(:,5); controlavg1], idx)
xticklabels({'Stroke S1', 'Stroke S5', 'Control'})
title('State 1')
ylim([0, 0.6])
nexttile;
boxplot([stroke_FO2(:,1);stroke_FO2(:,5); controlavg2], idx)
title('State 2')
xticklabels({'Stroke S1', 'Stroke S5', 'Control'})
ylim([0, 0.6])
nexttile;
boxplot([stroke_FO3(:,1);stroke_FO3(:,5); controlavg3], idx)
title('State 3')
xticklabels({'Stroke S1', 'Stroke S5', 'Control'})
ylim([0, 0.6])
nexttile;
boxplot([stroke_FO4(:,1);stroke_FO4(:,5); controlavg4], idx)
title('State 4')
xticklabels({'Stroke S1', 'Stroke S5', 'Control'})
ylim([0, 0.6])

%% state 1 vs state 5
tiledlayout(2,2,'padding','none')
nexttile;
violinplot([stroke_FO1(:,1),stroke_FO1(:,5)])
[h, p(1), ci, stats]=ttest(stroke_FO1(:,1),stroke_FO1(:,5))

nexttile;
violinplot([stroke_FO2(:,1),stroke_FO2(:,5)])
[h, p(2), ci, stats]=ttest(stroke_FO2(:,1),stroke_FO2(:,5))

nexttile;
violinplot([stroke_FO3(:,1),stroke_FO3(:,5)])
[h, p(3), ci, stats]=ttest(stroke_FO3(:,1),stroke_FO3(:,5))

nexttile;
violinplot([stroke_FO4(:,1),stroke_FO4(:,5)])
[h, p(4), ci, stats]=ttest(stroke_FO4(:,1),stroke_FO4(:,5))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

%% dwell time

% stroke vs control

idx=[ones(23,1);zeros(24,1)]
% session 1 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,1,i);dwell_avg_control(:,1,i)], idx)
    [h, p, ci, stats]=ttest2(dwell_avg_stroke(:,1,i),dwell_avg_control(:,1,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end

% session 2 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,2,i);dwell_avg_control(:,2,i)], idx)
    [h, p, ci, stats]=ttest2(dwell_avg_stroke(:,2,i),dwell_avg_control(:,2,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end

% session 3 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,3,i);dwell_avg_control(:,3,i)], idx)
    [h, p, ci, stats]=ttest2(dwell_avg_stroke(:,3,i),dwell_avg_control(:,3,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end

% session 4 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,4,i);dwell_avg_control(:,4,i)], idx)
    [h, p, ci, stats]=ttest2(dwell_avg_stroke(:,4,i),dwell_avg_control(:,4,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end


% session 5 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,5,i);dwell_avg_control(:,5,i)], idx)
    [h, p, ci, stats]=ttest2(dwell_avg_stroke(:,5,i),dwell_avg_control(:,5,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end

%% 
tiledlayout(4,5,'padding','none')
for i=1:4
    for j= 1:5
    nexttile;
    idx=[zeros(23,1);ones(24,1)];
    boxplot([dwell_avg_stroke(:,j,i); dwell_avg_control(:,j,i)], idx)
    ylim([0, 5])
    title(['state: ', num2str(i)])
    end
end

%% 
tiledlayout(4,1,'padding','none')
for i=1:4
    nexttile;
    violinplot([dwell_avg_control(:,1,i),dwell_avg_control(:,2,i),dwell_avg_control(:,3,i),dwell_avg_control(:,4,i),dwell_avg_control(:,5,i)])
    title(['state: ', num2str(i)])
end

% calculate control avg dwell time across sessions
for i=1:4
    dwell_control_sessavg{i}=mean(dwell_avg_control(:,:,i),2);
end

idx=[ones(23,1); ones(23,1)*2;ones(23,1)*3;ones(23,1)*4;ones(23,1)*5;ones(24,1)*6];
tiledlayout(4,1,'padding','none')
for i=1:4
    nexttile;
    violinplot([dwell_avg_stroke(:,1,i);dwell_avg_stroke(:,2,i);dwell_avg_stroke(:,3,i);dwell_avg_stroke(:,4,i);dwell_avg_stroke(:,5,i); dwell_control_sessavg{i}], idx)
    title(['state: ', num2str(i)])
    ylim([0 6])
end

%% stroke state 1 vs state 5 - dwell time
tiledlayout(2,2,'padding','none')
nexttile;
violinplot([dwell_avg_stroke(:,1,1),dwell_avg_stroke(:,5,1)])
[h, p(1), ci, stats]=ttest(dwell_avg_stroke(:,1,1),dwell_avg_stroke(:,5,1))

nexttile;
violinplot([dwell_avg_stroke(:,1,2),dwell_avg_stroke(:,5,2)])
[h, p(2), ci, stats]=ttest(dwell_avg_stroke(:,1,2),dwell_avg_stroke(:,5,2))

nexttile;
violinplot([dwell_avg_stroke(:,1,3),dwell_avg_stroke(:,5,3)])
[h, p(3), ci, stats]=ttest(dwell_avg_stroke(:,1,3),dwell_avg_stroke(:,5,3))

nexttile;
violinplot([dwell_avg_stroke(:,1,4),dwell_avg_stroke(:,5,4)])
[h, p(4), ci, stats]=ttest(dwell_avg_stroke(:,1,4),dwell_avg_stroke(:,5,4))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')



%% plot single subjects?
violinplot([dwell_avg_stroke(:,1,1),dwell_avg_stroke(:,5,1)])
hold on;
plot([dwell_avg_stroke(:,1,1)';dwell_avg_stroke(:,5,1)'])


%amount of reduction over 5 months in state 1 corr with recovery?
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

fm15=fm_5-fm_1

diff_s1s5_state1=dwell_avg_stroke(:,1,1)-dwell_avg_stroke(:,5,1)
diff_s1s5_state4=dwell_avg_stroke(:,1,4)-dwell_avg_stroke(:,5,4)

diff_s1s5_state1=stroke_FO1(:,1)-stroke_FO1(:,5)
diff_s1s5_state4=stroke_FO4(:,1)-stroke_FO4(:,5)

plot(fm15,diff_s1s5_state1, '*')

plot(fm15,diff_s1s5_state4, '*')

[rho, p] =corr(fm15,diff_s1s5_state1, 'rows', 'complete')
[rho, p] =corr(fm15,diff_s1s5_state4, 'rows', 'complete')

%% stroke sess 1 vs sess 5 vs control

% stroke sess 1 vs control
[h, p1c(1), ci, stats]=ttest2(dwell_avg_stroke(:,1,1),dwell_avg_control(:,1,1))
[h, p1c(2), ci, stats]=ttest2(dwell_avg_stroke(:,1,2),dwell_avg_control(:,1,2))
[h, p1c(3), ci, stats]=ttest2(dwell_avg_stroke(:,1,3),dwell_avg_control(:,1,3))
[h, p1c(4), ci, stats]=ttest2(dwell_avg_stroke(:,1,4),dwell_avg_control(:,1,4))

% stroke sess 5 vs control
[h, p5c(1), ci, stats]=ttest2(dwell_avg_stroke(:,5,1),dwell_avg_control(:,5,1))
[h, p5c(2), ci, stats]=ttest2(dwell_avg_stroke(:,5,2),dwell_avg_control(:,5,2))
[h, p5c(3), ci, stats]=ttest2(dwell_avg_stroke(:,5,3),dwell_avg_control(:,5,3))
[h, p5c(4), ci, stats]=ttest2(dwell_avg_stroke(:,5,4),dwell_avg_control(:,5,4))

[h, ~, ~, p_adj] = fdr_bh(p5c, 0.05,'pdep')

