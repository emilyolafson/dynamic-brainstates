clear all;

baselineFM=load('baselineFM.mat');
baselineFM=baselineFM.basline;
save('baseline_scores.mat', 'baselineFM')

%load in functional mri scans from all 5 sessions
allsubs=load('ts_GSR_shen268_allsub_allsess.mat'); 
allsubs=allsubs.allsubs;

% control subjects have lower range? is it processing or biology?
% ctl=cell2mat(allsubs{33,1});
% stroke=cell2mat(allsubs{2,1});

allsubs_fmri=reshape(allsubs', [1 235]); % shorten scans to 12 minutes each 


%% K-means clustering
allsubs_fmri=allsubs_fmri(find(~cellfun(@isempty,allsubs_fmri)))
allsubs_fmri=allsubs_fmri';
np= [];

% concatenate timeseries into one NxP matrix; P = rois, N =
% subAtimseries,subBsimeries, etc concatenated so they are each 235 TR's
% (12 minutes)

for i=1:229 %including controls
    a=allsubs_fmri{i};
    b=cell2mat(a);
    %b=normalize(b);
    leng(i)=size(b,1);
    b=b(1:173,:); %truncated to shortest scan length (173 TRs)
    np=[np;b];
end

np_controlstroke=np;
distanceMethod= 'sqeuclidean'
nreps = 50;
totSum=[];
% should take about 5 mins to run for each cluster #
clear cluster_output
clear sumd
for i=4
    [cluster_output(:,i),~,sumd]=kmeans(np_controlstroke,i,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',1000);
    totSum(i)=sum(sumd);
end

best_number_of_clusters = 4;

%% elbow plots

nClusters=4; % pick/set number of clusters we're going to test
totSum=zeros(nClusters,1);  % preallocate the results
avgDist=zeros(nClusters,1); % preallocate the results

for i=1:6
  disp(num2str(i))
  [~,~,sumd]=kmeans(np,i); %sumd= within-cluster sums of point-to-centroid distances
  totSum(i)=sum(sumd); % Inertia
  avgDist(i)=mean(sumd); % Distortion
end

plot(avgDist)
title('Average within-cluster sums of point-to-centroid distances')
xlabel('Number of clusters')
set(gca, 'FontSize', 14)

%% plot
% compute centroids and plot
best_number_of_clusters =4;

centroids = GET_CENTROIDS(np,cluster_output(:,best_number_of_clusters),best_number_of_clusters);
min(centroids(:,1))
max(centroids(:,1))
min(centroids(:,2))
max(centroids(:,2))
min(centroids(:,3))
max(centroids(:,3))
min(centroids(:,4))
max(centroids(:,4))
centroids_norm = GET_CENTROIDS(np,cluster_output_norm(:,best_number_of_clusters),best_number_of_clusters);

% name clusters based on alignment with Yeo resting state networks
clusterNames = [1:best_number_of_clusters];

f = figure;
imagesc(centroids); title('Centroids'); xticks(1:best_number_of_clusters); xticklabels(clusterNames);
colormap('parula'); axis square; colorbar; set(gca,'FontSize',12); COLOR_TICK_LABELS(true,false,best_number_of_clusters);

xticklabels(clusterNames);  xtickangle(90);
COLOR_TICK_LABELS(true,true,best_number_of_clusters);
f.PaperUnits = 'inches';
f.PaperSize = [10 10];
f.PaperPosition = [0 0 4 2];

output=corr(centroids, centroids_norm)
imagesc(output)
%% once the number of clusters is determined, find the partition that best represents your clustering (i.e. is most similar to every other randomly initialized clustering)
numClusters = best_number_of_clusters;

numClusters=4;
parts4 = NaN(size(np,1),50);
% maybe takes a while?
for i=1:50
    [parts4(:,i),~,sumd]=kmeans(np,numClusters,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',1000);
end

numClusters=5;
parts5 = NaN(size(np,1),50);
for i=1:50
    [parts5(:,i),~,sumd]=kmeans(np,numClusters,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',1000);
end
%% calculate adjusted mutual information for every pair of partitions
% Running AMI for k=4
ami_results4 = NaN(50,50);
for i=1:50
    for j=1:50
        ami_results4(i,j) = ami(parts4(:,i),parts4(:,j));
    end
end

[m,ind] = max(sum(ami_results4,1)); %ind corresponds to the partition which has the highest mutual information with all other partitions
partition = parts4(:,ind); % take partition that has most agreement with all other for further analysis

%running AMI for k=5 

ami_results5 = NaN(50,50);
for i=1:50
    for j=1:50
        ami_results5(i,j) = ami(parts5(:,i),parts5(:,j));
    end
end
[m,ind] = max(sum(ami_results5,1)); %ind corresponds to the partition which has the highest mutual information with all other partitions
partition = parts5(:,ind); % take partition that has most agreement with all other for further analysis


% plot

%figure for k=4
f = figure;
imagesc(ami_results4); title('Adjusted Mutal Information between Partitions'); colorbar;
axis square; set(gca,'FontSize',8); 
f.PaperUnits = 'inches';
f.PaperSize = [4 2];
f.PaperPosition = [0 0 4 2];

saveas(f,fullfile(savedir,['AMI_k',num2str(4),'.pdf']));

%figure for k=5
f = figure;
imagesc(ami_results5); title('Adjusted Mutal Information between Partitions'); colorbar;
axis square; set(gca,'FontSize',8); 
f.PaperUnits = 'inches';
f.PaperSize = [4 2];
f.PaperPosition = [0 0 4 2];

saveas(f,fullfile(savedir,['AMI_k',num2str(5),'.pdf']));


%% count number of each cluster per scan
k=4;
numClusters=k;
partition = cluster_output(:,k); 
A=reshape(partition,173,[]);
indices=[1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3  1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2  1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5];
subjects=[1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 6 6 6 6 7 7 7 7 7 8 8 8 8 8 9 9 9 9 9 10 10 10 10 10 11 11 11 11 11 12 12 12 13 13 13 13 13 14 14 14 14 14 15 15 15 15 15 16 16 16 16 16 17 17 17 17 17 18 18 18 18 18 19 19 19 19 19 20 20 21 21 21 21 21 22 22 22 22 22 23 23 23 23 23 24 24 24 24 24 25 25 25 25 25 26 26 26 26 26 27 27 27 27 27 28 28 28 28 28 29 29 29 29 29 30 30 30 30 30 31 31 31 31 31 32 32 32 32 32 33 33 33 33 33 34 34 34 34 34 35 35 35 35 35 36 36 36 36 36 37 37 37 37 37 38 38 38 38 38 39 39 39 39 39 40 40 40 40 40 41 41 41 41 41 42 42 42 42 42 43 43 43 43 43 44 44 44 44 44 45 45 45 45 45 46 46 46 46 46 47 47 47 47 47];
clear newmat

for i=1:229
    newmat(:,subjects(i),indices(i))=A(:,i);
end

clear count
% count(X,Y,Z)
% -- first dimension X: contains the #TRs in each of the 4 states; each row
% is a different state)
% -- second dimension Y: subject # (subjects 1-23 are stroke; 24-47 are
%    controls)
% -- third dimension Z: session # (all controls have 5 sessions. some
%    stroke subjects have less than 5 sessions)

% so to get the # of TRs in each of the 4 states for subject 1 in session
% 2, we would write:
% temp = count(:, 1, 2) (the : means include all rows)
% when we look at the temp variable we get:
% temp =
% 
%    116
%     20
%     15
%     22
% which means they spent 116 TR's in state 1, 20 in state 2, etc.

for b=1:47
    [count(:,b,1),~] = hist(newmat(:,b,1),1:numClusters);
    [count(:,b,2),~] = hist(newmat(:,b,2),1:numClusters);
    [count(:,b,3),~] = hist(newmat(:,b,3),1:numClusters);
    [count(:,b,4),~] = hist(newmat(:,b,4),1:numClusters);
    [count(:,b,5),~] = hist(newmat(:,b,5),1:numClusters);
end

%% Calculate Fractional Occupancy
TR=173;
stroke=count(:,1:23,1);
stroke_FO = stroke/TR

state1_stroke=(stroke_FO(1,:))
state2_stroke=(stroke_FO(2,:))
state3_stroke=(stroke_FO(3,:))
state4_stroke=(stroke_FO(4,:))
%state5_stroke=(stroke_FO(5,:))

controls=count(:,24:47,1)
controls_FO = controls/TR; 

state1_controls=(controls_FO(1,:))
state2_controls=(controls_FO(2,:))
state3_controls=(controls_FO(3,:))
state4_controls=(controls_FO(4,:))
%state5_controls=(controls_FO(5,:))

subplot(2,3,1)
group=[zeros(1,23),ones(1,24)]
boxplot([state1_stroke state1_controls],group)
title({'State 1 Group Difference', ' in Fractional Occupancy'})
xticklabels({'Stroke Patients','Controls'})
ylabel('Fractional Occupancy')
set(gca,'FontSize', 15)

subplot(2,3,2)
group=[zeros(1,23),ones(1,24)]
boxplot([state2_stroke state2_controls],group)
title({'State 2 Group Difference', ' in Fractional Occupancy'})
xticklabels({'Stroke Patients','Controls'})
ylabel('Fractional Occupancy')
set(gca,'FontSize', 15)

subplot(2,3,3)
group=[zeros(1,23),ones(1,24)]
boxplot([state3_stroke state3_controls],group)
title({'State 3 Group Difference', ' in Fractional Occupancy'})
xticklabels({'Stroke Patients','Controls'})
ylabel('Fractional Occupancy')
set(gca,'FontSize', 15)

subplot(2,3,4)
group=[zeros(1,23),ones(1,24)]
boxplot([state4_stroke state4_controls],group)
title({'State 4 Group Difference', ' in Fractional Occupancy'})
xticklabels({'Stroke Patients','Controls'})
ylabel('Fractional Occupancy')
set(gca,'FontSize', 15)

subplot(2,3,5)
group=[zeros(1,23),ones(1,24)]
boxplot([state5_stroke state5_controls],group)
title({'State 5 Group Difference', ' in Fractional Occupancy'})
xticklabels({'Stroke Patients','Controls'})
ylabel('Fractional Occupancy')
set(gca,'FontSize', 15)

[h,p,ci,stats]=ttest2(state1_stroke,state1_controls)
[h,p,ci,stats]=ttest2(state2_stroke,state2_controls)
[h,p,ci,stats]=ttest2(state3_stroke,state3_controls)
[h,p,ci,stats]=ttest2(state4_stroke,state4_controls)
[h,p,ci,stats]=ttest2(state5_stroke,state5_controls)

group=[zeros(1,23),ones(1,24)]
boxplot([state2_stroke state2_controls],group)
title('State 2 Group Difference in Fractional Occupancy')
xticklabels({'Stroke Patients','Controls'})
ylabel('Fractional Occupancy')
set(gca,'FontSize', 15)


[h,p,ci,stats]=ttest2(state1_stroke,state1_controls)
[h,p,ci,stats]=ttest2(state2_stroke,state2_controls)


TR=173;
stroke=countS2(:,1:23); 
stroke_FO = stroke/TR;

state4_stroke=(stroke_FO(4,:))

controls=countS2(:,24:47);
controls_FO = controls/TR; 
state4_controls=(controls_FO(4,:))
%state5_controls=(controls_FO(5,:))

group=[zeros(1,23),ones(1,24)]
boxplot([state4_stroke state4_controls],group)
title('State 4 Group Difference in Fractional Occupancy')
xticklabels({'Stroke Patients','Controls'})
ylabel('Fractional Occupancy')
set(gca,'FontSize', 15)

[h,p,ci,stats]=ttest2(state4_stroke,state4_controls)
[h,p,ci,stats]=ttest2(state5_stroke,state5_controls)


group=[zeros(1,23),ones(1,24)]
boxplot([state5_stroke state5_controls],group)
title('State 5 Group Difference in Fractional Occupancy')
xticklabels({'Stroke Patients','Controls'})
ylabel('Fractional Occupancy')
set(gca,'FontSize', 15)


[h,p,ci,stats]=ttest2(state4_stroke,state4_controls)
[h,p,ci,stats]=ttest2(state5_stroke,state5_controls)



%% Correlation between baseline impairment and FO time
fig=figure(1);
fig.Position=[0 0 1100 1200]

subplot(2,2,1)
plot(baselineFM', state1_stroke, '*r');
[rho,p]=corr(baselineFM, state1_stroke')
b=polyfit(baselineFM, state1_stroke,1);
a=polyval(b,baselineFM);
hold on;
plot(baselineFM, a, '-r')
title('Baseline Fugl Meyer vs FO in State 1')
set(gca, 'FontSize', 14)
xlabel('Fugl-Meyer score')
ylabel('Dwell Time (seconds)')
%ylim([2.5 5.5])
text(5, 5.25, ['correlation =' num2str(round(rho,3))], 'FontSize', 14)
text(5, 4.75, ['p=', num2str(round(p,3))], 'FontSize', 14)

subplot(2,2,2)
plot(baselineFM', state2_stroke, '*r');
[rho,p]=corr(baselineFM, state2_stroke');
b=polyfit(baselineFM, state2_stroke,1);
a=polyval(b,baselineFM);
hold on;
plot(baselineFM, a, '-r')
title('Baseline Fugl Meyer vs FO in State 2')
set(gca, 'FontSize', 14)
xlabel('Fugl-Meyer score')
ylabel('Dwell Time (seconds)')
%ylim([2.5 5.5])
text(5, 5.25, ['correlation =' num2str(round(rho,3))], 'FontSize', 14)
text(5, 4.75, ['p=', num2str(round(p,3))], 'FontSize', 14)

subplot(2,2,3)
plot(baselineFM, state3_stroke, '*b');
[rho,p]=corr(baselineFM, state3_stroke');
b=polyfit(baselineFM, state3_stroke,1);
a=polyval(b,baselineFM);
hold on;
plot(baselineFM, a, '-b')
title('Baseline Fugl Meyer vs Fractional Occupancy in State 3')
set(gca, 'FontSize', 14)
xlabel('Fugl-Meyer score')
ylabel('Fractional occupancy')
%ylim([0.42 0.58])
text(5, 0.565, ['correlation =' num2str(round(rho,3))], 'FontSize', 14)
text(5, 0.545, ['p=', num2str(round(p,3))], 'FontSize', 14)


subplot(2,2,4)
plot(baselineFM, state4_stroke, '*b');
[rho,p]=corr(baselineFM, state4_stroke')
b=polyfit(baselineFM, state4_stroke,1);
a=polyval(b,baselineFM);
hold on;
plot(baselineFM, a, '-b')
title('Baseline Fugl Meyer vs Fractional Occupancy in State 4')
set(gca, 'FontSize', 14)
xlabel('Fugl-Meyer score')
ylabel('Fractional occupancy')
%ylim([0.42 0.58])
text(5, 0.565, ['correlation =' num2str(round(rho,3))], 'FontSize', 14)
text(5, 0.545, ['p=', num2str(round(p,3))], 'FontSize', 14)



%% Calculate Dwell Time
dwell = []; %zeros(numClusters,220,89);
sumApp=0;
numClusters=4
for c=1:numClusters
    for b=1:47
        appear=find(newmat(:,b,1)==c);
        [maxIndex,~]=size(appear);
        appear=vertcat(appear,zeros(1,1)); %so that while statement will allow us to count the last index
        s=1;
        i=1;
        a=1;
        sumApp=sumApp+maxIndex;
        while a<maxIndex+1 
            if appear(a+1,1)==appear(a,1)+1 %if the next index value follows the current one in numerical order, then we are still dwelling in this state
                s=s+1;
                dwell(c,i,b)=s;
                a=a+1;
            else
                dwell(c,i,b)=s;
                i=i+1;
                a=a+1;
                s=1;
            end
        end
        %if sum(dwell(c,:,b)) ~= maxIndex
         %  disp(['Warning! Cluster ',num2str(c),' and Scan ',num2str(b),' sum does not match appearance count']);
        %end
    end
end

stroke_dwell=squeeze(dwell(:,:,1:23));
control_dwell=squeeze(dwell(:,:,24:47));

stroke_dwell_state1=stroke_dwell(1,:);
control_dwell_state1=control_dwell(1,:);

for i=1:23
    sub=stroke_dwell(1,:,i);
    dwell_avg_state1_stroke(i)=mean(sub(sub>0));
end
for i=1:24
    sub=control_dwell(1,:,i);
    dwell_avg_state1_control(i)=mean(sub(sub>0),2);
end
for i=1:23
    sub=stroke_dwell(2,:,i);
    dwell_avg_state2_stroke(i)=mean(sub(sub>0));
end
for i=1:24
    sub=control_dwell(2,:,i);
    dwell_avg_state2_control(i)=mean(sub(sub>0),2);
end
for i=1:23
    sub=stroke_dwell(3,:,i);
    dwell_avg_state3_stroke(i)=mean(sub(sub>0));
end
for i=1:24
    sub=control_dwell(3,:,i);
    dwell_avg_state3_control(i)=mean(sub(sub>0),2);
end
for i=1:23
    sub=stroke_dwell(4,:,i);
    dwell_avg_state4_stroke(i)=mean(sub(sub>0));
end
for i=1:24
    sub=control_dwell(4,:,i);
    dwell_avg_state4_control(i)=mean(sub(sub>0),2);
end
for i=1:23
    sub=stroke_dwell(5,:,i);
    dwell_avg_state5_stroke(i)=mean(sub(sub>0));
end
for i=1:24
    sub=control_dwell(5,:,i);
    dwell_avg_state5_control(i)=mean(sub(sub>0),2);
end

[h,p,ci,stats]=ttest2(dwell_avg_state1_stroke,dwell_avg_state1_control)
[h,p,ci,stats]=ttest2(dwell_avg_state2_stroke,dwell_avg_state2_control)

[h,p,ci,stats]=ttest2(dwell_avg_state3_stroke,dwell_avg_state3_control)
[h,p,ci,stats]=ttest2(dwell_avg_state4_stroke,dwell_avg_state4_control)
[h,p,ci,stats]=ttest2(dwell_avg_state5_stroke,dwell_avg_state5_control)


%% Correlation between baseline impairment and dwell time
fig=figure(1);
fig.Position=[0 0 1100 1200]

subplot(2,2,1)
plot(baselineFM', dwell_avg_state3_stroke, '*r');
[rho,p]=corr(baselineFM, dwell_avg_state3_stroke')
b=polyfit(baselineFM, dwell_avg_state3_stroke,1);
a=polyval(b,baselineFM);
hold on;
plot(baselineFM, a, '-r')
title('Baseline Fugl Meyer vs Dwell Time in State 1')
set(gca, 'FontSize', 14)
xlabel('Fugl-Meyer score')
ylabel('Dwell Time (seconds)')
%ylim([2.5 5.5])
text(5, 5.25, ['correlation =' num2str(round(rho,3))], 'FontSize', 14)
text(5, 4.75, ['p=', num2str(round(p,3))], 'FontSize', 14)

subplot(2,2,2)
plot(baselineFM', dwell_avg_state2_stroke, '*r');
[rho,p]=corr(baselineFM, dwell_avg_state2_stroke');
b=polyfit(baselineFM, dwell_avg_state2_stroke,1);
a=polyval(b,baselineFM);
hold on;
plot(baselineFM, a, '-r')
title('Baseline Fugl Meyer vs Dwell Time in State 2')
set(gca, 'FontSize', 14)
xlabel('Fugl-Meyer score')
ylabel('Dwell Time (seconds)')
%ylim([2.5 5.5])
text(5, 5.25, ['correlation =' num2str(round(rho,3))], 'FontSize', 14)
text(5, 4.75, ['p=', num2str(round(p,3))], 'FontSize', 14)

subplot(2,2,3)
plot(baselineFM, state1_stroke, '*b');
[rho,p]=corr(baselineFM, state1_stroke');
b=polyfit(baselineFM, state1_stroke,1);
a=polyval(b,baselineFM);
hold on;
plot(baselineFM, a, '-b')
title('Baseline Fugl Meyer vs Fractional Occupancy in State 1')
set(gca, 'FontSize', 14)
xlabel('Fugl-Meyer score')
ylabel('Fractional occupancy')
%ylim([0.42 0.58])
text(5, 0.565, ['correlation =' num2str(round(rho,3))], 'FontSize', 14)
text(5, 0.545, ['p=', num2str(round(p,3))], 'FontSize', 14)


subplot(2,2,4)
plot(baselineFM, state2_stroke, '*b');
[rho,p]=corr(baselineFM, state2_stroke')
b=polyfit(baselineFM, state2_stroke,1);
a=polyval(b,baselineFM);
hold on;
plot(baselineFM, a, '-b')
title('Baseline Fugl Meyer vs Fractional Occupancy in State 2')
set(gca, 'FontSize', 14)
xlabel('Fugl-Meyer score')
ylabel('Fractional occupancy')
%ylim([0.42 0.58])
text(5, 0.565, ['correlation =' num2str(round(rho,3))], 'FontSize', 14)
text(5, 0.545, ['p=', num2str(round(p,3))], 'FontSize', 14)


%% Correlation between realized recovery and dwell times.
potentialrecovery=100-basline;
truerecovery=change;

realizedrecovery=truerecovery./potentialrecovery;
realizedrecovery(10)=0;
realizedrecovery(14)=0;
realizedrecovery(15)=0;
save('/Users/emilyolafson/Documents/Thesis/SUB1_23_data/realizedrecovery.mat', 'realizedrecovery')

fig=figure(2);
fig.Position=[0 0 1100 1200]

subplot(2,2,1)
plot(realizedrecovery', dwell_avg_state1_stroke, '*r');
[rho,p]=corr(realizedrecovery, dwell_avg_state1_stroke');
b=polyfit(realizedrecovery, dwell_avg_state1_stroke,1);
a=polyval(b,realizedrecovery);
hold on;
plot(realizedrecovery, a, '-r')
title('Realized recovery vs Dwell Time in State 1')
set(gca, 'FontSize', 14)
xlabel('Realized recovery')
ylabel('Dwell Time (seconds)')
ylim([2.5 5.5])
text(5, 5.25, ['correlation =' num2str(round(rho,3))], 'FontSize', 14)
text(5, 4.75, ['p=', num2str(round(p,3))], 'FontSize', 14)

subplot(2,2,2)
plot(realizedrecovery', dwell_avg_state2_stroke, '*r');
[rho,p]=corr(realizedrecovery, dwell_avg_state2_stroke');
b=polyfit(realizedrecovery, dwell_avg_state2_stroke,1);
a=polyval(b,realizedrecovery);
hold on;
plot(realizedrecovery, a, '-r')
title('Realized recovery vs Dwell Time in State 2')
set(gca, 'FontSize', 14)
xlabel('Realized recovery')
ylabel('Dwell Time (seconds)')
ylim([2.5 5.5])
text(5, 5.25, ['correlation =' num2str(round(rho,3))], 'FontSize', 14)
text(5, 4.75, ['p=', num2str(round(p,3))], 'FontSize', 14)

subplot(2,2,3)
plot(realizedrecovery, state1_stroke, '*b');
[rho,p]=corr(realizedrecovery, state1_stroke');
b=polyfit(realizedrecovery, state1_stroke,1);
a=polyval(b,realizedrecovery);
hold on;
plot(realizedrecovery, a, '-b')
title('Realized recovery vs Fractional Occupancy in State 1')
set(gca, 'FontSize', 14)
xlabel('Realized recovery')
ylabel('Fractional occupancy')
ylim([0.42 0.58])
text(5, 0.565, ['correlation =' num2str(round(rho,3))], 'FontSize', 14)
text(5, 0.545, ['p=', num2str(round(p,3))], 'FontSize', 14)


subplot(2,2,4)
plot(realizedrecovery, state2_stroke, '*b');
[rho,p]=corr(realizedrecovery, state2_stroke')
b=polyfit(realizedrecovery, state2_stroke,1);
a=polyval(b,realizedrecovery);
hold on;
plot(realizedrecovery, a, '-b')
title('Realized recovery vs Fractional Occupancy in State 2')
set(gca, 'FontSize', 14)
xlabel('Realized recovery')
ylabel('Fractional occupancy')
ylim([0.42 0.58])
text(5, 0.565, ['correlation =' num2str(round(rho,3))], 'FontSize', 14)
text(5, 0.545, ['p=', num2str(round(p,3))], 'FontSize', 14)
