
baselineFM=load('baselineFM.mat');
baselineFM=baselineFM.basline;
save('baseline_scores.mat', 'baselineFM')

%load in functional mri scans from all 5 sessions
allsubs=load('ts_GSR_shen268_allsub_allsess.mat'); 
allsubs=allsubs.allsubs;
allsubs_fmri=reshape(allsubs, [1 235]);


%% Determine the optimal number of clusters using variance explained
allsubs_fmri=allsubs_fmri(find(~cellfun(@isempty,allsubs_fmri)))
allsubs_fmri=allsubs_fmri';

numscans1_11=[4 3 3 3 3 3 3 3 3 3 3;4 3 3 3 3 3 3 4 3 3 3;3 3 3 3 3 3 3 3 3 3 3;4 4 3 3 3 3 3 3 3 3 3 ; 3 3 3 3 3 0 3 3 3 3 3];
numscans12_23=[3 3 3 3 3 2 2 2 2 2 2 2; 3 3 3 4 3 2 2 2 2 2 2 2 ;3 3 3 3 3 2 2 2 0 2 2 2;0 3 3 3 3 2 2 2 0 2 2 2;0 3 4 3 3 2 2 2 0 2 2 2];
nscans = [numscans1_11, numscans12_23];
nscan = ones(5, 24)+4;
nscans = [nscans, nscan];

np= [];

for i=1:229
    a=allsubs_fmri{i};
    b=cell2mat(a);
    leng(i)=size(b,1);
    b=b(4:173,:);
    np=[np;b];
end
size(b)

distanceMethod= 'sqeuclidean'
nreps = 50;
totSum=[];

% should take about 5 mins to run for each cluster #
clear cluster_output
clear sumd
for i=3:5
    [cluster_output(:,i),~,sumd]=kmeans(np,i,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',1000);
    totSum(i)=sum(sumd);
end

best_number_of_clusters = 3 %define

%% Visualize cluster centroids
% compute centroids and plot
centroids = GET_CENTROIDS(np,cluster_output(:,best_number_of_clusters),best_number_of_clusters);

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
partition=parts;
clear A
A(:,1)=partition(1:47)
A(:,2)=partition(48:94)
A(:,3)=[partition(95:114);0;partition(115:140)]
A(:,4)=[partition(141:152);0;partition(153:159);0;partition(160:185)]
A(:,5)=[partition(186:191);0;partition(192:196);0;partition(197:203);0;partition(204:229)]


count = zeros(4,47);
for b=1:170
    [count(:,b),~] = hist(A(:,b),1:4);
end
%% Calculate Fractional Occupancy
TR=146;
stroke=count(:,1:23); 
stroke_FO = stroke/TR;

%% Calculate Fractional Occupancy
TR=146;
stroke=count(:,1:23); 
stroke_FO = stroke/TR;

state1_stroke=(stroke_FO(1,:))
state2_stroke=(stroke_FO(2,:))

controls=count(:,24:47);
controls_FO = controls/TR; 
state1_controls=(controls_FO(1,:))
state2_controls=(controls_FO(2,:))
group=[zeros(1,23),ones(1,24)]
boxplot([state2_stroke state2_controls],group)
title('State 2 Group Difference in Fractional Occupancy')
xticklabels({'Stroke Patients','Controls'})
ylabel('Fractional Occupancy')
set(gca,'FontSize', 15)

[h,p,ci,stats]=ttest2(state1_stroke,state1_controls)
[h,p,ci,stats]=ttest2(state2_stroke,state2_controls)

%% Calculate Dwell Time
dwell = []; %zeros(numClusters,220,89);
sumApp=0;

for c=1:numClusters
    for b=1:47
        appear=find(A(:,b)==c);
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
        if sum(dwell(c,:,b)) ~= maxIndex
           disp(['Warning! Cluster ',num2str(c),' and Scan ',num2str(b),' sum does not match appearance count']);
        end
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

[h,p,ci,stats]=ttest2(dwell_avg_state1_stroke,dwell_avg_state1_control)
[h,p,ci,stats]=ttest2(dwell_avg_state2_stroke,dwell_avg_state2_control)

histogram(dwell_avg_state1_stroke)
hold on;
histogram(dwell_avg_state1_control)

boxplot(dwell_avg_state1_stroke)

%% Correlation between baseline impairment and dwell time
fig=figure(1);
fig.Position=[0 0 1100 1200]

subplot(2,2,1)
plot(baselineFM', dwell_avg_state1_stroke, '*r');
[rho,p]=corr(baselineFM, dwell_avg_state1_stroke');
b=polyfit(baselineFM, dwell_avg_state1_stroke,1);
a=polyval(b,baselineFM);
hold on;
plot(baselineFM, a, '-r')
title('Baseline Fugl Meyer vs Dwell Time in State 1')
set(gca, 'FontSize', 14)
xlabel('Fugl-Meyer score')
ylabel('Dwell Time (seconds)')
ylim([2.5 5.5])
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
ylim([2.5 5.5])
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
ylim([0.42 0.58])
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
ylim([0.42 0.58])
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
