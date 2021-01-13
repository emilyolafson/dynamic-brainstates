clear all; close all;clc
basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);
addpath(genpath('code'))
%% set inputs
numClusters = %define;
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
load(fullfile(savedir,['Partition',num2str(numClusters),'.mat']))
TR = 220;


%% count number of each cluster per scan
A=reshape(partition,220,[]);
count = zeros(numClusters,89);
for b=1:89
    [count(:,b),~] = hist(A(:,b),1:numClusters);
%     for a=1:220
%         if A(a,b)==1
%             count(1,b)=count(1,b)+1;
%         end
%         if A(a,b)==2
%             count(2,b)=count(2,b)+1;
%         end
%         if A(a,b)==3
%             count(3,b)=count(3,b)+1;
%         end
%         if A(a,b)==4
%             count(4,b)=count(4,b)+1;
%         end
%         if A(a,b)==5
%             count(5,b)=count(5,b)+1;
%         end
%         if A(a,b)==6
%             count(6,b)=count(6,b)+1;
%         end
%         if A(a,b)==7
%             count(7,b)=count(7,b)+1;
%         end
%         if A(a,b)==8
%             count(8,b)=count(8,b)+1;
%         end
%     end
%     for i=1:numClusters
%         if count(i,:)==0
%             disp(['Warning! Scan ',num2str(b),' is missing cluster ',num2str(i)]);
%         end
%     end
%     
%     count(7,b) = sum(count(1:6,b));
%     tot = count(7,b);
%     if tot ~= 220
%         disp(['Warning! Scan ',num2str(b),' total does not equal 200! ',num2str(i)]);
%     end
end

%% Calculate Fractional Occupancy
LSD=count(:,[1 3:45]); %skip subj 2 condition 1 (did not take placebo scan)
LSDfo = LSD/TR;

PL=count(:,46:89);
PLfo = PL/TR;

%% Calculate Dwell Time and Appearance Rate
dwell = []; %zeros(numClusters,220,89);
sumApp=0;
for c=1:numClusters
    for b=1:89
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


LSD = dwell(:,:,[1,3:45]); %skip subj 2, condition one (no corresponding placebo scan)
LSD_count=zeros(6,3,44); %rows represent each cluster,column 1=total time, 2=total number appearances,3=avg dwell time
LSDdt=[]; 
LSDar=[];
[~,maxIndex,~]=size(dwell);
for c=1:numClusters
    for b=1:44
        LSD_count(c,1,b)=LSD_count(c,1,b)+sum(LSD(c,:,b));
        a=1;
        while a<=maxIndex
            if LSD(c,a,b)~= 0
                LSD_count(c,2,b)=LSD_count(c,2,b)+1;
                a=a+1;
            else
                a=a+1;
            end
        end
        LSDdt(c,b)=LSD_count(c,1,b)/LSD_count(c,2,b);
        LSDar(c,b)=LSD_count(c,2,b)/7.33; %appearance rate per minute = tot. appear / 7 min 20 s scan
    end
end


        
        
PL=dwell(:,:,46:89);
PL_count=zeros(6,3,44); %rows represent each cluster,column 1=total time, 2=total number appearances,3=avg dwell time
PLdt=[];
PLar=[];
[~,maxIndex,~]=size(dwell);
for c=1:numClusters
    for b=1:44
        PL_count(c,1,b)=PL_count(c,1,b)+sum(PL(c,:,b));
        a=1;
        while a<=maxIndex
            if PL(c,a,b)~= 0
                PL_count(c,2,b)=PL_count(c,2,b)+1;
                a=a+1;
            else
                a=a+1;
            end
        end
        PLdt(c,b)=PL_count(c,1,b)/PL_count(c,2,b);
        PLar(c,b)=PL_count(c,2,b)/7.33; %appearance rate per minute = tot. appear / 7 min 20 s scan
    end
end

%% F-test for equality of variances

%combine pre+post conditions into one group and during music into another

dvar = zeros(6,3,2); %get a quick measure of %differences in variance for the two groups

for i=1:numClusters
    dvar(i,1,1) = (var(LSDfo(i,[1:14 30:44]))-var(PLfo(i,[1:14 30:44])))/(var(LSDfo(i,[1:14 30:44])))*100;
    dvar(i,2,1) = (var(LSDdt(i,[1:14 30:44]))-var(PLdt(i,[1:14 30:44])))/(var(LSDdt(i,[1:14 30:44])))*100;
    dvar(i,3,1) = (var(LSDar(i,[1:14 30:44]))-var(PLar(i,[1:14 30:44])))/(var(LSDar(i,[1:14 30:44])))*100;
    
    dvar(i,1,2) = (var(LSDfo(i,15:29))-var(PLfo(i,15:29)))/(var(LSDfo(i,15:29)))*100;
    dvar(i,2,2) = (var(LSDdt(i,15:29))-var(PLdt(i,15:29)))/(var(LSDdt(i,15:29)))*100;
    dvar(i,3,2) = (var(LSDar(i,15:29))-var(PLar(i,15:29)))/(var(LSDar(i,15:29)))*100;
end

pvaluesPP = zeros(6,3);
pvaluesdur = zeros(6,3);

%

for i=1:numClusters
    [~,pvaluesPP(i,1)] = vartest2(LSDfo(i,[1:14 30:44]),PLfo(i,[1:14 30:44]));
    [~,pvaluesPP(i,2)] = vartest2(LSDdt(i,[1:14 30:44]),PLdt(i,[1:14 30:44]));
    [~,pvaluesPP(i,3)] = vartest2(LSDar(i,[1:14 30:44]),PLar(i,[1:14 30:44]));
    
    [~,pvaluesdur(i,1)] = vartest2(LSDfo(i,15:29),PLfo(i,15:29));
    [~,pvaluesdur(i,2)] = vartest2(LSDdt(i,15:29),PLdt(i,15:29));
    [~,pvaluesdur(i,3)] = vartest2(LSDar(i,15:29),PLar(i,15:29));
end

pvaluesPP = reshape(mafdr(reshape(pvaluesPP,1,18),"BHFDR",1),6,3);
pvaluesdur = reshape(mafdr(reshape(pvaluesdur,1,18),"BHFDR",1),6,3);

%% Detect outliers

outPP = NaN(6,29,6);
outdur = NaN(6,15,6);

for i=1:numClusters
    outPP(i,:,1) = isoutlier(LSDfo(i,[1:14 30:44]));
    outPP(i,:,2) = isoutlier(LSDdt(i,[1:14 30:44]));
    outPP(i,:,3) = isoutlier(LSDar(i,[1:14 30:44]));
    outPP(i,:,4) = isoutlier(PLfo(i,[1:14 30:44]));
    outPP(i,:,5) = isoutlier(PLdt(i,[1:14 30:44]));
    outPP(i,:,6) = isoutlier(PLar(i,[1:14 30:44]));
    
    outdur(i,:,1) = isoutlier(LSDfo(i,15:29));
    outdur(i,:,2) = isoutlier(LSDdt(i,15:29));
    outdur(i,:,3) = isoutlier(LSDar(i,15:29));
    outdur(i,:,4) = isoutlier(PLfo(i,15:29));
    outdur(i,:,5) = isoutlier(PLdt(i,15:29));
    outdur(i,:,6) = isoutlier(PLar(i,15:29));
    
end

sumPP(1,:) = sum(sum(outPP(:,:,1:3),3),1);
sumPP(2,:) = sum(sum(outPP(:,:,4:6),3),1);
sumdur(1,:) = sum(sum(outdur(:,:,1:3),3),1);
sumdur(2,:) = sum(sum(outdur(:,:,4:6),3),1);
    

%% Save 

clusters=char(clusterNames);
save(fullfile(savedir,['ViolinData',num2str(numClusters),'.mat']), 'LSDfo', 'PLfo', 'LSDdt', 'PLdt', 'LSDar', 'PLar', 'clusters')


