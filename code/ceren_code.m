%% Visualize cluster centroids
% compute centroids and plot
best_number_of_clusters =5 %define

centroids = GET_CENTROIDS(np,cluster_output(:,5),best_number_of_clusters);
atlasblobs_list=load('gummibrain/atlasblobs_saved.mat');
atlasblobs_list=atlasblobs_list.atlasblobs_list;

whichatlas={'shen268'};
cmap=colormap(plasma) %matplotlib colormap matlab toolbox

% gummibrain code (need to install from kjamison repo) https://github.com/kjamison/atlasblobs
results_dir='/Users/emilyolafson/GIT/dynamic-brainstates/results/centroids_5states/'
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
f=figure('Position', [0 0 900 800]);
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
    set(ax,'FontSize',12);
    %title(overallNames{K},'Color',clusterColors(K,:),'FontSize',8);
end


%% Count number of each cluster per scan
% Calculate fractional occupancy & 
%partition=cluster_output2(:,4)
clear A
partition=cluster_output(:,best_number_of_clusters)

clear stroke_count
clear control_count
clear *count_state*
clear appear
clear *FO*

z=1;
l=1;
k=best_number_of_clusters

    %leng = $ of TR's
    %n = number of subjects
for i=1:n
    for j=1
        clusters{i,j}=partition(z:z+leng-1);
        counts=hist(clusters{i,j}, 1:5);
        stroke_count_state1{i,j}=counts(1);
        stroke_count_state2{i,j}=counts(2);
        stroke_count_state3{i,j}=counts(3);
        stroke_count_state4{i,j}=counts(4);
        stroke_count_state5{i,j}=counts(5);

        stroke_FO1{i,j}=stroke_count_state1{i,j}/leng;
        stroke_FO2{i,j}=stroke_count_state2{i,j}/leng;
        stroke_FO3{i,j}=stroke_count_state3{i,j}/leng;
        stroke_FO4{i,j}=stroke_count_state4{i,j}/leng;
        stroke_FO5{i,j}=stroke_count_state5{i,j}/leng;

        dwell=zeros(size(clusters{i,j},1),1)
        c=1;
        for k=1:best_number_of_clusters
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
        z=z+leng;
        l=l+1;
    end
end

%% Calculate transition probabilities
%GET_TRANS_PROBS(partition,subjInd,numClusters)

% need to provide subjInd, which is a long vector with length [leng * n]
% subjInd: integer vector, subject index for partition

subjID=[];
for i=1:n
    subjID = [subjID; ones(leng, 1)*i];
end

[transProbs,transitionProbabilityMats,numTransitions] = GET_TRANS_PROBS(partition,subjID, 5);  

for i=1:n
    vec=i
    tmp=transitionProbabilityMats(vec,:,:);
    for j=1 
        for p=1:best_number_of_clusters
            tran_prob{i,j}=squeeze(tmp(j,:,:));
        end
    end
end


%% get appearance rates
clear *appearance*
for i=1:n
    vector=clusters{i,j};
    appear1=0;
    appear2=0;
    appear3=0;
    appear4=0;
    appear5=0;

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
        elseif vector(index)==5
            disp('New appearance: state 5')
            appear5=appear5+1;
            try
                while vector(index+c)==5
                    disp('Still in state 5')
                    c=c+1;
                end
            end
            disp(['Leaving state 5'])
            index=index+c;
            continue;
        end
    end
    stroke_appearance1{i,j}=appear1/leng*3/60);
    stroke_appearance2{i,j}=appear2/leng*3/60);
    stroke_appearance3{i,j}=appear3/leng*3/60);
    stroke_appearance4{i,j}=appear4/leng*3/60);
    stroke_appearance5{i,j}=appear5/leng*3/60);
end


% Set data to 0 when subjects missing fmri

%stroke
stroke_count_state1=cell2mat(stroke_count_state1)
stroke_count_state2=cell2mat(stroke_count_state2)
stroke_count_state3=cell2mat(stroke_count_state3)
stroke_count_state4=cell2mat(stroke_count_state4)
stroke_count_state5=cell2mat(stroke_count_state5)

stroke_FO1=cell2mat(stroke_FO1)
stroke_FO2=cell2mat(stroke_FO2)
stroke_FO3=cell2mat(stroke_FO3)
stroke_FO4=cell2mat(stroke_FO4)
stroke_FO5=cell2mat(stroke_FO5)

stroke_appearance1=cell2mat(stroke_appearance1)
stroke_appearance2=cell2mat(stroke_appearance2)
stroke_appearance3=cell2mat(stroke_appearance3)
stroke_appearance4=cell2mat(stroke_appearance4)
stroke_appearance5=cell2mat(stroke_appearance5)

%controls
control_count_state1=cell2mat(control_count_state1)
control_count_state2=cell2mat(control_count_state2)
control_count_state3=cell2mat(control_count_state3)
control_count_state4=cell2mat(control_count_state4)
control_count_state5=cell2mat(control_count_state5)

control_FO1=cell2mat(control_FO1)
control_FO2=cell2mat(control_FO2)
control_FO3=cell2mat(control_FO3)
control_FO4=cell2mat(control_FO4)
control_FO5=cell2mat(control_FO5)

control_appearance1=cell2mat(control_appearance1)
control_appearance2=cell2mat(control_appearance2)
control_appearance3=cell2mat(control_appearance3)
control_appearance4=cell2mat(control_appearance4)
control_appearance5=cell2mat(control_appearance5)

