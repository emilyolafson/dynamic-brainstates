%% for visualization 

atlasblobs_list=load('./atlasblobs_saved.mat');
atlasblobs_list=atlasblobs_list.atlasblobs_list;

whichatlas={'shen268'}

%for testing just assign y-axis position (AP axis) as value for each ROI
roivals=atlasblobs_list(strcmpi({atlasblobs_list.atlasname},whichatlas)).roicenters(:,2);
roivals=rand(size(roivals));
%cmap=colormap_videen(256);
%cmap=cmap(129:256,:);
%cmap=hot(256);
%cmap=cmap(1:254,:);

%cmap=cmap(129:256,:);

%for volume
cmap=parula(268);
%for cmap_mask=613:1007
    %cmap(cmap_mask,:)=[0.5, 0.5, 0.5];
%end

%for HE
%for cmap_mask=167:753
%    cmap(cmap_mask,:)=[0.5, 0.5, 0.5];
%end
%'clim', [-8.1, 8.1],...

%whichatlas_all={'cc400','cc200','aal','ez','ho','tt','fs86'};
close all;
clc;

%cc400_data needs to be a 1x392 vector
%shen368 needs to be 1x268
data = x;
data_min=min(data);
data_max=max(data);
img=display_atlas_blobs(data,atlasblobs_list,...
    'atlasname',whichatlas,...
    'render',true,...
    'backgroundimage',true,...
    'crop',true,...
    'colormap',cmap,...
    'alpha', data);
%rescale(mean_fi(2,:).^2)
%'alpha',rescale(pinv(rescale(abs(cc400_fi')))).^2
%'alpha',cc400_t_facealpha
%'alpha', rescale(abs(mean_fs86_size_featimp))
figure;
imshow(img);
c=colorbar('SouthOutside', 'fontsize', 20);
c.Label.String='Colorbar Label';
set(gca,'colormap',cmap);
caxis([data_min data_max]);
%title('TT', 'fontsize',18);
annotation('textbox',[.835 .38 .1 .2],'String','RH','EdgeColor','none', 'fontsize',20,'color','white')
annotation('textbox',[.12 .38 .1 .2],'String','LH','EdgeColor','none', 'fontsize',20,'color','white')
annotation('textbox',[.45 .9 .1 .04],'String','Lateral','EdgeColor','none', 'fontsize',20,'color','white', 'horizontalalignment','center','backgroundcolor','black')
annotation('textbox',[.45 .21 .1 .04],'String','Medial','EdgeColor','none', 'fontsize',20,'color','white', 'horizontalalignment','center', 'backgroundcolor','black')


%imwrite(img,sprintf('~/Downloads/%s_t.png',whichatlas));