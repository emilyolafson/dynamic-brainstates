%% plot results of permutation pvalue testing of transition probabilities
% generate data with calculate_transition_probabilities.m

figdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/figures/'
resdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/'


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

% visualize results
tiledlayout(1,5)
nexttile;imagesc(p1); caxis([0 0.05]);nexttile;imagesc(p2);caxis([0 0.05]);nexttile;imagesc(p3);caxis([0 0.05]);nexttile;imagesc(p4);caxis([0 0.05]);nexttile;imagesc(p5);caxis([0 0.05]);
colormap 'hot'
