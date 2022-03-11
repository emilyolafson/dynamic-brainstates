%% plot results of permutation pvalue testing of transition probabilities
% generate data with calculate_transition_probabilities.m

figdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/shen268/k4/'
close all;
figure('Position', [0 0 1500 300])

tiledlayout(1,5,'padding', 'none')
% session 1
nexttile;imagesc(stats1); colorbar;title('Session 1'); 
signif=(allpvals{1}<0.05)-(allpvals_adj{1}<0.05)
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.05, '*', 'FontSize', 30, 'Color', 'white')
end
signif=allpvals{1}<0.05;
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.05, '*', 'FontSize', 30, 'Color', 'white')
end
caxis([-3, 3])
yticks(1:4)
ylabel('Current state')
xlabel('Next state')
set(gca, 'FontSize', 15)

% session 2
nexttile;imagesc(stats2); colorbar;title('Session 2'); 
signif=allpvals_adj{2}<0.05;
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.05, '*', 'FontSize', 30, 'Color', 'white')
end
signif=(allpvals{2}<0.05)-(allpvals_adj{2}<0.05)
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.05, '*', 'FontSize', 30, 'Color', 'white')
end
caxis([-3, 3])
yticks(1:4)
ylabel('Current state')
xlabel('Next state')
set(gca, 'FontSize', 15)

% session 3
nexttile;imagesc(stats3); colorbar;title('Session 3'); 
signif=allpvals_adj{3}<0.05;
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.05, '*', 'FontSize', 30, 'Color', 'white')
end
signif=(allpvals{3}<0.05)-(allpvals_adj{3}<0.05)
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.05, '*', 'FontSize', 30, 'Color', 'white')
end
caxis([-3, 3])
yticks(1:4)
ylabel('Current state')
xlabel('Next state')
set(gca, 'FontSize', 15)

% session 4
nexttile;imagesc(stats4); colorbar;title('Session 4'); 
signif=allpvals_adj{4}<0.05;
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.05, '*', 'FontSize', 30, 'Color', 'white')
end
signif=(allpvals{4}<0.05)-(allpvals_adj{4}<0.05)
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.05, '*', 'FontSize', 30, 'Color', 'white')
end
caxis([-3, 3])
yticks(1:4)
ylabel('Current state')
xlabel('Next state')
set(gca, 'FontSize', 15)

% session 5
nexttile;imagesc(stats5); colorbar;title('Session 5'); 
signif=allpvals_adj{5}<0.05
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.1,rows(i)+0.05, '**', 'FontSize', 30, 'Color', 'white')
end
signif=(allpvals{5}<0.05)-(allpvals_adj{5}<0.05)
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.05, '*', 'FontSize', 30, 'Color', 'white')
end
caxis([-3, 3])
yticks(1:4)
ylabel('Current state')
xlabel('Next state')
colormap(flipud(brewermap([],'PuOr')))
set(gca, 'FontSize', 15)


saveas(gcf,strcat(figdir, 'transitionprobs_strokecontrol_withrealpersist.png'))



% plot transition probabilities - averaged across all
% with persist
figure('Position', [0 0 1000 500])
tiledlayout(1,2)
nexttile;imagesc(stats_all); colorbar; 
colormap(flipud(brewermap([],'PuOr')))

signif=p_adj_all<0.05
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.1,rows(i)+0.05, '**', 'FontSize', 30, 'Color', 'white')
end
signif=(p_all<0.05)-(p_adj_all<0.05)
[rows,cols]=find(signif)
for i=1:size(rows,1)
    text(cols(i)-0.01,rows(i)+0.05, '*', 'FontSize', 30, 'Color', 'white')
end
caxis([-3, 3])
yticks(1:4)
ylabel('Current state')
xlabel('Next state')
set(gca, 'FontSize', 20)
xticks(1:4)

saveas(gcf,strcat(figdir, 'transitionprobs_strokecontrol_withrealpersist_avg.png'))























