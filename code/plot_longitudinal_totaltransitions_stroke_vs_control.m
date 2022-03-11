
figdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/figures/'
resdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/'
stroke_meanappearance1=mean([stroke_appearance1(:,1),stroke_appearance2(:,1),stroke_appearance3(:,1),stroke_appearance4(:,1)],2)
stroke_meanappearance2=mean([stroke_appearance1(:,2),stroke_appearance2(:,2),stroke_appearance3(:,2),stroke_appearance4(:,2)],2)
stroke_meanappearance3=mean([stroke_appearance1(:,3),stroke_appearance2(:,3),stroke_appearance3(:,3),stroke_appearance4(:,3)],2)
stroke_meanappearance4=mean([stroke_appearance1(:,4),stroke_appearance2(:,4),stroke_appearance3(:,4),stroke_appearance4(:,4)],2)

control_meanappearance1=mean([control_appearance1(:,1),control_appearance2(:,1),control_appearance3(:,1),control_appearance4(:,1)],2)
control_meanappearance2=mean([control_appearance1(:,2),control_appearance2(:,2),control_appearance3(:,2),control_appearance4(:,2)],2)
control_meanappearance3=mean([control_appearance1(:,3),control_appearance2(:,3),control_appearance3(:,3),control_appearance4(:,3)],2)
control_meanappearance4=mean([control_appearance1(:,4),control_appearance2(:,4),control_appearance3(:,4),control_appearance4(:,4)],2)

close all;

figure('Position', [0 0 900 300])
idx=[zeros(23, 1); ones(24,1); ...
    ones(23, 1)*2; ones(24,1)*3;...
    ones(23, 1)*4; ones(24,1)*5;...
    ones(23, 1)*6; ones(24,1)*7
    ]

violin=violinplot([stroke_meanappearance1;control_meanappearance1;stroke_meanappearance2;control_meanappearance2;stroke_meanappearance3;control_meanappearance3;stroke_meanappearance4;control_meanappearance4], idx)

violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]
violin(3).ViolinColor=[0.4 0.4 0.8]
violin(4).ViolinColor=[1 0.4 0.4]
violin(5).ViolinColor=[0.4 0.4 0.8]
violin(6).ViolinColor=[1 0.4 0.4]
violin(7).ViolinColor=[0.4 0.4 0.8]
violin(8).ViolinColor=[1 0.4 0.4]

xticks(1:2:12)
xticklabels({'Session 1','','Session 2','','Session 3', '','Session 4'})
set(gca, 'FontSize', 15)
ylabel('Average # of transitions (to all states)')

[h, p(1), ci, stats]=ttest2(stroke_meanappearance1,control_meanappearance1)
[h, p(2), ci, stats]=ttest2(stroke_meanappearance2,control_meanappearance3)
[h, p(3), ci, stats]=ttest2(stroke_meanappearance3,control_meanappearance4)
[h, p(4), ci, stats]=ttest2(stroke_meanappearance4,control_meanappearance4)

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

saveas(gcf,strcat(figdir, 'Avg_number_transitions_overtime_stroke_vs_control.png'))

