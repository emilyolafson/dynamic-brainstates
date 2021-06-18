
figdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/5states/figures/'
resdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/5states/'

%% Analysis - Fractional occupancy and dwell time between stroke & control
% FO differences across all sessions, stroke and controls;
% mean over time
close all;
clear p
figure('Position', [0 0 900 300])
idx=[zeros(23, 1); ones(24,1); ...
    ones(23, 1)*2; ones(24,1)*3;...
    ones(23, 1)*4; ones(24,1)*5;...
    ones(23, 1)*6; ones(24,1)*7
    ]

violin=violinplot([mean(stroke_FO1, 2, 'omitnan')',mean(control_FO1,  2, 'omitnan')',...
    mean(stroke_FO2, 2, 'omitnan')',mean(control_FO2, 2, 'omitnan')',...
    mean(stroke_FO3,  2, 'omitnan')',mean(control_FO3, 2, 'omitnan')',...
    mean(stroke_FO4,  2, 'omitnan')',mean(control_FO4,  2, 'omitnan')'
    ], idx)
ylabel('Fractional Occupancy')
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]
violin(3).ViolinColor=[0.4 0.4 0.8]
violin(4).ViolinColor=[1 0.4 0.4]
violin(5).ViolinColor=[0.4 0.4 0.8]
violin(6).ViolinColor=[1 0.4 0.4]
violin(7).ViolinColor=[0.4 0.4 0.8]
violin(8).ViolinColor=[1 0.4 0.4]
xticks(1:8)

yt = get(gca, 'YTick');
xt = get(gca, 'XTick');

[h, p(1), ci, stats]=ttest2(mean(stroke_FO1, 2, 'omitnan'),mean(control_FO1, 2))
[h, p(2), ci, stats]=ttest2(mean(stroke_FO2,  2, 'omitnan'),mean(control_FO2, 2))
[h, p(3), ci, stats]=ttest2(mean(stroke_FO3,  2, 'omitnan'),mean(control_FO3, 2))
[h, p(4), ci, stats]=ttest2(mean(stroke_FO4,  2, 'omitnan'),mean(control_FO4, 2))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end

old=0

for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.1, max(yt)*1.15, '*k');
    end
end
xticklabels({'State 1','    ','State 2','     ','State 3', '','State 4'})
set(gca, 'FontSize', 15)



save(strcat(resdir, 'p_adjusted_Timeavg_FO_stroke_vs_control.mat'), 'p_adj'))
save(strcat(resdir, 'p_unc_Timeavg_FO_stroke_vs_control.mat'), 'p')
saveas(gcf,strcat(figdir, 'FO_Timeavg_stroke_vs_control.png'))

%% Dwell time differences across all sessions, stroke and controls;
close all;
% mean over time
figure('Position', [0 0 900 300])
idx=[zeros(23, 1); ones(24,1); ...
    ones(23, 1)*2; ones(24,1)*3;...
    ones(23, 1)*4; ones(24,1)*5;...
    ones(23, 1)*6; ones(24,1)*7
    ]

violin=violinplot([mean(dwell_avg_stroke(:,:,1), 2,'omitnan')',mean(dwell_avg_control(:,:,1), 2)',...
    mean(dwell_avg_stroke(:,:,2), 2,'omitnan')',mean(dwell_avg_control(:,:,2), 2)',...
    mean(dwell_avg_stroke(:,:,3), 2,'omitnan')',mean(dwell_avg_control(:,:,3), 2)',...
    mean(dwell_avg_stroke(:,:,4), 2,'omitnan')',mean(dwell_avg_control(:,:,4), 2)'
    ], idx)

violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]
violin(3).ViolinColor=[0.4 0.4 0.8]
violin(4).ViolinColor=[1 0.4 0.4]
violin(5).ViolinColor=[0.4 0.4 0.8]
violin(6).ViolinColor=[1 0.4 0.4]
violin(7).ViolinColor=[0.4 0.4 0.8]
violin(8).ViolinColor=[1 0.4 0.4]
xticklabels({'State 1','    ','State 2','     ','State 3', '','State 4'})
set(gca, 'FontSize', 15)
xticks(1:8)

ylabel('Dwell time (TRs)')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
[h, p(3), ci, stats]=ttest2(mean(dwell_avg_stroke(:,:,3),2, 'omitnan'),mean(dwell_avg_control(:,:,3),2))
[h, p(4), ci, stats]=ttest2(mean(dwell_avg_stroke(:,:,4),2, 'omitnan'),mean(dwell_avg_control(:,:,4),2))
[h, p(1), ci, stats]=ttest2(mean(dwell_avg_stroke(:,:,1),2, 'omitnan'),mean(dwell_avg_control(:,:,1), 2))
[h, p(2), ci, stats]=ttest2(mean(dwell_avg_stroke(:,:,2),2, 'omitnan'),mean(dwell_avg_control(:,:,2), 2))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end

old=0

for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.1, max(yt)*1.15, '*k');
    end
end



save(strcat(resdir, 'p_adjusted_Timeavg_DT_stroke_vs_control.mat'), 'p_adj')
save(strcat(resdir, 'p_unc_Timeavg_DR_stroke_vs_control.mat'), 'p')
saveas(gcf,strcat(figdir, 'DT_Timeavg_stroke_vs_control.png'))

%% Appearance rate
close all;
% mean over time
figure('Position', [0 0 900 300])
idx=[zeros(23, 1); ones(24,1); ...
    ones(23, 1)*2; ones(24,1)*3;...
    ones(23, 1)*4; ones(24,1)*5;...
    ones(23, 1)*6; ones(24,1)*7
    ]

violin=violinplot([mean(stroke_appearance1,2, 'omitnan')',mean(control_appearance1,2)',...
    mean(stroke_appearance2,2, 'omitnan')',mean(control_appearance2,2)',...
    mean(stroke_appearance3,2, 'omitnan')',mean(control_appearance3,2)',...
    mean(stroke_appearance4,2, 'omitnan')',mean(control_appearance4,2)'], idx)

ylabel('Appearance Rate')
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]
violin(3).ViolinColor=[0.4 0.4 0.8]
violin(4).ViolinColor=[1 0.4 0.4]
violin(5).ViolinColor=[0.4 0.4 0.8]
violin(6).ViolinColor=[1 0.4 0.4]
violin(7).ViolinColor=[0.4 0.4 0.8]
violin(8).ViolinColor=[1 0.4 0.4]

xticks(1:8)
xticklabels({'State 1','    ','State 2','     ','State 3', '','State 4'})
set(gca, 'FontSize', 15)
clear p
[h, p(1), ci, stats]=ttest2(mean(stroke_appearance1,2, 'omitnan'),mean(control_appearance1,2))
[h, p(2), ci, stats]=ttest2(mean(stroke_appearance2,2, 'omitnan'),mean(control_appearance2,2))
[h, p(3), ci, stats]=ttest2(mean(stroke_appearance3,2, 'omitnan'),mean(control_appearance3,2))
[h, p(4), ci, stats]=ttest2(mean(stroke_appearance4, 2,'omitnan'),mean(control_appearance4,2))
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end

old=0

for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.1, max(yt)*1.15, '*k');
    end
end


save(strcat(resdir, 'p_adjusted_Timeavg_AR_stroke_vs_control.mat'), 'p_adj')
save(strcat(resdir, 'p_unc_Timeavg_R_stroke_vs_control.mat'), 'p')
saveas(gcf,strcat(figdir, 'AR_Timeavg_stroke_vs_control.png'))

