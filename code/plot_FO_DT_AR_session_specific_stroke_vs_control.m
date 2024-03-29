
figdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/shen268/k2/'

%% calculate significant differences between stroke and controls in FO at each itme point.
clear p
[h,p(1,1),~,~]=ttest2(stroke_FO1(:,1),control_FO1(:,1))
[h,p(2,1),~,~]=ttest2(stroke_FO2(:,1),control_FO2(:,1))
[h,p(3,1),~,~]=ttest2(stroke_FO3(:,1),control_FO3(:,1))
[h,p(4,1),~,~]=ttest2(stroke_FO4(:,1),control_FO4(:,1))

[h,p(1,2),~,~]=ttest2(stroke_FO1(:,2),control_FO1(:,2))
[h,p(2,2),~,~]=ttest2(stroke_FO2(:,2),control_FO2(:,2))
[h,p(3,2),~,~]=ttest2(stroke_FO3(:,2),control_FO3(:,2))
[h,p(4,2),~,~]=ttest2(stroke_FO4(:,2),control_FO4(:,2))

[h,p(1,3),~,~]=ttest2(stroke_FO1(:,3),control_FO1(:,3))
[h,p(2,3),~,~]=ttest2(stroke_FO2(:,3),control_FO2(:,3))
[h,p(3,3),~,~]=ttest2(stroke_FO3(:,3),control_FO3(:,3))
[h,p(4,3),~,~]=ttest2(stroke_FO4(:,3),control_FO4(:,3))

[h,p(1,4),~,~]=ttest2(stroke_FO1(:,4),control_FO1(:,4))
[h,p(2,4),~,~]=ttest2(stroke_FO2(:,4),control_FO2(:,4))
[h,p(3,4),~,~]=ttest2(stroke_FO3(:,4),control_FO3(:,4))
[h,p(4,4),~,~]=ttest2(stroke_FO4(:,4),control_FO4(:,4))

[h,p(1,5),~,~]=ttest2(stroke_FO1(:,5),control_FO1(:,5))
[h,p(2,5),~,~]=ttest2(stroke_FO2(:,5),control_FO2(:,5))
[h,p(3,5),~,~]=ttest2(stroke_FO3(:,5),control_FO3(:,5))
[h,p(4,5),~,~]=ttest2(stroke_FO4(:,5),control_FO4(:,5))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

%% Plot difference between stroke & control subjects at each time point
%% Fractional Occupancy

figure('Position', [0 0 500 1500])
tiledlayout(5, 1, 'padding', 'none')
nexttile;

idx=[zeros(23, 1); ones(24,1); ...
    ones(23, 1)*2; ones(24,1)*3;...
    ones(23, 1)*4; ones(24,1)*5;...
    ones(23, 1)*6; ones(24,1)*7
    ]
s=1
violin=violinplot([stroke_FO1(:,s);control_FO1(:,s);...
    stroke_FO2(:,s);control_FO2(:,s);...
    stroke_FO3(:,s);control_FO3(:,s);...
    stroke_FO4(:,s);control_FO4(:,s)
    ], idx)
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
ylabel('Fractional Occupancy')
title('Session 1')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ylim([0 .6])

old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end
old=0

for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end

nexttile;
s=2;
violin=violinplot([stroke_FO1(:,s);control_FO1(:,s);...
    stroke_FO2(:,s);control_FO2(:,s);...
    stroke_FO3(:,s);control_FO3(:,s);...
    stroke_FO4(:,s);control_FO4(:,s)
    ], idx)
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
ylabel('Fractional Occupancy')
title('Session 2')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ylim([0 .6])


old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end

nexttile;
s=3;
violin=violinplot([stroke_FO1(:,s);control_FO1(:,s);...
    stroke_FO2(:,s);control_FO2(:,s);...
    stroke_FO3(:,s);control_FO3(:,s);...
    stroke_FO4(:,s);control_FO4(:,s)
    ], idx)
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

ylabel('Fractional Occupancy')
title('Session 3')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ylim([0 .6])

old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end


nexttile;
s=4;
violin=violinplot([stroke_FO1(:,s);control_FO1(:,s);...
    stroke_FO2(:,s);control_FO2(:,s);...
    stroke_FO3(:,s);control_FO3(:,s);...
    stroke_FO4(:,s);control_FO4(:,s)
    ], idx)
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
ylabel('Fractional Occupancy')
title('Session 4')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ylim([0 .6])


old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end


nexttile;
s=5;
violin=violinplot([stroke_FO1(:,s);control_FO1(:,s);...
    stroke_FO2(:,s);control_FO2(:,s);...
    stroke_FO3(:,s);control_FO3(:,s);...
    stroke_FO4(:,s);control_FO4(:,s)
    ], idx)
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
ylabel('Fractional Occupancy')
title('Session 5')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ylim([0 .6])

old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end

saveas(gcf,strcat(figdir, 'FO_allsessions_stroke_vs_control.png'))

%% Dwell time 
clear p
for i=1:4
    for j=1:5
    [~, p(i, j), ~,~]=ttest2(dwell_avg_stroke(:,j,i),dwell_avg_control(:,j,i))
    end
end

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

figure('Position', [0 0 500 1500])
tiledlayout(5, 1, 'padding', 'none')
nexttile;

idx=[zeros(23, 1); ones(24,1); ...
    ones(23, 1)*2; ones(24,1)*3;...
    ones(23, 1)*4; ones(24,1)*5;...
    ones(23, 1)*6; ones(24,1)*7
    ]
s=1
violin=violinplot([dwell_avg_stroke(:,s,1);dwell_avg_control(:,s,1);...
    dwell_avg_stroke(:,s,2);dwell_avg_control(:,s,2);...
    dwell_avg_stroke(:,s,3);dwell_avg_control(:,s,3);...
    dwell_avg_stroke(:,s,4);dwell_avg_control(:,s,4)
    ], idx)
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
ylabel('Dwell Time (s)')
title('Session 1')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ylim([1 5])

old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end

old=0

for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end

nexttile;
s=2;
violin=violinplot([dwell_avg_stroke(:,s,1);dwell_avg_control(:,s,1);...
    dwell_avg_stroke(:,s,2);dwell_avg_control(:,s,2);...
    dwell_avg_stroke(:,s,3);dwell_avg_control(:,s,3);...
    dwell_avg_stroke(:,s,4);dwell_avg_control(:,s,4)
    ], idx)
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
ylabel('Dwell Time (s)')
title('Session 2')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ylim([1 5])


old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end

nexttile;
s=3;
violin=violinplot([dwell_avg_stroke(:,s,1);dwell_avg_control(:,s,1);...
    dwell_avg_stroke(:,s,2);dwell_avg_control(:,s,2);...
    dwell_avg_stroke(:,s,3);dwell_avg_control(:,s,3);...
    dwell_avg_stroke(:,s,4);dwell_avg_control(:,s,4)
    ], idx)
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

ylabel('Dwell Time (s)')
title('Session 3')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ylim([1 5])

old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end


nexttile;
s=4;
violin=violinplot([dwell_avg_stroke(:,s,1);dwell_avg_control(:,s,1);...
    dwell_avg_stroke(:,s,2);dwell_avg_control(:,s,2);...
    dwell_avg_stroke(:,s,3);dwell_avg_control(:,s,3);...
    dwell_avg_stroke(:,s,4);dwell_avg_control(:,s,4)
    ], idx)
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
ylabel('Dwell Time (s)')
title('Session 4')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ylim([1 5])


old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end


nexttile;
s=5;
violin=violinplot([dwell_avg_stroke(:,s,1);dwell_avg_control(:,s,1);...
    dwell_avg_stroke(:,s,2);dwell_avg_control(:,s,2);...
    dwell_avg_stroke(:,s,3);dwell_avg_control(:,s,3);...
    dwell_avg_stroke(:,s,4);dwell_avg_control(:,s,4)
    ], idx)
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
ylabel('Dwell Time (s)')
title('Session 5')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ylim([1 5])

old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end

saveas(gcf,strcat(figdir, 'DT_allsessions_stroke_vs_control.png'))

%% Appearance rate - session-specific comparison
for i=1:5
    [h,p(1,i),~,~]=ttest2(stroke_appearance1(:,i),control_appearance1(:,i))
    [h,p(2,i),~,~]=ttest2(stroke_appearance2(:,i),control_appearance2(:,i))
    [h,p(3,i),~,~]=ttest2(stroke_appearance3(:,i),control_appearance3(:,i))
    [h,p(4,i),~,~]=ttest2(stroke_appearance4(:,i),control_appearance4(:,i))
end

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

figure('Position', [0 0 500 1500])
tiledlayout(5, 1, 'padding', 'none')
nexttile;

idx=[zeros(23, 1); ones(24,1); ...
    ones(23, 1)*2; ones(24,1)*3;...
    ones(23, 1)*4; ones(24,1)*5;...
    ones(23, 1)*6; ones(24,1)*7
    ]
s=1
violin=violinplot([stroke_appearance1(:,s);control_appearance1(:,s);...
    stroke_appearance2(:,s);control_appearance2(:,s);...
    stroke_appearance3(:,s);control_appearance3(:,s);...
    stroke_appearance4(:,s);control_appearance4(:,s)
    ], idx)
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
ylabel('Appearance Rate')
title('Session 1')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ylim([1 5])
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end
old=0

for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end

nexttile;
s=2;
violin=violinplot([stroke_appearance1(:,s);control_appearance1(:,s);...
    stroke_appearance2(:,s);control_appearance2(:,s);...
    stroke_appearance3(:,s);control_appearance3(:,s);...
    stroke_appearance4(:,s);control_appearance4(:,s)
    ], idx)
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
ylabel('Appearance Rate')
title('Session 2')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ylim([1 5])


old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end

nexttile;
s=3;
violin=violinplot([stroke_appearance1(:,s);control_appearance1(:,s);...
    stroke_appearance2(:,s);control_appearance2(:,s);...
    stroke_appearance3(:,s);control_appearance3(:,s);...
    stroke_appearance4(:,s);control_appearance4(:,s)
    ], idx)
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]
violin(3).ViolinColor=[0.4 0.4 0.8]
violin(4).ViolinColor=[1 0.4 0.4]
violin(5).ViolinColor=[0.4 0.4 0.8]
violin(6).ViolinColor=[1 0.4 0.4]
violin(7).ViolinColor=[0.4 0.4 0.8]
violin(8).ViolinColor=[1 0.4 0.4]
xticks(1:8)
ylim([1 5])

xticklabels({'State 1','    ','State 2','     ','State 3', '','State 4'})
set(gca, 'FontSize', 15)

ylabel('Appearance Rate')
title('Session 3')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');

old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end


nexttile;
s=4;
violin=violinplot([stroke_appearance1(:,s);control_appearance1(:,s);...
    stroke_appearance2(:,s);control_appearance2(:,s);...
    stroke_appearance3(:,s);control_appearance3(:,s);...
    stroke_appearance4(:,s);control_appearance4(:,s)
    ], idx)
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
ylabel('Appearance Rate')
title('Session 4')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ylim([1 5])


old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end


nexttile;
s=5;
violin=violinplot([stroke_appearance1(:,s);control_appearance1(:,s);...
    stroke_appearance2(:,s);control_appearance2(:,s);...
    stroke_appearance3(:,s);control_appearance3(:,s);...
    stroke_appearance4(:,s);control_appearance4(:,s)
    ], idx)
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
ylabel('Appearance Rate')
title('Session 5')
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ylim([1 5])

old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p(i,s)<0.05
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k')
    end
end
old=0
for i=1:4
    xtickindex=i+old
    old=i;
    signif=p_adj(i,s)<0.05;
    if signif
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1])), max(yt)*1.15, '*k');
        plot(xt([xtickindex xtickindex+1]), [1 1]*max(yt)*1.1, '-k',  mean(xt([xtickindex xtickindex+1]))+0.2, max(yt)*1.15, '*k');
    end
end

saveas(gcf,strcat(figdir, 'AR_allsessions_stroke_vs_control.png'))
close all;
