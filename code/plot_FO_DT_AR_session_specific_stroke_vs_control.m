%% Plot difference between stroke & control subjects at each time point

%%Fractional occupancy - session-specific comparison
% compared to controls
ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_FO1(:,1);control_FO1(:,1)],ids)
nexttile;
violinplot([stroke_FO2(:,1);control_FO2(:,1)],ids)
nexttile;
violinplot([stroke_FO3(:,1);control_FO3(:,1)],ids)
nexttile;
violinplot([stroke_FO4(:,1);control_FO4(:,1)],ids)

[h,p(1),~,~]=ttest2(stroke_FO1(:,1),control_FO1(:,1))
[h,p(2),~,~]=ttest2(stroke_FO2(:,1),control_FO2(:,1))
[h,p(3),~,~]=ttest2(stroke_FO3(:,1),control_FO3(:,1))
[h,p(4),~,~]=ttest2(stroke_FO4(:,1),control_FO4(:,1))

% plot session w sig diffs
close all;
figure()
ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violin=violinplot([stroke_FO1(:,2);control_FO1(:,2)],ids)
title('State 1')
xticklabels({"Stroke", "Control"})
set(gca, 'FontSize', 15)
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]
nexttile;
violin=violinplot([stroke_FO2(:,2);control_FO2(:,2)],ids)
title('State 2')
xticklabels({"Stroke", "Control"})
set(gca, 'FontSize', 15)
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]

nexttile;
violin=violinplot([stroke_FO3(:,2);control_FO3(:,2)],ids)
title('State 3')
xticklabels({"Stroke", "Control"})
set(gca, 'FontSize', 15)
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]

nexttile;
violin=violinplot([stroke_FO4(:,2);control_FO4(:,2)],ids)
title('State 4')
xticklabels({"Stroke", "Control"})
sgtitle('Fractial occupancy - 2 weeks post-stroke')
set(gca, 'FontSize', 15)
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]

[h,p(5),~,~]=ttest2(stroke_FO1(:,2),control_FO1(:,2))
[h,p(6),~,~]=ttest2(stroke_FO2(:,2),control_FO2(:,2))
[h,p(7),~,~]=ttest2(stroke_FO3(:,2),control_FO3(:,2))
[h,p(8),~,~]=ttest2(stroke_FO4(:,2),control_FO4(:,2))


figure()
ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_FO1(:,3);control_FO1(:,3)],ids)
nexttile;
violinplot([stroke_FO2(:,3);control_FO2(:,3)],ids)
nexttile;
violinplot([stroke_FO3(:,3);control_FO3(:,3)],ids)
nexttile;
violinplot([stroke_FO4(:,3);control_FO4(:,3)],ids)
[h,p(9),~,~]=ttest2(stroke_FO1(:,3),control_FO1(:,3))
[h,p(10),~,~]=ttest2(stroke_FO2(:,3),control_FO2(:,3))
[h,p(11),~,~]=ttest2(stroke_FO3(:,3),control_FO3(:,3))
[h,p(12),~,~]=ttest2(stroke_FO4(:,3),control_FO4(:,3))


ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_FO1(:,4);control_FO1(:,4)],ids)
nexttile;
violinplot([stroke_FO2(:,4);control_FO2(:,4)],ids)
nexttile;
violinplot([stroke_FO3(:,4);control_FO3(:,4)],ids)
nexttile;
violinplot([stroke_FO4(:,4);control_FO4(:,4)],ids)

[h,p(13),~,~]=ttest2(stroke_FO1(:,4),control_FO1(:,4))
[h,p(14),~,~]=ttest2(stroke_FO2(:,4),control_FO2(:,4))
[h,p(15),~,~]=ttest2(stroke_FO3(:,4),control_FO3(:,4))
[h,p(16),~,~]=ttest2(stroke_FO4(:,4),control_FO4(:,4))


ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_FO1(:,5);control_FO1(:,5)],ids)
nexttile;
violinplot([stroke_FO2(:,5);control_FO2(:,5)],ids)
nexttile;
violinplot([stroke_FO3(:,5);control_FO3(:,5)],ids)
nexttile;
violinplot([stroke_FO4(:,5);control_FO4(:,5)],ids)

[h,p(17),~,~]=ttest2(stroke_FO1(:,5),control_FO1(:,5))
[h,p(18),~,~]=ttest2(stroke_FO2(:,5),control_FO2(:,5))
[h,p(19),~,~]=ttest2(stroke_FO3(:,5),control_FO3(:,5))
[h,p(20),~,~]=ttest2(stroke_FO4(:,5),control_FO4(:,5))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')


%% Dwell time - specific sessions
% stroke vs control

idx=[ones(23,1);zeros(24,1)]
% session 1 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,1,i);dwell_avg_control(:,1,i)], idx)
    [h, p(i), ci, stats]=ttest2(dwell_avg_stroke(:,1,i),dwell_avg_control(:,1,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end

% session 2 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,2,i);dwell_avg_control(:,2,i)], idx)
    [h, p(i+4), ci, stats]=ttest2(dwell_avg_stroke(:,2,i),dwell_avg_control(:,2,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end

% session 3 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,3,i);dwell_avg_control(:,3,i)], idx)
    [h, p(i+8), ci, stats]=ttest2(dwell_avg_stroke(:,3,i),dwell_avg_control(:,3,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end

% session 4 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,4,i);dwell_avg_control(:,4,i)], idx)
    [h, p(i+12), ci, stats]=ttest2(dwell_avg_stroke(:,4,i),dwell_avg_control(:,4,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end

% session 5 stroke; state 1-4
tiledlayout(2,2,'padding', 'none')
for i=1:4
    disp(i)
    nexttile
    violinplot([dwell_avg_stroke(:,5,i);dwell_avg_control(:,5,i)], idx)
    [h, p(i+16), ci, stats]=ttest2(dwell_avg_stroke(:,5,i),dwell_avg_control(:,5,i))
    title(['state: ', num2str(i)], ['p = ', num2str(p), ', tstat = ', num2str(stats.tstat)])
end

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

%% Appearance rate - session-specific comparison
% compared to controls
ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_appearance1(:,1);control_appearance1(:,1)],ids)
nexttile;
violinplot([stroke_appearance2(:,1);control_appearance2(:,1)],ids)
nexttile;
violinplot([stroke_appearance3(:,1);control_appearance3(:,1)],ids)
nexttile;
violinplot([stroke_appearance4(:,1);control_appearance4(:,1)],ids)

[h,p(1),~,~]=ttest2(stroke_appearance1(:,1),control_appearance1(:,1))
[h,p(2),~,~]=ttest2(stroke_appearance2(:,1),control_appearance2(:,1))
[h,p(3),~,~]=ttest2(stroke_appearance3(:,1),control_appearance3(:,1))
[h,p(4),~,~]=ttest2(stroke_appearance4(:,1),control_appearance4(:,1))

% plot session w sig diffs
close all;
figure()
ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violin=violinplot([stroke_appearance1(:,2);control_appearance1(:,2)],ids)
title('State 1')
xticklabels({"Stroke", "Control"})
set(gca, 'FontSize', 15)
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]
nexttile;
violin=violinplot([stroke_appearance2(:,2);control_appearance2(:,2)],ids)
title('State 2')
xticklabels({"Stroke", "Control"})
set(gca, 'FontSize', 15)
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]

nexttile;
violin=violinplot([stroke_appearance3(:,2);control_appearance3(:,2)],ids)
title('State 3')
xticklabels({"Stroke", "Control"})
set(gca, 'FontSize', 15)
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]

nexttile;
violin=violinplot([stroke_appearance4(:,2);control_appearance4(:,2)],ids)
title('State 4')
xticklabels({"Stroke", "Control"})
sgtitle('Fractial occupancy - 2 weeks post-stroke')
set(gca, 'FontSize', 15)
violin(1).ViolinColor=[0.4 0.4 0.8]
violin(2).ViolinColor=[1 0.4 0.4]

[h,p(5),~,~]=ttest2(stroke_appearance1(:,2),control_appearance1(:,2))
[h,p(6),~,~]=ttest2(stroke_appearance2(:,2),control_appearance2(:,2))
[h,p(7),~,~]=ttest2(stroke_appearance3(:,2),control_appearance3(:,2))
[h,p(8),~,~]=ttest2(stroke_appearance4(:,2),control_appearance4(:,2))


figure()
ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_appearance1(:,3);control_appearance1(:,3)],ids)
nexttile;
violinplot([stroke_appearance2(:,3);control_appearance2(:,3)],ids)
nexttile;
violinplot([stroke_appearance3(:,3);control_appearance3(:,3)],ids)
nexttile;
violinplot([stroke_appearance4(:,3);control_appearance4(:,3)],ids)
[h,p(9),~,~]=ttest2(stroke_appearance1(:,3),control_appearance1(:,3))
[h,p(10),~,~]=ttest2(stroke_appearance2(:,3),control_appearance2(:,3))
[h,p(11),~,~]=ttest2(stroke_appearance3(:,3),control_appearance3(:,3))
[h,p(12),~,~]=ttest2(stroke_appearance4(:,3),control_appearance4(:,3))


ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_appearance1(:,4);control_appearance1(:,4)],ids)
nexttile;
violinplot([stroke_appearance2(:,4);control_appearance2(:,4)],ids)
nexttile;
violinplot([stroke_appearance3(:,4);control_appearance3(:,4)],ids)
nexttile;
violinplot([stroke_appearance4(:,4);control_appearance4(:,4)],ids)

[h,p(13),~,~]=ttest2(stroke_appearance1(:,4),control_appearance1(:,4))
[h,p(14),~,~]=ttest2(stroke_appearance2(:,4),control_appearance2(:,4))
[h,p(15),~,~]=ttest2(stroke_appearance3(:,4),control_appearance3(:,4))
[h,p(16),~,~]=ttest2(stroke_appearance4(:,4),control_appearance4(:,4))


ids=[zeros(1,23), ones(1,24)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([stroke_appearance1(:,5);control_appearance1(:,5)],ids)
nexttile;
violinplot([stroke_appearance2(:,5);control_appearance2(:,5)],ids)
nexttile;
violinplot([stroke_appearance3(:,5);control_appearance3(:,5)],ids)
nexttile;
violinplot([stroke_appearance4(:,5);control_appearance4(:,5)],ids)

[h,p(17),~,~]=ttest2(stroke_appearance1(:,5),control_appearance1(:,5))
[h,p(18),~,~]=ttest2(stroke_appearance2(:,5),control_appearance2(:,5))
[h,p(19),~,~]=ttest2(stroke_appearance3(:,5),control_appearance3(:,5))
[h,p(20),~,~]=ttest2(stroke_appearance4(:,5),control_appearance4(:,5))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')


