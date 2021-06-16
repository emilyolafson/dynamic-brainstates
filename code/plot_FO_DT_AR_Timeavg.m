%% Analysis - Fractional occupancy and dwell time between stroke & control
% FO differences across all sessions, stroke and controls;
% mean over time
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

xticks(1:2:12)
xticklabels({'State 1','State 2','State 3', 'State 4'})
set(gca, 'FontSize', 15)

[h, p(1), ci, stats]=ttest2(mean(stroke_FO1, 2, 'omitnan'),mean(control_FO1, 2))
[h, p(2), ci, stats]=ttest2(mean(stroke_FO2,  2, 'omitnan'),mean(control_FO2, 2))
[h, p(3), ci, stats]=ttest2(mean(stroke_FO3,  2, 'omitnan'),mean(control_FO3, 2))
[h, p(4), ci, stats]=ttest2(mean(stroke_FO4,  2, 'omitnan'),mean(control_FO4, 2))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

%% Dwell time differences across all sessions, stroke and controls;

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
xticks(1:2:12)
xticklabels({'State 1','State 2','State 3', 'State 4'})
set(gca, 'FontSize', 15)
ylabel('Dwell time (TRs)')

[h, p(3), ci, stats]=ttest2(mean(dwell_avg_stroke(:,:,3),2, 'omitnan'),mean(dwell_avg_control(:,:,3),2))
[h, p(4), ci, stats]=ttest2(mean(dwell_avg_stroke(:,:,4),2, 'omitnan'),mean(dwell_avg_control(:,:,4),2))
[h, p(1), ci, stats]=ttest2(mean(dwell_avg_stroke(:,:,1),2, 'omitnan'),mean(dwell_avg_control(:,:,1), 2))
[h, p(2), ci, stats]=ttest2(mean(dwell_avg_stroke(:,:,2),2, 'omitnan'),mean(dwell_avg_control(:,:,2), 2))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')


%% Appearance rate

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

xticks(1:2:12)
xticklabels({'State 1','State 2','State 3', 'State 4'})
set(gca, 'FontSize', 15)
clear p
[h, p(1), ci, stats]=ttest2(mean(stroke_appearance1,2, 'omitnan'),mean(control_appearance1,2))
[h, p(2), ci, stats]=ttest2(mean(stroke_appearance2,2, 'omitnan'),mean(control_appearance2,2))
[h, p(3), ci, stats]=ttest2(mean(stroke_appearance3,2, 'omitnan'),mean(control_appearance3,2))
[h, p(4), ci, stats]=ttest2(mean(stroke_appearance4, 2,'omitnan'),mean(control_appearance4,2))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')
