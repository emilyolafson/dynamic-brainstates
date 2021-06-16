
figdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/figures/'
resdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/'


%% Fractional occupancy
idx=[ones(23,1); ones(23,1)*2;ones(23,1)*3;ones(23,1)*4;ones(23,1)*5];
tiledlayout(1,4,'padding','none')


nexttile;
viol=violinplot([stroke_FO1(:,1);stroke_FO1(:,2);stroke_FO1(:,3);stroke_FO1(:,4);stroke_FO1(:,5)], idx)
title('State 1')
ylim([0, 0.6])
viol(1).ViolinColor=[0.5 0.5 0.8]
viol(2).ViolinColor=[0.35 0.35 0.8]
viol(3).ViolinColor=[0.30 0.30 1]
viol(4).ViolinColor=[0.25 0.25 1]
viol(5).ViolinColor=[0.15 0.15 1]

ylabel('Fractional Occupancy')

nexttile;
viol=violinplot([stroke_FO2(:,1);stroke_FO2(:,2);stroke_FO2(:,3);stroke_FO2(:,4);stroke_FO2(:,5)], idx)
title('State 2')
ylim([0, 0.6])
viol(1).ViolinColor=[0.5 0.5 0.8]
viol(2).ViolinColor=[0.35 0.35 0.8]
viol(3).ViolinColor=[0.30 0.30 1]
viol(4).ViolinColor=[0.25 0.25 1]
viol(5).ViolinColor=[0.15 0.15 1]

ylabel('Fractional Occupancy')


nexttile;
viol=violinplot([stroke_FO3(:,1);stroke_FO3(:,2);stroke_FO3(:,3);stroke_FO3(:,4);stroke_FO3(:,5)], idx)
title('State 3')
ylim([0, 0.6])
viol(1).ViolinColor=[0.5 0.5 0.8]
viol(2).ViolinColor=[0.35 0.35 0.8]
viol(3).ViolinColor=[0.30 0.30 1]
viol(4).ViolinColor=[0.25 0.25 1]
viol(5).ViolinColor=[0.15 0.15 1]
ylabel('Fractional Occupancy')

nexttile;
viol=violinplot([stroke_FO4(:,1);stroke_FO4(:,2);stroke_FO4(:,3);stroke_FO4(:,4);stroke_FO4(:,5)], idx)
title('State 4')
ylim([0, 0.6])
viol(1).ViolinColor=[0.5 0.5 0.8]
viol(2).ViolinColor=[0.35 0.35 0.8]
viol(3).ViolinColor=[0.30 0.30 1]
viol(4).ViolinColor=[0.25 0.25 1]
viol(5).ViolinColor=[0.15 0.15 1]

ylabel('Fractional Occupancy')


%% Dwell time
idx=[ones(23,1); ones(23,1)*2;ones(23,1)*3;ones(23,1)*4;ones(23,1)*5];
tiledlayout(1,4,'padding','none')
for i=1:4
    nexttile;
    viol=violinplot([dwell_avg_stroke(:,1,i);dwell_avg_stroke(:,2,i);dwell_avg_stroke(:,3,i);dwell_avg_stroke(:,4,i);dwell_avg_stroke(:,5,i)], idx)
    title(['State ', num2str(i)])
    ylabel('Dwell Time')
    viol(1).ViolinColor=[0.5 0.5 0.8]
    viol(2).ViolinColor=[0.35 0.35 0.8]
    viol(3).ViolinColor=[0.30 0.30 1]
    viol(4).ViolinColor=[0.25 0.25 1]
    viol(5).ViolinColor=[0.15 0.15 1]
    p1=plot(mean([dwell_avg_stroke(:,1,i),dwell_avg_stroke(:,2,i),dwell_avg_stroke(:,3,i),dwell_avg_stroke(:,4,i),dwell_avg_stroke(:,5,i)], 'omitnan')','b-' ,'LineWidth',3)
    ylim([0 4])
end

%% Appearance rate

idx=[ones(23,1); ones(23,1)*2;ones(23,1)*3;ones(23,1)*4;ones(23,1)*5];
tiledlayout(1,4,'padding','none')

nexttile;
viol=violinplot([stroke_appearance1(:,1);stroke_appearance1(:,2);stroke_appearance1(:,3);stroke_appearance1(:,4);stroke_appearance1(:,5)], idx)
title('State 1')
viol(1).ViolinColor=[0.5 0.5 0.8]
viol(2).ViolinColor=[0.35 0.35 0.8]
viol(3).ViolinColor=[0.30 0.30 1]
viol(4).ViolinColor=[0.25 0.25 1]
viol(5).ViolinColor=[0.15 0.15 1]
plot(mean([stroke_appearance1(:,1),stroke_appearance1(:,2),stroke_appearance1(:,3),stroke_appearance1(:,4),stroke_appearance1(:,5)], 'omitnan')', '-b', 'LineWidth', 3)
ylim([0 4.5])

ylabel('Appearance Rate')

nexttile;
viol=violinplot([stroke_appearance2(:,1);stroke_appearance2(:,2);stroke_appearance2(:,3);stroke_appearance2(:,4);stroke_appearance2(:,5)], idx)
title('State 2')
viol(1).ViolinColor=[0.5 0.5 0.8]
viol(2).ViolinColor=[0.35 0.35 0.8]
viol(3).ViolinColor=[0.30 0.30 1]
viol(4).ViolinColor=[0.25 0.25 1]
viol(5).ViolinColor=[0.15 0.15 1]
plot(mean([stroke_appearance2(:,1),stroke_appearance2(:,2),stroke_appearance2(:,3),stroke_appearance2(:,4),stroke_appearance2(:,5)], 'omitnan')', '-b', 'LineWidth', 3)
ylim([0 4.5])

ylabel('Appearance Rate')


nexttile;
viol=violinplot([stroke_appearance3(:,1);stroke_appearance3(:,2);stroke_appearance3(:,3);stroke_appearance3(:,4);stroke_appearance3(:,5)], idx)
title('State 3')
viol(1).ViolinColor=[0.5 0.5 0.8]
viol(2).ViolinColor=[0.35 0.35 0.8]
viol(3).ViolinColor=[0.30 0.30 1]
viol(4).ViolinColor=[0.25 0.25 1]
viol(5).ViolinColor=[0.15 0.15 1]
plot(mean([stroke_appearance3(:,1),stroke_appearance3(:,2),stroke_appearance3(:,3),stroke_appearance3(:,4),stroke_appearance3(:,5)], 'omitnan')', '-b', 'LineWidth', 3)
ylim([0 4.5])

ylabel('Appearance Rate')

nexttile;
viol=violinplot([stroke_appearance4(:,1);stroke_appearance4(:,2);stroke_appearance4(:,3);stroke_appearance4(:,4);stroke_appearance4(:,5)], idx)
title('State 4')
viol(1).ViolinColor=[0.5 0.5 0.8]
viol(2).ViolinColor=[0.35 0.35 0.8]
viol(3).ViolinColor=[0.30 0.30 1]
viol(4).ViolinColor=[0.25 0.25 1]
viol(5).ViolinColor=[0.15 0.15 1]
plot(mean([stroke_appearance4(:,1),stroke_appearance4(:,2),stroke_appearance4(:,3),stroke_appearance4(:,4),stroke_appearance4(:,5)], 'omitnan')', '-b', 'LineWidth', 3)
ylabel('Appearance Rate')
ylim([0 4.5])
