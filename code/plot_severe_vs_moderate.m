%% Severe vs moderate stroke subjects
idx_severe=baselineFM < 50
idx_mod=baselineFM >= 50
sev=sum(idx_severe)
mod=sum(idx_mod)


%% Fractional occupancy - moderate to severe
% across all sessions
severe_FO1=stroke_FO1(idx_severe,:)
mod_FO1=stroke_FO1(idx_mod,:)
severe_FO2=stroke_FO2(idx_severe,:)
mod_FO2=stroke_FO2(idx_mod,:)
severe_FO3=stroke_FO3(idx_severe,:)
mod_FO3=stroke_FO3(idx_mod,:)
severe_FO4=stroke_FO4(idx_severe,:)
mod_FO4=stroke_FO4(idx_mod,:)

ids=[zeros(1,5), ones(1,18)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([mean(severe_FO1(:,:),2, 'omitnan');mean(mod_FO1(:,:),2, 'omitnan')],ids)
nexttile;
violinplot([mean(severe_FO2(:,:),2, 'omitnan');mean(mod_FO2(:,:), 2, 'omitnan')],ids)
nexttile;
violinplot([mean(severe_FO3(:,:),2, 'omitnan');mean(mod_FO3(:,:), 2, 'omitnan')],ids)
nexttile;
violinplot([mean(severe_FO4(:,:),2, 'omitnan');mean(mod_FO4(:,:), 2, 'omitnan')],ids)
set(gcf, 'Position', [0 0 1000 1000])

[h, p(1), ~, ~]=ttest2(mean(severe_FO1(:,:), 2, 'omitnan'), mean(mod_FO1(:,:), 2, 'omitnan'))
[h, p(2), ~, ~]=ttest2(mean(severe_FO2(:,:),2, 'omitnan'), mean(mod_FO2(:,:),2, 'omitnan'))
[h, p(3), ~, ~]=ttest2(mean(severe_FO3(:,:),2, 'omitnan'), mean(mod_FO3(:,:),2, 'omitnan'))
[h, p(4), ~, ~]=ttest2(mean(severe_FO4(:,:),2, 'omitnan'), mean(mod_FO4(:,:), 2, 'omitnan'))



%% Dwell time - moderate vs severe
% across single sessions
severe_FO1=dwell_avg_stroke(idx_severe,:,1)
mod_FO1=dwell_avg_stroke(idx_mod,:,1)
severe_FO2=dwell_avg_stroke(idx_severe,:,2)
mod_FO2=dwell_avg_stroke(idx_mod,:,2)
severe_FO3=dwell_avg_stroke(idx_severe,:,3)
mod_FO3=dwell_avg_stroke(idx_mod,:,3)
severe_FO4=dwell_avg_stroke(idx_severe,:,4)
mod_FO4=dwell_avg_stroke(idx_mod,:,4)

%across all sessions
ids=[zeros(1,11), ones(1,12)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([reshape(severe_FO1(:,:),[],1);reshape(mod_FO1(:,:), [],1)],ids)
nexttile;
violinplot([reshape(severe_FO2(:,:),[],1);reshape(mod_FO2(:,:), [],1)],ids)
nexttile;
violinplot([reshape(severe_FO3(:,:),[],1);reshape(mod_FO3(:,:), [],1)],ids)
nexttile;
violinplot([reshape(severe_FO4(:,:),[],1);reshape(mod_FO4(:,:), [],1)],ids)

[h, p(1), ~, ~]=ttest2(reshape(severe_FO1(:,:),[],1), reshape(mod_FO1(:,:), [],1))
[h, p(2), ~, ~]=ttest2(reshape(severe_FO2(:,:),[],1), reshape(mod_FO2(:,:), [],1))
[h, p(3), ~, ~]=ttest2(reshape(severe_FO3(:,:),[],1), reshape(mod_FO3(:,:), [],1))
[h, p(4), ~, ~]=ttest2(reshape(severe_FO4(:,:),[],1), reshape(mod_FO4(:,:), [],1))

%mean across all sessions
ids=[zeros(1,11), ones(1,12)]
tiledlayout(2,2,'padding', 'none')
nexttile;
violinplot([mean(severe_FO1(:,:),2, 'omitnan');mean(mod_FO1(:,:),2, 'omitnan')],ids)
nexttile;
violinplot([mean(severe_FO2(:,:),2, 'omitnan');mean(mod_FO2(:,:), 2, 'omitnan')],ids)
nexttile;
violinplot([mean(severe_FO3(:,:),2, 'omitnan');mean(mod_FO3(:,:), 2, 'omitnan')],ids)
nexttile;
violinplot([mean(severe_FO4(:,:),2, 'omitnan');mean(mod_FO4(:,:), 2, 'omitnan')],ids)

[h, p(1), ~, ~]=ttest2(mean(severe_FO1(:,:), 2, 'omitnan'), mean(mod_FO1(:,:), 2, 'omitnan'))
[h, p(2), ~, ~]=ttest2(mean(severe_FO2(:,:),2, 'omitnan'), mean(mod_FO2(:,:),2, 'omitnan'))
[h, p(3), ~, ~]=ttest2(mean(severe_FO3(:,:),2, 'omitnan'), mean(mod_FO3(:,:),2, 'omitnan'))
[h, p(4), ~, ~]=ttest2(mean(severe_FO4(:,:),2, 'omitnan'), mean(mod_FO4(:,:), 2, 'omitnan'))



%% Appearance rate - moderate vs severe
% across single sessions
severe_FO1=stroke_appearance1(idx_severe,:)
mod_FO1=stroke_appearance1(idx_mod,:)
severe_FO2=stroke_appearance2(idx_severe,:)
mod_FO2=stroke_appearance2(idx_mod,:)
severe_FO3=stroke_appearance3(idx_severe,:)
mod_FO3=stroke_appearance3(idx_mod,:)
severe_FO4=stroke_appearance4(idx_severe,:)
mod_FO4=stroke_appearance4(idx_mod,:)

ids=[zeros(1,sev), ones(1,mod)]
tiledlayout(1,4,'padding', 'none')
nexttile;
violinplot([severe_FO1(:,1);mod_FO1(:,1)],ids)
nexttile;
violinplot([severe_FO2(:,1);mod_FO2(:,1)],ids)
nexttile;
violinplot([severe_FO3(:,1);mod_FO3(:,1)],ids)
nexttile;
violinplot([severe_FO4(:,1);mod_FO4(:,1)],ids)

[h, p, ~, ~]=ttest2(severe_FO1(:,1), mod_FO1(:,1))
[h, p, ~, ~]=ttest2(severe_FO2(:,1), mod_FO2(:,1))
[h, p, ~, ~]=ttest2(severe_FO3(:,1), mod_FO3(:,1))
[h, p, ~, ~]=ttest2(severe_FO4(:,1), mod_FO4(:,1))

% across all sessions
ids=[zeros(1,sev), ones(1,mod)]
tiledlayout(1,4,'padding', 'none')
nexttile;
viol=violinplot([mean(severe_FO1(:,:),2, 'omitnan');mean(mod_FO1(:,:), 2, 'omitnan')],ids)
viol(1).ViolinColor=[0.4 0.4 0.8]
viol(2).ViolinColor=[1 0.4 0.4]
ylim([1, 5])
title('State 1')
xticklabels({"Severe", "Moderate"})
nexttile;

viol=violinplot([mean(severe_FO2(:,:),2, 'omitnan');mean(mod_FO2(:,:),2, 'omitnan')],ids)
viol(1).ViolinColor=[0.4 0.4 0.8]
viol(2).ViolinColor=[1 0.4 0.4]
ylim([1, 5])
title('State 2')
xticklabels({"Severe", "Moderate"})

nexttile;

viol=violinplot([mean(severe_FO3(:,:),2, 'omitnan');mean(mod_FO3(:,:), 2, 'omitnan')],ids)
viol(1).ViolinColor=[0.4 0.4 0.8]
viol(2).ViolinColor=[1 0.4 0.4]
ylim([1, 5])
title('State 3')
xticklabels({"Severe", "Moderate"})

nexttile;

viol=violinplot([mean(severe_FO4(:,:),2, 'omitnan');mean(mod_FO4(:,:) ,2, 'omitnan')],ids)
viol(1).ViolinColor=[0.4 0.4 0.8]
viol(2).ViolinColor=[1 0.4 0.4]
ylim([1, 5])
title('State 4')
xticklabels({"Severe", "Moderate"})

set(gcf, 'Position', [0 0 1000 1000])

[h, p(1), ~, ~]=ttest2(mean(severe_FO1(:,:),2, 'omitnan'), mean(mod_FO1(:,:), 2, 'omitnan'))
[h, p(2), ~, ~]=ttest2(mean(severe_FO2(:,:),2, 'omitnan'), mean(mod_FO2(:,:), 2, 'omitnan'))
[h, p(3), ~, ~]=ttest2(mean(severe_FO3(:,:),2, 'omitnan'), mean(mod_FO3(:,:), 2, 'omitnan'))
[h, p(4), ~, ~]=ttest2(mean(severe_FO4(:,:),2, 'omitnan'), mean(mod_FO4(:,:), 2, 'omitnan'))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')
