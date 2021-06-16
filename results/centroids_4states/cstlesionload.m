% CST lesion load and state 3 dynamics 

cst=load('/Users/emilyolafson/GIT/stroke-graph-matching/project/data/CSTlesionload.mat')
cst=cst.lesionload;
cst=cell2mat(cst)

cst2=load('/Users/emilyolafson/GIT/dynamic-brainstates/lesionload_CSTR_CSTL.mat')
cst2=cst2.ll
rightlesions=[1;1;0;0;1;0;1;1;1;1;0;0;0;1;1;0;0;1;1;1;1;0;1]

cstright=cst2(:,1)
cstleft=cst2(:,2)

corr(cst, fm_1, 'rows', 'complete')
plot(cst, fm_1)

idx_severe=cst > 0.01
idx_mod=cst <=  0.01
sev=sum(idx_severe)
mod=sum(idx_mod)

idx_severe=cstleft > 0.001
idx_mod=cstleft <=  0.0001
sev=sum(idx_severe)
mod=sum(idx_mod)

cstleftvals=cstleft(idx_severe)
cstrightvals=cstright(idx_mod)


severe_FO1=dwell_avg_stroke(idx_severe,:,1)
mod_FO1=dwell_avg_stroke(idx_mod,:,1)
severe_FO2=dwell_avg_stroke(idx_severe,:,2)
mod_FO2=dwell_avg_stroke(idx_mod,:,2)
severe_FO3=dwell_avg_stroke(idx_severe,:,3)
mod_FO3=dwell_avg_stroke(idx_mod,:,3)
severe_FO4=dwell_avg_stroke(idx_severe,:,4)
mod_FO4=dwell_avg_stroke(idx_mod,:,4)


[rho,p(2)]=corr(cstrightvals,mean(mod_FO3(:,:), 2, 'omitnan'), 'Type', 'Pearson')
[rho,p(1)]=corr(cstleftvals,mean(severe_FO3(:,:), 2, 'omitnan'), 'Type', 'Pearson')
[rho,p]=corr(cst,mean(dwell_avg_stroke(:,:,3), 2, 'omitnan'), 'Type', 'Pearson')


figure
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);

scatter(ax1,cstleftvals,mean(severe_FO3(:,:),2, 'omitnan'), 'ro', 'filled')
xlabel(ax1,'CST lesion load')
ylabel(ax1, 'Dwell time state 3')
title(ax1,'Affected motor system')

scatter(ax2,cstrightvals,mean(mod_FO3(:,:), 2, 'omitnan'), 'bo', 'filled')
xlabel(ax2,'CST lesion load')
ylabel(ax2, 'Dwell time state 3')
title(ax2,'Unaffected motor system')

h1 = lsline(ax1);
h1.Color = 'r';

h1 = lsline(ax2);
h1.Color = 'b';
hold on

severe_FO1=stroke_FO1(idx_severe,:)
mod_FO1=stroke_FO1(idx_mod,:)
severe_FO2=stroke_FO2(idx_severe,:)
mod_FO2=stroke_FO2(idx_mod,:)
severe_FO3=stroke_FO3(idx_severe,:)
mod_FO3=stroke_FO3(idx_mod,:)
severe_FO4=stroke_FO4(idx_severe,:)
mod_FO4=stroke_FO4(idx_mod,:)


severe_FO1=stroke_appearance1(idx_severe,:)
mod_FO1=stroke_appearance1(idx_mod,:)
severe_FO2=stroke_appearance2(idx_severe,:)
mod_FO2=stroke_appearance2(idx_mod,:)
severe_FO3=stroke_appearance3(idx_severe,:)
mod_FO3=stroke_appearance3(idx_mod,:)
severe_FO4=stroke_appearance4(idx_severe,:)
mod_FO4=stroke_appearance4(idx_mod,:)

[rho,p]=corr(cst,dwell_avg_stroke(:,1,2))
[rho,p]=corr(cst,dwell_avg_stroke(:,1,2))

[rho,p]=corr(cstleft,mean(stroke_FO3(:,:),2, 'omitnan'), 'Type', 'Pearson')
[rho,p]=corr(cstleft,mean(dwell_avg_stroke(:,:,4), 2, 'omitnan'), 'Type', 'Pearson')

cstleft(6)=cstright(6)
plot(cstleft,stroke_appearance3(:,1), '*r')
plot(cstleft,mean(stroke_appearance2(:,:), 2, 'omitnan'), '*r')
plot(cstleft,mean(stroke_FO3(:,:),2, 'omitnan'), '*r')

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

[h, p(1), ~, ~]=ttest2(mean(severe_FO1(:,:),2, 'omitnan'), mean(mod_FO1(:,:), 2, 'omitnan'))
[h, p(2), ~, ~]=ttest2(mean(severe_FO2(:,:),2, 'omitnan'), mean(mod_FO2(:,:), 2, 'omitnan'))
[h, p(3), ~, ~]=ttest2(mean(severe_FO3(:,:),2, 'omitnan'), mean(mod_FO3(:,:), 2, 'omitnan'))
[h, p(4), ~, ~]=ttest2(mean(severe_FO4(:,:),2, 'omitnan'), mean(mod_FO4(:,:), 2, 'omitnan'))

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

