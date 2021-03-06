
figdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/figures/'
resdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/'

%% Severe vs moderate stroke subjects
%lesion load
lesionload=load('/Users/emilyolafson/GIT/dynamic-brainstates/data/lesionload_CSTR_CSTL.mat')
lesionload=lesionload.ll;
cst_r=lesionload(:,1)
cst_l=lesionload(:,2)
r=[1,1,0,0,1,0,1,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,1]

for i=1:23
    if r(i)==1
        cstall(i)=lesionload(i,1);
    else
        cstall(i)=lesionload(i,2);
    end
end

rh=[1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
stroke_dom=~r.*rh

idx_severe=logical(stroke_dom')
idx_mod=~stroke_dom'
sev=sum(idx_severe)
mod=sum(idx_mod)

% correlation between impairment and lesion load only in subjects whose
% dominant hemisphere was impacted
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

for i=1:5
    [h,p(1,i),~,~]=ttest2(severe_FO1(:,i),mod_FO1(:,i))
    [h,p(2,i),~,~]=ttest2(severe_FO2(:,i),mod_FO2(:,i))
    [h,p(3,i),~,~]=ttest2(severe_FO3(:,i),mod_FO3(:,i))
    [h,p(4,i),~,~]=ttest2(severe_FO4(:,i),mod_FO4(:,i))
end

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')
s=1
[rho, p(1)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')
s=2
[rho, p(1)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')
s=3
[rho, p(1)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')
s=4
[rho, p(1)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')
s=5
[rho, p(1)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')



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

for i=1:5
    [h,p(1,i),~,~]=ttest2(severe_FO1(:,i),mod_FO1(:,i))
    [h,p(2,i),~,~]=ttest2(severe_FO2(:,i),mod_FO2(:,i))
    [h,p(3,i),~,~]=ttest2(severe_FO3(:,i),mod_FO3(:,i))
    [h,p(4,i),~,~]=ttest2(severe_FO4(:,i),mod_FO4(:,i))
end
[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')
clear p
s=1
[rho, p(1,s)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2,s)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3,s)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4,s)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')
s=2
[rho, p(1,s)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2,s)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3,s)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4,s)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')
s=3
[rho, p(1,s)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2,s)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3,s)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4,s)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')
s=4
[rho, p(1,s)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2,s)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3,s)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4,s)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')
s=5
[rho, p(1,s)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2,s)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3,s)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4,s)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')

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

for i=1:5
    [h,p(1,i),~,~]=ttest2(severe_FO1(:,i),mod_FO1(:,i))
    [h,p(2,i),~,~]=ttest2(severe_FO2(:,i),mod_FO2(:,i))
    [h,p(3,i),~,~]=ttest2(severe_FO3(:,i),mod_FO3(:,i))
    [h,p(4,i),~,~]=ttest2(severe_FO4(:,i),mod_FO4(:,i))
end
[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

clear p
s=1
[rho, p(1,s)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2,s)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3,s)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4,s)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')
s=2
[rho, p(1,s)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2,s)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3,s)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4,s)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')
s=3
[rho, p(1,s)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2,s)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3,s)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4,s)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')
s=4
[rho, p(1,s)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2,s)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3,s)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4,s)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')
s=5
[rho, p(1,s)]=corr(severe_FO1(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(2,s)]=corr(severe_FO2(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(3,s)]=corr(severe_FO3(:,s), cstall(idx_severe)', 'rows', 'complete')
[rho, p(4,s)]=corr(severe_FO4(:,s), cstall(idx_severe)', 'rows', 'complete')

%% Plot
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
ylabel('Fractional Occupancy')
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
ylabel('Fractional Occupancy')
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

ylabel('Fractional Occupancy')
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
ylabel('Fractional Occupancy')
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
ylabel('Fractional Occupancy')
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



% dwell time
