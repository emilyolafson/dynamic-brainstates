
figdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/5states/figures/'
resdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/5states.'
load('partitions_k4_50reps.mat')

fm_dir=strcat('/Users/emilyolafson/GIT/stroke-graph-matching/data/');
fuglmeyer=readtable(strcat(fm_dir, 'fuglmeyer_allpts.csv'));
fm_1=fuglmeyer.Var2;
fm_2=fuglmeyer.Var3;
fm_3=fuglmeyer.Var4;
fm_4=fuglmeyer.Var5;
fm_5=fuglmeyer.Var6;

fm_1(23)=NaN;
fm_1(22)=NaN;
fm_3(20)=NaN;
fm_4(12)=NaN;
fm_4(20)=NaN;
fm_5(20)=NaN;
fm_5(6)=NaN;
%% Severe vs moderate stroke subjects
%lesion load
lesionload=load('/Users/emilyolafson/GIT/dynamic-brainstates/data/lesionload_CSTR_CSTL.mat')
lesionload=lesionload.ll;

r=[1,1,0,0,1,0,1,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,1]

dom_cst=[3,4, 6, 11, 12, 13,16, 17, 22, 7]

tmp=zeros(1,23);
tmp(dom_cst)=1;

idx_dom=logical(tmp);
idx_ndom=~logical(tmp);

dom=sum(idx_dom)
ndom=sum(idx_ndom)

%% Dwell time - moderate vs severe
% across single sessions
fm_1(22)=0
fm_1(23)=0
fm41=fm_4-fm_1
fm31=fm_3-fm_1;
fm51=fm_5-fm_1;
fm21=fm_2-fm_1;
fm45=fm_5-fm_4

ChangeDT(idx_dom)
ChangeDT=dwell_avg_stroke(:,5,5)-dwell_avg_stroke(:,1,5);
ChangeFM=fm51;
Dom_Affected=idx_dom;

for i=1:23
    if Dom_Affected(i)==1
        DomAffectedFactor(i)="Dominant CST"
    else
        DomAffectedFactor(i)="Non-dominant CST"
    end
end

dataset=[ChangeDT,ChangeFM,DomAffectedFactor', fm_1, fm_2, fm_3, fm_4, fm_5];

ds=array2table(dataset)
ds.Properties.VariableNames={'ChangeDT', 'ChangeFM', 'Dom_Affected_Factor', 'FM1', 'FM2', 'FM3', 'FM4', 'FM5'};

writetable(ds, '/Users/emilyolafson/GIT/dynamic-brainstates/data/datatable_k5.csv')

%regression in R

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')

c=1
figure('Position', [0 0 500 1000])
tiledlayout(2, 1, 'padding', 'none')
nexttile;
scatter(diffDT2_sev,fm41(idx_dom), 'ko', 'filled')
title('Dwell Time')
xlabel('Change in DT (3 month - 1 week)')
ylabel('Change in Fugl-Meyer score (3 month - 1 week)')
message = sprintf('Correlation = %.4f \n p(FDR) = %0.4f', rhos(c), p_adj(c))
text(-0.50, 60, message, 'FontSize', 15)
%ylim([0 80])
set(gca, 'FontSize', 15)
lsline

nexttile;c=c+1;
scatter(diffDT2_mod,fm41(idx_ndom), 'ko', 'filled')
title('Dwell Time')
xlabel('Change in DT (3 month - 1 week)')
ylabel('Change in Fugl-Meyer score (3 month - 1 week)')
message = sprintf('Correlation = %.4f \n p(FDR) = %0.4f', rhos(c), p_adj(c))
text(-0.65, 60, message, 'FontSize', 15)
lsline
%ylim([0 80])
set(gca, 'FontSize', 15)


severe_FO2=stroke_FO2(idx_severe,1)
severe_FO2=dwell_avg_stroke(idx_severe,:,2)

[rho, p]=corr(mean(severe_FO2, 2, 'omitnan'), cstall')



