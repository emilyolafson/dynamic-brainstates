
figdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/figures/'
resdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/'
load('partitions_k4_50reps.mat')


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

dom


rh=[1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
stroke_dom=~r.*rh %people whose dominant hand was lesioned

idx_severe=logical(stroke_dom') %severe=cst of rh affected
idx_mod=~stroke_dom'
sev=sum(idx_severe)
mod=sum(idx_mod)

[h, p, ~, stats]=ttest2(cstall(idx_severe),cstall(idx_mod))

histogram(cstall(idx_severe), sev);
hold on;
histogram(cstall(idx_mod),mod)

[rho,p]=corr(cstall(idx_severe)', diffFO2_sev, 'rows','complete')

idx_severe=cstall>0.009
idx_mod=cstall<=0.009
sev=sum(idx_severe)
mod=sum(idx_mod)

%% Dwell time - moderate vs severe
% across single sessions
fm41=fm_4-fm_1
fm31=fm_3-fm_1;
fm51=fm_5-fm_1;
fm21=fm_2-fm_1;

s=5
severe_FO2=dwell_avg_stroke(idx_severe,:,2)
mod_FO2=dwell_avg_stroke(idx_mod,:,2)
diffDT2_sev=severe_FO2(:,s)-severe_FO2(:,1);
diffDT2_mod=mod_FO2(:,s)-mod_FO2(:,1);

severe_FO2=stroke_FO2(idx_severe,:)
mod_FO2=stroke_FO2(idx_mod,:)
diffFO2_sev=severe_FO2(:,s)-severe_FO2(:,1)
diffFO2_mod=mod_FO2(:,s)-mod_FO2(:,1);

severe_FO2=stroke_appearance2(idx_severe,:)
mod_FO2=stroke_appearance2(idx_mod,:)
diffAR2_sev=severe_FO2(:,s)-severe_FO2(:,1);
diffAR2_mod=mod_FO2(:,s)-mod_FO2(:,1);

clear p

[rhos(1), p(1)]=corr(diffFO2_sev, fm51(idx_severe), 'rows', 'complete')
[rhos(2), p(2)]=corr(diffDT2_sev, fm51(idx_severe), 'rows', 'complete')
[rhos(3), p(3)]=corr(diffAR2_sev, fm51(idx_severe), 'rows', 'complete')
[rhos(4), p(4)]=corr(diffFO2_mod, fm51(idx_mod), 'rows', 'complete')
[rhos(5), p(5)]=corr(diffDT2_mod, fm51(idx_mod), 'rows', 'complete')
[rhos(6), p(6)]=corr(diffAR2_mod, fm51(idx_mod), 'rows', 'complete')

[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')


c=1
figure('Position', [0 0 1600 1000])
tiledlayout(2, 3, 'padding', 'none')
nexttile;
scatter(diffFO2_sev,fm41(idx_severe), 'ko', 'filled')
title('Fractional Occupancy')
xlabel('Change in FO (3 month - 1 week)')
ylabel('Change in Fugl-Meyer score (3 month - 1 week)')
message = sprintf('Correlation = %.4f \n p(FDR) = %0.4f', rhos(c), p_adj(c))
text(-0.10, 60, message, 'FontSize', 15)
%ylim([0 80])
set(gca, 'FontSize', 15)
lsline

nexttile;c=c+1;
scatter(diffDT2_sev,fm41(idx_severe), 'ko', 'filled')
title('Dwell Time')
xlabel('Change in DT (3 month - 1 week)')
ylabel('Change in Fugl-Meyer score (3 month - 1 week)')
message = sprintf('Correlation = %.4f \n p(FDR) = %0.4f', rhos(c), p_adj(c))
text(-0.50, 60, message, 'FontSize', 15)
%ylim([0 80])
set(gca, 'FontSize', 15)

lsline
nexttile;c=c+1;
scatter(diffAR2_sev,fm41(idx_severe), 'ko', 'filled')
title('Appearance Rate')
xlabel('Change in AR (3 month - 1 week)')
ylabel('Change in Fugl-Meyer score (3 month - 1 week)')
message = sprintf('Correlation = %.4f \n p(FDR) = %0.4f', rhos(c), p_adj(c))
text(-0.60, 60, message, 'FontSize', 15)
lsline
%ylim([0 80])
set(gca, 'FontSize', 15, 'FontSize', 15)

nexttile;c=c+1;
scatter(diffFO2_mod,fm41(idx_mod), 'ko', 'filled')
title('Fractional Occupancy')
xlabel('Change in FO (3 month - 1 week)')
ylabel('Change in Fugl-Meyer score (3 month - 1 week)')
message = sprintf('Correlation = %.4f \n p(FDR) = %0.4f', rhos(c), p_adj(c))
text(-0.05, 60, message, 'FontSize', 15)
lsline
%ylim([0 80])
set(gca, 'FontSize', 15)

nexttile;c=c+1;
scatter(diffDT2_mod,fm41(idx_mod), 'ko', 'filled')
title('Dwell Time')
xlabel('Change in DT (3 month - 1 week)')
ylabel('Change in Fugl-Meyer score (3 month - 1 week)')
message = sprintf('Correlation = %.4f \n p(FDR) = %0.4f', rhos(c), p_adj(c))
text(-0.140, 60, message, 'FontSize', 15)
lsline
%ylim([0 80])
set(gca, 'FontSize', 15)

nexttile;c=c+1;
scatter(diffAR2_mod,fm41(idx_mod), 'ko', 'filled')
title('Appearance Rate')
xlabel('Change in AR (3 month - 1 week)')
ylabel('Change in Fugl-Meyer score (3 month - 1 week)')
message = sprintf('Correlation = %.4f \n p(FDR) = %0.4f', rhos(c), p_adj(c))
text(-0.40, 65, message, 'FontSize', 15)
lsline
%ylim([0 80])
set(gca, 'FontSize', 15)


saveas(gcf,strcat(figdir, 'changeFM_3month_1week.png'))


severe_FO2=stroke_FO2(idx_severe,1)
severe_FO2=dwell_avg_stroke(idx_severe,:,2)

[rho, p]=corr(mean(severe_FO2, 2, 'omitnan'), cstall')



