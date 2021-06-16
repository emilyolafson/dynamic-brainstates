%% relationship to impairment and recovery

figdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/figures/'
resdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/'


% load motor scores
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

%% Fractional Occupancy- related to impairment at same scan

[rho, p(1)]=corr(stroke_FO1(:,1), fm_1, 'rows', 'complete')
[rho, p(2)]=corr(stroke_FO2(:,1), fm_1, 'rows', 'complete')
[rho, p(3)]=corr(stroke_FO3(:,1), fm_1, 'rows', 'complete')
[rho, p(4)]=corr(stroke_FO4(:,1), fm_1, 'rows', 'complete')

[rho, p(5)]=corr(stroke_FO1(:,2), fm_2, 'rows', 'complete')
[rho, p(6)]=corr(stroke_FO2(:,2), fm_2, 'rows', 'complete')
[rho, p(7)]=corr(stroke_FO3(:,2), fm_2, 'rows', 'complete')
[rho, p(8)]=corr(stroke_FO4(:,2), fm_2, 'rows', 'complete')

[rho, p(9)]=corr(stroke_FO1(:,3), fm_3, 'rows', 'complete')
[rho, p(10)]=corr(stroke_FO2(:,3), fm_3, 'rows', 'complete')
[rho, p(11)]=corr(stroke_FO3(:,3), fm_3, 'rows', 'complete')
[rho, p(12)]=corr(stroke_FO4(:,3), fm_3, 'rows', 'complete')

[rho, p(13)]=corr(stroke_FO1(:,4), fm_4, 'rows', 'complete')
[rho, p(14)]=corr(stroke_FO2(:,4), fm_4, 'rows', 'complete')
[rho, p(15)]=corr(stroke_FO3(:,4), fm_4, 'rows', 'complete')
[rho, p(16)]=corr(stroke_FO4(:,4), fm_4, 'rows', 'complete')

[rho, p(17)]=corr(stroke_FO1(:,5), fm_5, 'rows', 'complete')
[rho, p(18)]=corr(stroke_FO2(:,5), fm_5, 'rows', 'complete')
[rho, p(19)]=corr(stroke_FO3(:,5), fm_5, 'rows', 'complete')
[rho, p(20)]=corr(stroke_FO4(:,5), fm_5, 'rows', 'complete')

% diff in parameters related to change in fm score?
fm12=fm_2-fm_1;
fm23=fm_3-fm_2
fm34=fm_4-fm_3;
fm45=fm_5-fm_4;

fo_12=stroke_FO1(:,2)-stroke_FO1(:,1);
fo_23=stroke_FO1(:,3)-stroke_FO1(:,2);
fo_34=stroke_FO1(:,4)-stroke_FO1(:,3);
fo_45=stroke_FO1(:,5)-stroke_FO1(:,4);

clear p
[rho,p(1)]=corr(fm12, fo_12, 'rows', 'complete')
[rho,p(2)]=corr(fm23, fo_23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, fo_34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, fo_45, 'rows', 'complete')

%%  Dwell time - related to impairment at same scan
[rho, p(1)]=corr(dwell_avg_stroke(:,1,1), fm_1, 'rows', 'complete')
[rho, p(2)]=corr(dwell_avg_stroke(:,1,2), fm_1, 'rows', 'complete')
[rho, p(3)]=corr(dwell_avg_stroke(:,1,3), fm_1, 'rows', 'complete')
[rho, p(4)]=corr(dwell_avg_stroke(:,1,4), fm_1, 'rows', 'complete')

[rho, p(5)]=corr(dwell_avg_stroke(:,2,1), fm_2, 'rows', 'complete')
[rho, p(6)]=corr(dwell_avg_stroke(:,2,2), fm_2, 'rows', 'complete')
[rho, p(7)]=corr(dwell_avg_stroke(:,2,3), fm_2, 'rows', 'complete')
[rho, p(8)]=corr(dwell_avg_stroke(:,2,4), fm_2, 'rows', 'complete')

[rho, p(9)]=corr(dwell_avg_stroke(:,3,1), fm_3, 'rows', 'complete')
[rho, p(10)]=corr(dwell_avg_stroke(:,3,2), fm_3, 'rows', 'complete')
[rho, p(11)]=corr(dwell_avg_stroke(:,3,3), fm_3, 'rows', 'complete')
[rho, p(12)]=corr(dwell_avg_stroke(:,3,4), fm_3, 'rows', 'complete')

[rho, p(13)]=corr(dwell_avg_stroke(:,4,1), fm_4, 'rows', 'complete')
[rho, p(14)]=corr(dwell_avg_stroke(:,4,2), fm_4, 'rows', 'complete')
[rho, p(15)]=corr(dwell_avg_stroke(:,4,3), fm_4, 'rows', 'complete')
[rho, p(16)]=corr(dwell_avg_stroke(:,4,4), fm_4, 'rows', 'complete')

[rho, p(17)]=corr(dwell_avg_stroke(:,5,1), fm_5, 'rows', 'complete')
[rho, p(18)]=corr(dwell_avg_stroke(:,5,2), fm_5, 'rows', 'complete')
[rho, p(19)]=corr(dwell_avg_stroke(:,5,3), fm_5, 'rows', 'complete')
[rho, p(20)]=corr(dwell_avg_stroke(:,5,4), fm_5, 'rows', 'complete')


%state 1
dt12=dwell_avg_stroke(:,2,1)-dwell_avg_stroke(:,1,1)
dt23=dwell_avg_stroke(:,3,1)-dwell_avg_stroke(:,2,1)
dt34=dwell_avg_stroke(:,4,1)-dwell_avg_stroke(:,3,1)
dt45=dwell_avg_stroke(:,5,1)-dwell_avg_stroke(:,4,1)

clear p
[rho,p(1)]=corr(fm12, dt12, 'rows', 'complete')
[rho,p(2)]=corr(fm23, dt23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, dt34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, dt45, 'rows', 'complete')

% state 2
dt12=dwell_avg_stroke(:,2,2)-dwell_avg_stroke(:,1,2)
dt23=dwell_avg_stroke(:,3,2)-dwell_avg_stroke(:,2,2)
dt34=dwell_avg_stroke(:,4,2)-dwell_avg_stroke(:,3,2)
dt45=dwell_avg_stroke(:,5,2)-dwell_avg_stroke(:,4,2)

clear p
[rho,p(1)]=corr(fm12, dt12, 'rows', 'complete')
[rho,p(2)]=corr(fm23, dt23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, dt34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, dt45, 'rows', 'complete')

% state 3
dt12=dwell_avg_stroke(:,2,3)-dwell_avg_stroke(:,1,3)
dt23=dwell_avg_stroke(:,3,3)-dwell_avg_stroke(:,2,3)
dt34=dwell_avg_stroke(:,4,3)-dwell_avg_stroke(:,3,3)
dt45=dwell_avg_stroke(:,5,3)-dwell_avg_stroke(:,4,3)

clear p
[rho,p(1)]=corr(fm12, dt12, 'rows', 'complete', 'Type', 'Spearman')
[rho,p(2)]=corr(fm23, dt23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, dt34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, dt45, 'rows', 'complete')

close all
figure()
scatter(fm12, dt12, 'ok', 'filled')
xlabel('Change in Fugl-Meyer scores 1-2 weeks post-stroke')
ylabel('Change in State 3 dwell time 1-2 weeks post-stroke')
text(30, -0.4, 'Corr = 0.53, p = 0.013 (unc)', 'FontSize', 14)
set(gca, 'FontSize', 14)
lsline
title('State 3 change in dwell time vs. change in FM')
[h, ~, ~, p_adj] = fdr_bh(p, 0.05,'pdep')


% state 4
dt12=dwell_avg_stroke(:,2,4)-dwell_avg_stroke(:,1,4)
dt23=dwell_avg_stroke(:,3,4)-dwell_avg_stroke(:,2,4)
dt34=dwell_avg_stroke(:,4,4)-dwell_avg_stroke(:,3,4)
dt45=dwell_avg_stroke(:,5,4)-dwell_avg_stroke(:,4,4)

clear p
[rho,p(1)]=corr(fm12, dt12, 'rows', 'complete')
[rho,p(2)]=corr(fm23, dt23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, dt34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, dt45, 'rows', 'complete')

%% Appearance rate
clear p
[rho, p(1)]=corr(stroke_appearance1(:,1), fm_1, 'rows', 'complete')
[rho, p(2)]=corr(stroke_appearance2(:,1), fm_1, 'rows', 'complete')
[rho, p(3)]=corr(stroke_appearance3(:,1), fm_1, 'rows', 'complete')
[rho, p(4)]=corr(stroke_appearance4(:,1), fm_1, 'rows', 'complete')

[rho, p(5)]=corr(stroke_appearance1(:,2), fm_2, 'rows', 'complete')
[rho, p(6)]=corr(stroke_appearance2(:,2), fm_2, 'rows', 'complete')
[rho, p(7)]=corr(stroke_appearance3(:,2), fm_2, 'rows', 'complete')
[rho, p(8)]=corr(stroke_appearance4(:,2), fm_2, 'rows', 'complete')

[rho, p(9)]=corr(stroke_appearance1(:,3), fm_3, 'rows', 'complete')
[rho, p(10)]=corr(stroke_appearance2(:,3), fm_3, 'rows', 'complete')
[rho, p(11)]=corr(stroke_appearance3(:,3), fm_3, 'rows', 'complete')
[rho, p(12)]=corr(stroke_appearance4(:,3), fm_3, 'rows', 'complete')

[rho, p(13)]=corr(stroke_appearance1(:,4), fm_4, 'rows', 'complete')
[rho, p(14)]=corr(stroke_appearance2(:,4), fm_4, 'rows', 'complete')
[rho, p(15)]=corr(stroke_appearance3(:,4), fm_4, 'rows', 'complete')
[rho, p(16)]=corr(stroke_appearance4(:,4), fm_4, 'rows', 'complete')

[rho, p(17)]=corr(stroke_appearance1(:,5), fm_5, 'rows', 'complete')
[rho, p(18)]=corr(stroke_appearance2(:,5), fm_5, 'rows', 'complete')
[rho, p(19)]=corr(stroke_appearance3(:,5), fm_5, 'rows', 'complete')
[rho, p(20)]=corr(stroke_appearance4(:,5), fm_5, 'rows', 'complete')

% diff in parameters related to change in fm score?
fm12=fm_2-fm_1;
fm23=fm_3-fm_2
fm34=fm_4-fm_3;
fm45=fm_5-fm_4;

clear p
fo_12=stroke_appearance1(:,2)-stroke_appearance1(:,1);
fo_23=stroke_appearance1(:,3)-stroke_appearance1(:,2);
fo_34=stroke_appearance1(:,4)-stroke_appearance1(:,3);
fo_45=stroke_appearance1(:,5)-stroke_appearance1(:,4);

[rho,p(1)]=corr(fm12, fo_12, 'rows', 'complete')
[rho,p(2)]=corr(fm23, fo_23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, fo_34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, fo_45, 'rows', 'complete')

clear p
fo_12=stroke_appearance2(:,2)-stroke_appearance2(:,1);
fo_23=stroke_appearance2(:,3)-stroke_appearance2(:,2);
fo_34=stroke_appearance2(:,4)-stroke_appearance2(:,3);
fo_45=stroke_appearance2(:,5)-stroke_appearance2(:,4);

[rho,p(1)]=corr(fm12, fo_12, 'rows', 'complete')
[rho,p(2)]=corr(fm23, fo_23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, fo_34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, fo_45, 'rows', 'complete')

clear p
fo_12=stroke_appearance3(:,2)-stroke_appearance3(:,1);
fo_23=stroke_appearance3(:,3)-stroke_appearance3(:,2);
fo_34=stroke_appearance3(:,4)-stroke_appearance3(:,3);
fo_45=stroke_appearance3(:,5)-stroke_appearance3(:,4);

[rho,p(1)]=corr(fm12, fo_12, 'rows', 'complete')
[rho,p(2)]=corr(fm23, fo_23, 'rows', 'complete')
[rho,p(3)]=corr(fm34, fo_34, 'rows', 'complete')
[rho,p(4)]=corr(fm45, fo_45, 'rows', 'complete')


