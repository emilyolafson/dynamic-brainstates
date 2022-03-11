%

%% Analysis - difference between stroke & control transition probabilities
% perm test
nperms = 100000;
%pvals_twotail = PERM_TEST(tran_prob_stroke,tran_prob_control,nperms);
    % INPUT:
    % A and B: NxKxK stack of KxK transition matrices for N subjects and two conditions, A and B
    % nperms: number of permutations, pretty fast so 10000 is good to get a precise p-value

%% Transition probabilities: separate sessions
% across only session 1
allstroke=tran_prob_stroke(:,1)
allstroke_trans=[]
for i=1:size(allstroke,1)
    tmp=cell2mat(allstroke(i));
    if isempty(tmp)
        continue
    end
     allstroke_trans=cat(3, allstroke_trans,tmp)
end


allstroke_trans=permute(allstroke_trans, [3 1 2])
%format control
allctl=tran_prob_control(:,1)
tmp=cell2mat(allctl(i));
allcontrol_trans=[]
for i=1:size(allctl,1)
    tmp=cell2mat(allctl(i));
    if isempty(tmp)
        continue
    end
     allcontrol_trans=cat(3, allcontrol_trans,tmp);
end
allcontrol_trans=permute(allcontrol_trans, [3 1 2])
pvals_twotail_S1 = PERM_TEST(allstroke_trans,allcontrol_trans,nperms);
% only across session 2
allstroke=tran_prob_stroke(:,2)
allstroke_trans=[]
for i=1:size(allstroke,1)
    tmp=cell2mat(allstroke(i));
    if isempty(tmp)
        continue
    end
     allstroke_trans=cat(3, allstroke_trans,tmp)
end
allstroke_trans=permute(allstroke_trans, [3 1 2])
%format control
allctl=tran_prob_control(:,2)
tmp=cell2mat(allctl(i));
allcontrol_trans=[]
for i=1:size(allctl,1)
    tmp=cell2mat(allctl(i));
    if isempty(tmp)
        continue
    end
     allcontrol_trans=cat(3, allcontrol_trans,tmp);
end
allcontrol_trans=permute(allcontrol_trans, [3 1 2])
pvals_twotail_S2 = PERM_TEST(allstroke_trans,allcontrol_trans,nperms);
% session 3
allstroke=tran_prob_stroke(:,3)
allstroke_trans=[]
for i=1:size(allstroke,1)
    tmp=cell2mat(allstroke(i));
    if isempty(tmp)
        continue
    end
     allstroke_trans=cat(3, allstroke_trans,tmp)
end
allstroke_trans=permute(allstroke_trans, [3 1 2])
%format control
allctl=tran_prob_control(:,3)
tmp=cell2mat(allctl(i));
allcontrol_trans=[]
for i=1:size(allctl,1)
    tmp=cell2mat(allctl(i));
    if isempty(tmp)
        continue
    end
     allcontrol_trans=cat(3, allcontrol_trans,tmp);
end
allcontrol_trans=permute(allcontrol_trans, [3 1 2])
pvals_twotail_S3 = PERM_TEST(allstroke_trans,allcontrol_trans,nperms);
% session 4
allstroke=tran_prob_stroke(:,4)
allstroke_trans=[]
for i=1:size(allstroke,1)
    tmp=cell2mat(allstroke(i));
    if isempty(tmp)
        continue
    end
     allstroke_trans=cat(3, allstroke_trans,tmp)
end
allstroke_trans=permute(allstroke_trans, [3 1 2])
%format control
allctl=tran_prob_control(:,4)
tmp=cell2mat(allctl(i));
allcontrol_trans=[]
for i=1:size(allctl,1)
    tmp=cell2mat(allctl(i));
    if isempty(tmp)
        continue
    end
     allcontrol_trans=cat(3, allcontrol_trans,tmp);
end
allcontrol_trans=permute(allcontrol_trans, [3 1 2])
pvals_twotail_S4 = PERM_TEST(allstroke_trans,allcontrol_trans,nperms);


% session 5
allstroke=tran_prob_stroke(:,5)
allstroke_trans=[]
for i=1:size(allstroke,1)
    tmp=cell2mat(allstroke(i));
    if isempty(tmp)
        continue
    end
     allstroke_trans=cat(3, allstroke_trans,tmp)
end
allstroke_trans=permute(allstroke_trans, [3 1 2])
%format control
allctl=tran_prob_control(:,5)
tmp=cell2mat(allctl(i));
allcontrol_trans=[]
for i=1:size(allctl,1)
    tmp=cell2mat(allctl(i));
    if isempty(tmp)
        continue
    end
     allcontrol_trans=cat(3, allcontrol_trans,tmp);
end
allcontrol_trans=permute(allcontrol_trans, [3 1 2])
pvals_twotail_S5 = PERM_TEST(allstroke_trans,allcontrol_trans,nperms);

pvals_twotail_all=[]
pvals_twotail_all{1}=pvals_twotail_S1;
pvals_twotail_all{2}=pvals_twotail_S2;
pvals_twotail_all{3}=pvals_twotail_S3;
pvals_twotail_all{4}=pvals_twotail_S4;
pvals_twotail_all{5}=pvals_twotail_S5;
[h, ~, ~, p_adjs1] = fdr_bh(pvals_twotail_S1, 0.05,'pdep')
[h, ~, ~, p_adjs2] = fdr_bh(pvals_twotail_S2, 0.05,'pdep')
[h, ~, ~, p_adjs3] = fdr_bh(pvals_twotail_S3, 0.05,'pdep')
[h, ~, ~, p_adjs4] = fdr_bh(pvals_twotail_S4, 0.05,'pdep')
[h, ~, ~, p_adjs5] = fdr_bh(pvals_twotail_S5, 0.05,'pdep')

pvals_twotail_corr=[]
pvals_twotail_corr{1}=p_adjs1;
pvals_twotail_corr{2}=p_adjs2;
pvals_twotail_corr{3}=p_adjs3;
pvals_twotail_corr{4}=p_adjs4;
pvals_twotail_corr{5}=p_adjs5;

%%  get t-stats for transition probs for separate sessions

%test significance between transition probabilities for each session
%separately
s1stroke_trans=[]
s1control_trans=[]
% session 1
for a=1:4
    for b=1:4
        s1stroke=tran_prob_stroke(:,1)
        for i=1:size(s1stroke,1)
            tmp=cell2mat(s1stroke(i));
            s1stroke_trans(i)=tmp(a,b);
        end
        s1control=tran_prob_control(:,1)
        for i=1:size(s1control,1)
            tmp=cell2mat(s1control(i));
            s1control_trans(i)=tmp(a,b);
        end
        [h, p1(a,b),~, stat]=ttest2(s1stroke_trans, s1control_trans)
        stats1(a,b)=stat.tstat;

    end
end


% session 2
s1stroke_trans=[]
s1control_trans=[]
clear p2
for a=1:4
    for b=1:4
        s1stroke=tran_prob_stroke(:,2)
        for i=1:size(s1stroke,1)
            tmp=cell2mat(s1stroke(i));
            s1stroke_trans(i)=tmp(a,b);
        end
         s1control=tran_prob_control(:,2)
        for i=1:size(s1control,1)
            tmp=cell2mat(s1control(i));
            s1control_trans(i)=tmp(a,b);
        end
        [h, p2(a,b),~, stat]=ttest2(s1stroke_trans, s1control_trans)
         stats2(a,b)=stat.tstat;

    end
end

% session 3
s1stroke_trans=[]
s1control_trans=[]
clear p3
for a=1:4
    for b=1:4
        s1stroke=tran_prob_stroke(:,3)
        for i=1:size(s1stroke,1)
            if i==20
                continue
            end
            tmp=cell2mat(s1stroke(i));
            s1stroke_trans(i)=tmp(a,b);
        end
         s1control=tran_prob_control(:,3)
        for i=1:size(s1control,1)
            tmp=cell2mat(s1control(i));
            s1control_trans(i)=tmp(a,b);
        end
        [h, p3(a,b),~, stat]=ttest2(s1stroke_trans, s1control_trans)
         stats3(a,b)=stat.tstat;

    end
end

% session 4
s1stroke_trans=[]
s1control_trans=[]
clear p4
for a=1:4
    for b=1:4
        s1stroke=tran_prob_stroke(:,4)
        for i=1:size(s1stroke,1)
            if i==20
                continue;
            elseif i==12
                continue;
            end
            tmp=cell2mat(s1stroke(i));
            s1stroke_trans(i)=tmp(a,b);
        end
         s1control=tran_prob_control(:,4)
        for i=1:size(s1control,1)
            tmp=cell2mat(s1control(i));
            s1control_trans(i)=tmp(a,b);
        end
        [h, p4(a,b),~, stat]=ttest2(s1stroke_trans, s1control_trans)
         stats4(a,b)=stat.tstat;

    end
end

% session 5
s1stroke_trans=[]
s1control_trans=[]
clear p5
for a=1:4
    for b=1:4
        s1stroke=tran_prob_stroke(:,5)
        for i=1:size(s1stroke,1)
             if i==20
                continue;
            elseif i==12
                continue;
             elseif i==6
                 continue;
            end
            tmp=cell2mat(s1stroke(i));
            s1stroke_trans(i)=tmp(a,b);
        end
        s1control=tran_prob_control(:,5)
        for i=1:size(s1control,1)
            tmp=cell2mat(s1control(i));
            s1control_trans(i)=tmp(a,b);
        end
        [h, p5(a,b),~, stat]=ttest2(s1stroke_trans, s1control_trans)
        
        stats5(a,b)=stat.tstat;
    end
end

[h, ~, ~, p1_adj] = fdr_bh(p1, 0.05,'pdep')
[h, ~, ~, p2_adj] = fdr_bh(p2, 0.05,'pdep')
[h, ~, ~, p3_adj] = fdr_bh(p3, 0.05,'pdep')
[h, ~, ~, p4_adj] = fdr_bh(p4, 0.05,'pdep')
[h, ~, ~, p5_adj] = fdr_bh(p5, 0.05,'pdep')

alltstats{1}=stats1;
alltstats{2}=stats2;
alltstats{3}=stats3;
alltstats{4}=stats4;
alltstats{5}=stats5;

allpvals{1}=p1;
allpvals{2}=p2;
allpvals{3}=p3;
allpvals{4}=p4;
allpvals{5}=p5;

allpvals_adj{1}=p1_adj;
allpvals_adj{2}=p2_adj;
allpvals_adj{3}=p3_adj;
allpvals_adj{4}=p4_adj;
allpvals_adj{5}=p5_adj;





%% Transition probabilities: averaged across sessions
stest2=[]
for i=1:23
    transp=[]
    allstroke=tran_prob_stroke(i,:)
    if i==20 || i==12 || i==6
        continue
    end
    stest2(i)=allstroke{5}(4,3);
    for j=1:size(allstroke,2)
        transp=cat(3,transp,allstroke{j})
    end
    average_transp_str{i}=mean(transp,3)
end

for i=1:24
    transp=[]
    allcontrol=tran_prob_control(i,:)
    ctest2(i)=allcontrol{5}(4,3);
    for j=1:size(allcontrol,2)
        transp=cat(3,transp,allcontrol{j})
    end
    average_transp_ctl{i}=mean(transp,3)
end

figure()
histogram(stest);hold on; histogram(ctest)
figure()

histogram(stest2);hold on; histogram(ctest2)
title('s5')
% t-tests to find stroke-control differences
% format average_transp_str to 4x4xN
average_transp_str_form=[]
for j=1:size(average_transp_str,2)
    average_transp_str_form=cat(3,average_transp_str_form,average_transp_str{j})
end
average_transp_ctl_form=[]
for j=1:size(average_transp_ctl,2)
    average_transp_ctl_form=cat(3,average_transp_ctl_form,average_transp_ctl{j})
end

imagesc(squeeze(average_transp_ctl_form(:,:,1)))

p=ones(4,4)
%w persist
for a=1:4
    for b=1:4
        if isnan(squeeze(average_transp_str_form(a,b,:)))
            disp('isnan')
        elseif isnan(squeeze(average_transp_ctl_form(a,b,:)))
            disp('isnan2')
        end
        [h, p(a,b),~, stat]=ttest2(squeeze(average_transp_str_form(a,b,:)), squeeze(average_transp_ctl_form(a,b,:)))
        stats_all(a,b)=stat.tstat
    end
end

p_all=p;
[h, ~, ~, p_adj_all] = fdr_bh(p, 0.05,'pdep')




%% correlation between transition probabilities and fugl-meyer scores?

% Load Fugl-Meyer scores
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

fm_1(23)=0
fm_1(22)=0




%linear regression fugl meyer vs trans prob

% get predictors
counter=1;
transprobs_vec=[]
for a=1:4
    for b=1:4
        s1stroke=tran_prob_stroke(:,1)
        for i=1:size(s1stroke,1)
            tmp=cell2mat(s1stroke(i));
            s1stroke_trans(i)=tmp(a,b);
        end
        transprobs_vec(:,counter)=s1stroke_trans';
        counter=counter+1
    end
end

fitlm(transprobs_vec,fm_1)



%separately
s1stroke_trans=[]
s1control_trans=[]
% session 1
for a=1:4
    for b=1:4
        s1stroke=tran_prob_stroke(:,1)
        s5stroke=tran_prob_stroke_nopersist(:,5)
        for i=1:size(s1stroke,1)
            tmp=cell2mat(s1stroke(i));
            s1stroke_trans(i)=tmp(a,b);
        end
        one=s1stroke_trans
        [rho1(a,b), p1(a,b)]=corr(one', fm_1)
    end
end


rho2=[]
p2=[]
for a=1:4
    for b=1:4
        s1stroke=tran_prob_stroke(:,1)
        s5stroke=tran_prob_stroke(:,5)
        for i=1:size(s1stroke,1)
            tmp1=cell2mat(s1stroke(i));  
            s1stroke_trans(i)=tmp1(a,b);
        end
        for i=1:size(s1stroke,1)
             if i==20
                continue;
            elseif i==12
                continue;
             elseif i==6
                 continue;
            end
            tmp2=cell2mat(s5stroke(i));  
            s5stroke_trans(i)=tmp2(a,b);
        end
        one=s1stroke_trans(idx_dom)
        two=s5stroke_trans(idx_dom)

        [rho2(a,b), p2(a,b)]=corr((two-one)', fm_5(idx_dom)-fm_1(idx_dom), 'rows', 'complete')
    end
end
imagesc(rho2)




















