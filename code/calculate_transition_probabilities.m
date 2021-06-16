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