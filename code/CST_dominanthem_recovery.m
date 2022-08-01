
figdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/shen268/figures/'
resdir='/Users/emilyolafson/GIT/dynamic-brainstates/results/shen268/'
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
fm_1(22)=NaN;
fm_1(23)=NaN


fm41=fm_4-fm_1
fm31=fm_3-fm_1;
fm51=fm_5-fm_1;
fm21=fm_2-fm_1;
fm45=fm_5-fm_4
fmall=[fm_1,fm_2,fm_3,fm_4,fm_5]

sex=[0,1,1,1,1,0,1,1,1,1,1,0,1,1,0,0,1,1,0,0,1,1,0]
age=[54    57    59    48    63    34    60    61    56    68    62    74    55    54    73    62    60    56    64    55    42    40    56]

%% Severe vs moderate stroke subjects
%lesion load
lesionload=load('/Users/emilyolafson/GIT/dynamic-brainstates/data/lesionload_CSTR_CSTL.mat')
lesionload=lesionload.ll;

ll=max(lesionload')'

r=[1,1,0,0,1,0,1,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,1]
idx_lesion=ll>0.01

dom_cst=[3,4, 11, 12, 13,16, 17, 22, 7]
ndom_cst= [1,5,8,9,10,18,19,20,21,23]

tmp=zeros(1,23);
tmp(dom_cst)=1;
idx_dom=logical(tmp);


tmp=zeros(1,23);
tmp(ndom_cst)=1;
idx_ndom=logical(tmp);

dom=sum(idx_dom)
ndom=sum(idx_ndom)

for i=1:23
    if idx_dom(i)==0
        colours(i)='r'
    elseif idx_dom(i)==1
        colours(i)='k'
    end
end


%% Session 4 
ChangeDT=dwell_avg_stroke(:,4,1)-dwell_avg_stroke(:,1,1);
ChangeAR=stroke_appearance1(:,4)-stroke_appearance1(:,1);
ChangeFO=stroke_FO1(:,4)-stroke_FO1(:,1);
ChangeFM=fm41;

Dom_Affected=idx_dom
Ndom_affected=idx_ndom

clear DomAffectedFactor
for i=1:23
    if Dom_Affected(i)==1
        DomAffectedFactor(i)="Dominant CST"
        DomAffectedNumber(i)=1
    else
        DomAffectedFactor(i)="Non-dominant CST"
        DomAffectedNumber(i)=0
    end
end

dataset=[ChangeDT,ChangeAR,ChangeFO, ChangeFM,DomAffectedFactor',DomAffectedNumber', age', sex'];

ds=array2table(dataset)
ds.Properties.VariableNames={'ChangeDT', 'ChangeAR', 'ChangeFO','ChangeFM', 'Dom_Affected_Factor','Dom_Affected_Number', 'age', 'sex'};

writetable(ds, '/Users/emilyolafson/GIT/dynamic-brainstates/data/datatable_session4_20subs.csv')

%regression in R


%% Session 5

ChangeDT=dwell_avg_stroke(:,5,1)-dwell_avg_stroke(:,1,1);
ChangeAR=stroke_appearance1(:,5)-stroke_appearance1(:,1);
ChangeFO=stroke_FO1(:,5)-stroke_FO1(:,1);
ChangeFM=fm51;

Dom_Affected=idx_dom
Ndom_affected=idx_ndom

clear DomAffectedFactor
for i=1:23
    if Dom_Affected(i)==1
        DomAffectedFactor(i)="Dominant CST"
        DomAffectedNumber(i)=1
    else
        DomAffectedFactor(i)="Non-dominant CST"
        DomAffectedNumber(i)=0
    end
end

dataset=[ChangeDT,ChangeAR,ChangeFO, ChangeFM,DomAffectedFactor',DomAffectedNumber', age', sex'];

ds=array2table(dataset)
ds.Properties.VariableNames={'ChangeDT', 'ChangeAR', 'ChangeFO','ChangeFM', 'Dom_Affected_Factor',  'Dom_Affected_Number','age', 'sex'};


writetable(ds, '/Users/emilyolafson/GIT/dynamic-brainstates/data/datatable_session5_20subs.csv')

%regression in R

%% lesion load of CST

lesionload=load('/Users/emilyolafson/GIT/dynamic-brainstates/data/lesionload_CSTR_CSTL.mat')
lesionload=lesionload.ll;

tmp=zeros(1,23);
tmp(dom_cst)=1;

ll=max(lesionload')'


idx_dom=logical(ll>0.0085)
idx_dom(13)=true;
idx_ndom=~logical(idx_dom);

absdom=idx_dom

dom=sum(idx_dom)
ndom=sum(idx_ndom)

for i=1:23
    if idx_dom(i)==1
        AbsDamage(i)="HighDamage"
    else
        AbsDamage(i)="LowDamage"
    end
end





t=readtable('/Users/emilyolafson/GIT/dynamic-brainstates/data/datatable_state2_DT_session4.csv')

t.AbsDamage=AbsDamage'
t.LL=ll
t.DomAffectedFactor2=DomAffectedFactor'
t.age=age'
t.sex=sex'
t([2,14,15],:) = [];


writetable(t, '/Users/emilyolafson/GIT/dynamic-brainstates/data/datatable_state2_DT_session4_withabsD.csv')

t=readtable('/Users/emilyolafson/GIT/dynamic-brainstates/data/datatable_state2_DT_session5.csv')

t.AbsDamage=AbsDamage'
t.LL=ll
t.DomAffectedFactor2=DomAffectedFactor'
t.age=age'
t.sex=sex'
t([2,14,15],:) = [];
size(t)
writetable(t, '/Users/emilyolafson/GIT/dynamic-brainstates/data/datatable_state2_DT_session5_withabsD.csv')


t=readtable('/Users/emilyolafson/GIT/dynamic-brainstates/data/datatable_state4_FO_session4.csv')
t.AbsDamage=AbsDamage'
t.LL=ll
t.DomAffectedFactor2=DomAffectedFactor'
t.age=age'
t.sex=sex'
t([2,14,15],:) = [];

writetable(t, '/Users/emilyolafson/GIT/dynamic-brainstates/data/datatable_state4_FO_session4_withabsD.csv')


t=readtable('/Users/emilyolafson/GIT/dynamic-brainstates/data/datatable_state4_FO_session5.csv')
t.AbsDamage=AbsDamage'
t.LL=ll
t.DomAffectedFactor2=DomAffectedFactor'
t.age=age'
t.sex=sex'
t([2,14,15],:) = [];

writetable(t, '/Users/emilyolafson/GIT/dynamic-brainstates/data/datatable_state4_FO_session5_withabsD.csv')


%

