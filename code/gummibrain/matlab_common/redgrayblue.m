function C = redgrayblue(m)
%     n1  n2   m/2
%r = 0->0->.25->.5    .5->1->.5
%g = 0->0->.25->.5     .5->0->0
%b = .5->1->.75->.5   .5->0->0


p = [0 .1 .4 .5 .6 .9 1];
c = [0 0 .5; 0 0 1; .25 .25 .75; .5 .5 .5; .75 .25 .25; 1 0 0; .5 0 0];

C = GenerateColormap(p,c,m);