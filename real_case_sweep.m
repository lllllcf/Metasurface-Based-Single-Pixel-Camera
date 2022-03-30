%orignGraph = readmatrix('SJTU.xlsx');
 orignGraph = imread('SJTU.png');
 orignGraph = orignGraph(:,:,1);
 orignGraph = double(orignGraph);
 orignGraph = orignGraph.*-1+1;
% orignGraph = [orignGraph,zeros(50,10)];
% orignGraph = [zeros(50,10),orignGraph];
% orignGraph = [orignGraph;zeros(5,120)];
% orignGraph = [zeros(5,120);orignGraph];

f = 320e9; %Hz
c = 3e8;
lambda = c/f;
d_m = 580/1e6;
a_on = -9.95;
a_off = -24.96;
phase_on = -206.47;
phase_off = -123.06;
d = 2;

a_on = exp(a_on/20);
a_off = exp(a_off/20);
phase_on = phase_on/180*3.14;
phase_off = phase_off/180*3.14;

[x, y] = size(orignGraph);
Signal = orignGraph(:);
xy = x*y;

finalRes = [];

interval = 100;
for m = 100:interval:2500
    disp(m);
    % Graph: (x,y)  measure times: m
    % Signal = Base * Ksparse   Measure = Pattern * Signal
    % (x*y,1) (x*y, x*y)(x*y,1) (m,1)    (m,x*y)  (x*y,1)
    % Signal = (PB) K  solve b = Ax

    Pattern0 = randi(2,m,xy) - 1;
    Pattern = Pattern0.*(a_on-a_off)+a_off.*ones(size(Pattern0));
    phaseD = Pattern0.*(phase_on-phase_off)+phase_off.*ones(size(Pattern0));
    mD = (repmat(1:x,m,y)-1).*d_m;
    nD = (ceil(repmat(1:xy,m,1)./x)-1).*d_m;
    R = sqrt(d^2+(mD-(d_m*(x-1)/2)).^2+(nD-(d_m*(y-1)/2)).^2);
    Pattern = abs(Pattern.*cos(2.*3.14./lambda.*R+phaseD));
    Measure = Pattern * Signal;
    Base = idct(eye(xy,xy))';

    Ksparse0 = pinv(Pattern*Base)*Measure;
    Ksparse = l1eq_pd(Ksparse0, Pattern*Base, Measure);
    % (x0,A,b)

    Result = reshape(Base * Ksparse, x, y); 
    finalRes = [finalRes, corr2(orignGraph, Result)];
end

plot(100:interval:2500,finalRes);
