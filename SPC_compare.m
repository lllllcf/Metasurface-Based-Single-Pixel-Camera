clear, clc;
orignGraph = readmatrix('SJTU.xlsx');

interval = 80;
Sim = interval:interval:2500;

for m = interval:interval:2500
    disp(m);
    [x, y] = size(orignGraph);
    Signal = orignGraph(:);
    xy = x*y;

    Pattern = randi(2,m,xy) - 1;
    Measure = Pattern * Signal;
    Base = idct(eye(xy,xy))';

    Ksparse0 = pinv(Pattern*Base)*Measure;
    Ksparse = l1eq_pd(Ksparse0, Pattern*Base, Measure);

    Result = reshape(Base * Ksparse, x, y);
    
    Sim(m/interval) = corr2(orignGraph, Result);
end

plot(interval:interval:2500, Sim);
pp;
