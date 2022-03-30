clear, clc;
orignGraph = readmatrix('SJTU.xlsx');
subplot(2,1,1);
contourf(flipud(orignGraph),'LineStyle','none')
colorbar;

m = 800; 

guessMode = false; %仅对 0/1 有效的玩具 需要一半以上的m
guessValue = 0.5;

timeUseless = 4.622e-06*m*m - 0.002499*m + 7.285;
disp(timeUseless);

% Graph: (x,y)  measure times: m
% Signal = Base * Ksparse   Measure = Pattern * Signal
% (x*y,1) (x*y, x*y)(x*y,1) (m,1)    (m,x*y)  (x*y,1)
% Signal = (PB) K  solve b = Ax

[x, y] = size(orignGraph);
Signal = orignGraph(:);
xy = x*y;

d_0 = 0.1;
d_1 = 0.2;
Pattern = (randi(2,m,xy) - 1).*(d_1-d_0)+d_0;
Measure = Pattern * Signal;
Base = idct(eye(xy,xy))';

Ksparse0 = pinv(Pattern*Base)*Measure;
Ksparse = l1eq_pd(Ksparse0, Pattern*Base, Measure);
% (x0,A,b)

Result = reshape(Base * Ksparse, x, y);
if guessMode == true
    Result = (Result >= guessValue);
    
    %孤立点处理
    xpos = 2:x-1;
    ypos = 2:y-1;
    sumMatrix(xpos,ypos) = Result(xpos,ypos-1) + Result(xpos,ypos+1) + Result(xpos-1, ypos) + Result(xpos+1,ypos)...
        + Result(xpos-1,ypos-1) + Result(xpos-1,ypos+1) + Result(xpos+1,ypos-1) + Result(xpos+1,ypos+1);
    singleMatrix = (sumMatrix == 0);
    singleMatrix = logical(resize(singleMatrix, x, y));
    singleMatrix = (Result == 1) & (singleMatrix == 1);
    
    
    Result = int8(Result) - int8(singleMatrix);
end
subplot(2,1,2);
contourf(flipud(Result),'LineStyle','none')
colorbar;

disp(corr2(orignGraph, Result));

function [new] = resize(growthCell, x, y)
    useless1 = zeros(x - 1, 1);
    useless0 = zeros(1, y - 0);
    growthCell = [growthCell useless1];
    new = [growthCell; useless0];
end



