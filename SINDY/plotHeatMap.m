
%% plot heatmaps


%% average data
nWrongTermsSm       = mean(nWrongTermsS,3)'; % SINDy
modelErrorSm        = mean(modelErrorS,3)'; % SINDy
successSm           = mean(successS,3)'; % SINDy

nWrongTermsEmedm    = mean(nWrongTermsE,3)'; % Bragging SINDy
nWrongTermsEm       = mean(nWrongTermsE2,3)'; % Bagging SINDy
modelErrorEmedm     = mean(modelErrorE,3)'; % Bragging SINDy
modelErrorEm        = mean(modelErrorE2,3)'; % Bagging SINDy
successEmedm        = mean(successE,3)'; % Bragging SINDy
successEm           = mean(successE2,3)'; % Bagging SINDy

% without replacement
nWrongTermsWREmedm    = mean(nWrongTermsWRE,3)'; % Bragging SINDy
nWrongTermsWREm       = mean(nWrongTermsWRE2,3)'; % Bagging SINDy
modelErrorWREmedm     = mean(modelErrorWRE,3)'; % Bragging SINDy
modelErrorWREm        = mean(modelErrorWRE2,3)'; % Bagging SINDy
successWREmedm        = mean(successWRE,3)'; % Bragging SINDy
successWREm           = mean(successWRE2,3)'; % Bagging SINDy

nWrongTermsDEm      = mean(nWrongTermsDE,3)'; % Library Bagging 
nWrongTermsDDEm     = mean(nWrongTermsDDE,3)'; % Library Bagging Bagging SINDy
nWrongTermsDDEmedm  = mean(nWrongTermsDDE2,3)'; % Library Bagging Bragging SINDy
modelErrorDEm       = mean(modelErrorDE,3)'; % Library Bagging 
modelErrorDDEm      = mean(modelErrorDDE,3)'; % Library Bagging Bagging SINDy
modelErrorDDEmedm   = mean(modelErrorDDE2,3)'; % Library Bagging Bragging SINDy
successDEm      = mean(successDE,3)'; % Library Bagging 
successDDEm     = mean(successDDE,3)'; % Library Bagging Bagging SINDy
successDDEmedm  = mean(successDDE2,3)'; % Library Bagging Bragging SINDy

% Jackknife Resampling
nWrongTermsDEJKm      = mean(nWrongTermsDEJK,3)'; % Library Bagging 
modelErrorDEJKm       = mean(modelErrorDEJK,3)'; % Library Bagging 
successDEJKm      = mean(successDEJK,3)'; % Library Bagging 


%% labels
legendN = 'Simulation time = ';
xlabelN = 'noise level';
ylabelN = 'data length, s';

ylimnwt = [0 2];
ylimnwtTicks = 0:0.5:2;
ylimme = [0 0.15];
ylimmeTicks = 0:0.05:0.15;
ylimSR = [0 1];
ylimSRTicks = 0:0.25:1;

pX = epsL;
pL = tEndL;
pXp = pX(1:2:end);
pLp = pL;

%% colormaps
C1 = [0 119 187]/255;
C2 = [51 187 238]/255;
C3 = [0 153 136]/255;
C4 = [238 119 51]/255;
C5 = [204 51 17]/255;
C6 = [238 51 119]/255;
C7 = [187 187 187]/255;
C8 = [80 80 80]/255;
C9 = [140 140 140]/255;
C10 = [0 128 255]/255;
yellow = [255 255 153]/255;
        
mymap = [C2
    C10
    0 0 0
    C5
    yellow];

interpMap = interp1(1:5,mymap,linspace(1,5,200));
        
m = 100;
m1 = floor(m*0.5);
r = (0:m1-1)'/max(m1,1);
g = r;
r = [r; ones(m1+1,1)];
g = [g; 1; flipud(g)];
b = flipud(r);
redblue = [r g b]; 

% colormapN = summer;
% colormapN = redblue;
% colormapN = flip(interpMap);
% colormapN = flip(inferno(100));
% colormapN = flip(viridis(100));
colormapN = flip(plasma(100));
% colormapN = flip(turbo(100));
% colormapN = flip(magma(100));


%% PLOT
nDat = 2;
posCB = 'east';
nAll = 5;

sizeX = 1200;
sizeY = 300;

figure('Position', [10 10 sizeX sizeY])

ax0 = subplot(nDat,nAll,1);
c1 = colorbar;
c1.Label.String = 'model parameter error';
c1.Ticks = ylimmeTicks;
caxis(ylimme)
c1.Location = posCB;
colormap(ax0,colormapN) 
axis off

ax1 = subplot(nDat,nAll,1+1*nAll);
c2 = colorbar;
c2.Label.String = 'success rate';
c2.Ticks = ylimSRTicks;
c2.Limits = ylimSR;
caxis(ylimSR)
c2.Location = posCB;
colormap(ax1,flipud(colormapN))
axis off

ax2 = subplot(nDat,nAll,2);
imagesc(pX,pL,modelErrorSm,ylimme)
colormap(ax2,colormapN) 
yticks(pLp)
xticks(pXp)
xticklabels([])
yticklabels([])
ax = gca;
ax.YAxisLocation = 'right';
 title('SINDy')

ax3 = subplot(nDat,nAll,2+1*nAll);
imagesc(pX,pL,successSm,ylimSR)
colormap(ax3,flipud(colormapN))
yticks(pLp)
xticks(pXp)
yticklabels([])
xlabel(xlabelN)
ax = gca;
ax.YAxisLocation = 'right';

%% bagging  
ax2 = subplot(nDat,nAll,3);
imagesc(pX,pL,modelErrorEm,ylimme)
colormap(ax2,colormapN) 
yticks(pLp)
xticks(pXp)
xticklabels([])
yticklabels([])
ax = gca;
ax.YAxisLocation = 'right';
title('Bagging')

ax3 = subplot(nDat,nAll,3+1*nAll);
imagesc(pX,pL,successEm,ylimSR)
colormap(ax3,flipud(colormapN))
yticks(pLp)
xticks(pXp)
yticklabels([])
xlabel(xlabelN)
ax = gca;
ax.YAxisLocation = 'right';

%% bragging sindy
ax2 = subplot(nDat,nAll,4);
imagesc(pX,pL,modelErrorEmedm,ylimme)
colormap(ax2,colormapN) 
yticks(pLp)
xticks(pXp)
xticklabels([])
yticklabels([])
ax = gca;
ax.YAxisLocation = 'right';
title('Bragging')

ax3 = subplot(nDat,nAll,4+1*nAll);
imagesc(pX,pL,successEmedm,ylimSR)
colormap(ax3,flipud(colormapN))
yticks(pLp)
xticks(pXp)
yticklabels([])
xlabel(xlabelN)
ax = gca;
ax.YAxisLocation = 'right';

%% double library bagging
ax2 = subplot(nDat,nAll,5);
imagesc(pX,pL,modelErrorDDEm,ylimme)
colormap(ax2,colormapN) 
yticks(pLp)
xticks(pXp)
xticklabels([])
ylabel(ylabelN)
ax = gca;
ax.YAxisLocation = 'right';
title('Lib. Bagging')

ax3 = subplot(nDat,nAll,5+1*nAll);
imagesc(pX,pL,successDDEm,ylimSR)
colormap(ax3,flipud(colormapN))
yticks(pLp)
xticks(pXp)
ylabel(ylabelN)
xlabel(xlabelN)
ax = gca;
ax.YAxisLocation = 'right';


