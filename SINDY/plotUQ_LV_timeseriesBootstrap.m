function plotUQ_LV_timeseriesBootstrap(Xi,Xibs,Xistd,x0,tspan,polys,nUQ,pct,nE,tobs,xobs,options,lhpop)

rng(1,'twister')

iEE = 1;
che = 10; % check if sim doverges, then don't add this

for iE = 1:nUQ
    
    XiUQ = Xibs(:,:,iE);
    
    [~,xSINDYb]=ode45(@(t,x)sparseGalerkin(t,x,XiUQ,polys),tspan,x0,options);  % approximate

    if max(max(xSINDYb)) < che
        estimateUQ1(1:size(xSINDYb,1),iEE) = xSINDYb(:,1);
        estimateUQ2(1:size(xSINDYb,1),iEE) = xSINDYb(:,2);
        iEE = iEE + 1;
    end
end

% ensemble SINDy mean
[tspanSINDyXi,xSINDYXi]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polys),tspan,x0,options);  % approximate

pctL = (100-pct)/2;
pctH = (100+pct)/2;
prcTT1 = prctile(estimateUQ1,[pctL 50 pctH],2);
prcTT2 = prctile(estimateUQ2,[pctL 50 pctH],2);


%% plot combined

% color
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

lw = 1.0;
lwlb = 1.0;

Ct = C9;
Cs = C5;
Clb = C10;
ClbA = Clb;


cLynx = 'b';%Cs;
cHare = 'r';%Clb;
lslb = '-';

faceAlpha = 0.2;
ms = 8;
cGrey = 0.6*[1 1 1];

y0 = lhpop(1,1);

tTrue = tobs + y0;
tspanSINDyXi = tspanSINDyXi + y0;
tspan = tspan + y0;

sizeX = 600;
sizeY = 300;
figure('Position', [10 10 sizeX sizeY])

% ghost plot legend
plot([0 1],[-1000 -1000],'-','Color',cHare,'LineWidth',lwlb); hold on
plot([0 1],[-1000 -1000],'-','Color',cLynx,'LineWidth',lwlb); hold on
plot([0 1],[-1000 -1000],'x','Color',cGrey,'LineWidth',lw,'MarkerSize',ms); hold on
plot([0 1],[-1000 -1000],lslb,'Color',cGrey,'LineWidth',lwlb); hold on
patch([tspan'-1000; flipud(tspan'-1000)], [prcTT1(:,1); flipud(prcTT1(:,3))], cGrey, 'FaceAlpha',faceAlpha, 'EdgeColor','none'); hold on

% actual plot
n1 = std(lhpop(3,:)'); % scale output
n2 = std(lhpop(2,:)');
patch([tspan'; flipud(tspan')], [prcTT1(:,1)*n1; flipud(prcTT1(:,3)*n1)], cHare, 'FaceAlpha',faceAlpha, 'EdgeColor','none'); hold on
plot(tTrue,xobs(:,1)*n1,'x','Color',cHare,'LineWidth',lw,'MarkerSize',ms); hold on
plot(tspanSINDyXi,xSINDYXi(:,1)*n1,lslb,'Color',cHare,'LineWidth',lwlb); hold on
patch([tspan'; flipud(tspan')], [prcTT2(:,1)*n2; flipud(prcTT2(:,3)*n2)], cLynx, 'FaceAlpha',faceAlpha, 'EdgeColor','none'); hold on
plot(tTrue,xobs(:,2)*n2,'x','Color',cLynx,'LineWidth',lw,'MarkerSize',ms); hold on
plot(tspanSINDyXi,xSINDYXi(:,2)*n2,lslb,'Color',cLynx,'LineWidth',lwlb); hold on

ylim([-5 165])
yticks([0:40:160])
ylabel('Population, thousands')
box on
xlabel('Time, years')
xticks([1900:5:1920])
xlim([1899.5 1920.5])

legend({'Hare','Lynx','Obs. data','LB-SINDy',sprintf('%d%% conf.',pct)},'Location','NorthEast','NumColumns',1)

title('Observed vs. Library Bagging SINDy time series')
