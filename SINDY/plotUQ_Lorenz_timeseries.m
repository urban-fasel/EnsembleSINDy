function plotUQ_Lorenz_timeseries(Xi,Xistd,x0,tspan,Beta,polys,nUQ,pct,nE,XiS,Em,options)


estimateUQ1 = zeros(length(tspan),nUQ);
estimateUQ2 = zeros(length(tspan),nUQ);
estimateUQ3 = zeros(length(tspan),nUQ);

rng(1,'twister')

for iE = 1:nUQ

    if nE == 1
        XiUQ = normrnd(Xi,Xistd);
    else
        XiUQ = zeros(size(Xi,1),size(Xi,2),nE);
        for iE2 = 1:nE
            XiUQ(:,:,iE2) = normrnd(Xi,Xistd);
        end
        if Em == 1
            XiUQ = mean(XiUQ,3);
        end
    end
    
    if Em == 1
        [~,xSINDYb]=ode45(@(t,x)sparseGalerkin(t,x,XiUQ,polys),tspan,x0,options);  % approximate
    else
        for iE2 = 1:nE
            [~,xSINDYbiE]=ode45(@(t,x)sparseGalerkin(t,x,XiUQ(:,:,iE2),polys),tspan,x0,options);  % approximate
            xSINDYbiEOut(:,:,iE2) = xSINDYbiE;
        end
        xSINDYb = mean(xSINDYbiEOut,3);
    end

    estimateUQ1(1:size(xSINDYb,1),iE) = xSINDYb(:,1);
    estimateUQ2(1:size(xSINDYb,1),iE) = xSINDYb(:,2);
    estimateUQ3(1:size(xSINDYb,1),iE) = xSINDYb(:,3);

end


[tspanSINDyXi,xSINDYXi]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polys),tspan,x0,options);  % approximate

[tTrue,xTrue]=ode45(@(t,x) lorenz(t,x,Beta),tspan,x0,options);

% sindy
[tspanSINDy,xSINDY]=ode45(@(t,x)sparseGalerkin(t,x,XiS,polys),tspan,x0,options);  % approximate

pctL = (100-pct)/2;
pctH = (100+pct)/2;
prcTT1 = prctile(estimateUQ1,[pctL 50 pctH],2);
prcTT2 = prctile(estimateUQ2,[pctL 50 pctH],2);
prcTT3 = prctile(estimateUQ3,[pctL 50 pctH],2);

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
lwlbA = 1.0;

lss = '-'; % line style sindy
lslb = '-'; % line style library bagging
lslb1 = '-';

Ct = C9;
Cs = C5;
Clb = C10;
ClbA = Clb;

faceAlpha = 0.3;

sizeX = 600;
sizeY = 500;

figure('Position', [10 10 sizeX sizeY])

subplot(3,1,1)

% ghost plot for legend
plot([0 1],[-1000 -1000],'Color',Ct,'LineWidth',lw); hold on
plot([0 1],[-1000 -1000],lss,'Color',Cs,'LineWidth',lw); hold on
plot([0 1],[-1000 -1000],lslb,'Color',Clb,'LineWidth',lwlb); hold on


% actual plot
patch([tspan'; flipud(tspan')], [prcTT1(:,1); flipud(prcTT1(:,3))], ClbA, 'FaceAlpha',faceAlpha, 'EdgeColor','none'); hold on
plot(tTrue,xTrue(:,1),'Color',Ct,'LineWidth',lw); hold on
plot(tspanSINDy,xSINDY(:,1),lss,'Color',Cs,'LineWidth',lw); hold on
plot(tspanSINDyXi,xSINDYXi(:,1),lslb,'Color',Clb,'LineWidth',lwlb); hold on
ylim([-25 45])
xticklabels([])
ylabel('x')
box on
legend({'True dyn.','SINDy','LB-SINDy',sprintf('%d%% conf.',pct)},'Location','NorthEast','NumColumns',2)
% legend('boxoff')

subplot(3,1,2)
patch([tspan'; flipud(tspan')], [prcTT2(:,1); flipud(prcTT2(:,3))], ClbA, 'FaceAlpha',faceAlpha, 'EdgeColor','none'); hold on
plot(tTrue,xTrue(:,2),'Color',Ct,'LineWidth',lw); hold on
plot(tspanSINDy,xSINDY(:,2),lss,'Color',Cs,'LineWidth',lw); hold on
plot(tspanSINDyXi,xSINDYXi(:,2),lslb,'Color',Clb,'LineWidth',lwlb); hold on
ylim([-25 25])
ylabel('y')
xticklabels([])
box on

subplot(3,1,3)

% ghost plot for legend
plot([0 1],[-1000 -1000],'Color',Ct,'LineWidth',lw); hold on
plot([0 1],[-1000 -1000],lss,'Color',Cs,'LineWidth',lw); hold on
plot([0 1],[-1000 -1000],lslb,'Color',Clb,'LineWidth',lwlb); hold on

patch([tspan'; flipud(tspan')], [prcTT3(:,1); flipud(prcTT3(:,3))], ClbA, 'FaceAlpha',faceAlpha, 'EdgeColor','none'); hold on
plot(tTrue,xTrue(:,3),'Color',Ct,'LineWidth',lw); hold on
plot(tspanSINDy,xSINDY(:,3),lss,'Color',Cs,'LineWidth',lw); hold on
plot(tspanSINDyXi,xSINDYXi(:,3),lslb,'Color',Clb,'LineWidth',lwlb); hold on
ylim([0 50])
yticks(0:25:50)
ylabel('z')
xlabel('time, s')
box on

sgtitle('Forecast: SINDy vs. Library Bagging')

