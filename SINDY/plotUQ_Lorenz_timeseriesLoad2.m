function plotUQ_Lorenz_timeseriesLoad2(results,pct,nE,tspan,colors)


xSINDYXi = results.xSINDYXi;
tspanSINDyXi = results.tspanSINDyXi;
xTrue = results.xTrue;
tTrue = results.tTrue;
xSINDY = results.xSINDY;
tspanSINDy = results.tspanSINDy;

prcTT1 = results.prcTT1;
prcTT2 = results.prcTT2;
prcTT3 = results.prcTT3;

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

lw = 1.2;
lwlb = 1.2;
lwlbA = 1.0;

lst = '--'; % line style true
lss = ':'; % line style sindy
lslb = '-'; % line style library bagging
lslb1 = '-';

Ct = C9;%C7;
Cs = C5;
Clb = C10;
ClbA = Clb;

mymap3 = [255,237,160
    254,178,76
    240,59,32]./255;
Cs = mymap3(3,:);

% new colors

Clb = colors;
ClbA = Clb;
Cs = 'k';
Ct = 'k';


faceAlpha = 0.3;

sizeX = 440;
sizeY = 500;

fos = 16;
fos2 = 14;
fosT = 14;


figure('Position', [10 10 sizeX sizeY])

subplot(3,1,1)

% ghost plot for legend
plot([0 1],[-1000 -1000],lst,'Color',Ct,'LineWidth',lw); hold on
plot([0 1],[-1000 -1000],lss,'Color',Cs,'LineWidth',lw); hold on
plot([0 1],[-1000 -1000],lslb,'Color',Clb(1,:),'LineWidth',lwlb); hold on


% actual plot
patch([tspan'; flipud(tspan')], [prcTT1(:,1); flipud(prcTT1(:,3))], ClbA(1,:), 'FaceAlpha',faceAlpha, 'EdgeColor','none'); hold on
plot(tTrue,xTrue(:,1),lst,'Color',Ct,'LineWidth',lw); hold on
plot(tspanSINDy,xSINDY(:,1),lss,'Color',Cs,'LineWidth',lw); hold on
plot(tspanSINDyXi,xSINDYXi(:,1),lslb,'Color',Clb(1,:),'LineWidth',lwlb); hold on
% ylim([-25 45])
ylim([-25 25])
xlim([0 5])
xticklabels([])
set(gca,'ticklabelinterpreter','latex','FontSize',fos2)
ylabel('x','interpreter','latex','FontSize',fos)
box on

subplot(3,1,2)
patch([tspan'; flipud(tspan')], [prcTT2(:,1); flipud(prcTT2(:,3))], ClbA(2,:), 'FaceAlpha',faceAlpha, 'EdgeColor','none'); hold on
plot(tTrue,xTrue(:,2),lst,'Color',Ct,'LineWidth',lw); hold on
plot(tspanSINDy,xSINDY(:,2),lss,'Color',Cs,'LineWidth',lw); hold on
plot(tspanSINDyXi,xSINDYXi(:,2),lslb,'Color',Clb(2,:),'LineWidth',lwlb); hold on
ylim([-25 25])
xlim([0 5])
set(gca,'ticklabelinterpreter','latex','FontSize',fos2)
ylabel('y','interpreter','latex','FontSize',fos)
xticklabels([])
box on

subplot(3,1,3)

% ghost plot for legend
plot([0 1],[-1000 -1000],lst,'Color',Ct,'LineWidth',lw); hold on
plot([0 1],[-1000 -1000],lss,'Color',Cs,'LineWidth',lw); hold on
plot([0 1],[-1000 -1000],lslb,'Color',Clb(3,:),'LineWidth',lwlb); hold on

patch([tspan'; flipud(tspan')], [prcTT3(:,1); flipud(prcTT3(:,3))], ClbA(3,:), 'FaceAlpha',faceAlpha, 'EdgeColor','none'); hold on
plot(tTrue,xTrue(:,3),lst,'Color',Ct,'LineWidth',lw); hold on
plot(tspanSINDy,xSINDY(:,3),lss,'Color',Cs,'LineWidth',lw); hold on
plot(tspanSINDyXi,xSINDYXi(:,3),lslb,'Color',Clb(3,:),'LineWidth',lwlb); hold on
ylim([0 50])
xlim([0 5])
yticks(0:25:50)
set(gca,'ticklabelinterpreter','latex','FontSize',fos2)
ylabel('z','interpreter','latex','FontSize',fos)
xlabel('time, s','interpreter','latex','FontSize',fos)
box on


