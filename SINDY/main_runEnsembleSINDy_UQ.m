%%%%%%%%%%%%%%%%%%%
% 
% run sims for one noise level and data length to plot ensemble forecasting and UQ
%
%

clear all
close all
clc

%% sweep over a set of noise levels and data length to generate heatmap plots
% noise level
eps = 0.025;

% simulation time
tEnd = 10;

% at each noise level and simulation time, nTest different instantiations of noise are run (model errors and success rate are then averaged for plotting)
nTest1 = 1; % generate models nTest1 times for SINDy
nTest2 = 1; % generate models nTest times for ensemble SINDy


%% hyperparameters
% SINDy sparsifying hyperparameters
lambda = 0.2;

% ensemble hyperparameters
% data ensembling
nEnsembles = 100; % number of bootstraps (SINDy models using sampled data) in ensemble
ensembleT = 0.6; % Threshold model coefficient inclusion probability: set ensemble SINDy model coefficient to zero if inclusion probability is below ensembleT

% library
nEnsemble1P = 0.9; % percentage of full library that is sampled without replacement for library bagging
nEnsemble2 = 100; % number of bootstraps (SINDy models using sampled library terms) in ensemble
ensT = 0.4; % Threshold library term inclusion probabilty: cut library entries that occur less than ensT

% double bagging
nEnsemblesDD = 100; % number of models in ensemble for data bagging after library bagging


%% common parameters, true Lorenz system, signal power for noise calculation

% generate synthetic Lorenz system data
ode_params = {10, 8/3, 28}; 
x0 = [-8 7 27]';
n = length(x0); 

% set common params
polys = 1:5;
trigs = [];
common_params = {polys,trigs};
gamma = 0;
tol_ode = 1e-10;         % set tolerance (abs and rel) of ode45
options = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));
Beta = cell2mat(ode_params);

% time step
dt = 0.01;
tspan = dt:dt:tEnd;

% get true Lorenz system for comparison
true_nz_weights = getTrueWeights(ode_params,common_params,n);

% signal power for noise calculation
[~,x10]=ode45(@(t,x) lorenz(t,x,Beta),dt:dt:10,x0,options);
signal_power = rms(x10(:));


%% general parameters

% smooth data using golay filter 
sgolayON = 1;

% generate data
[t,x]=ode45(@(t,x) lorenz(t,x,Beta),tspan,x0,options);

% set rnd number for reproduction
rng(1,'twister')

% add noise
sigma = eps*signal_power;
noise = normrnd(0,sigma,size(x));
xobs = x + noise;

% data before smoothing for plotting
xobsPlotE = xobs;

% smooth data
if sgolayON 
    order = 3;
    framelen = 5;
    xobs = sgolayfilt(xobs,order,framelen);
end

% build library
Theta_0 = build_theta(xobs,common_params);


%% SINDy
% sindy with central difference differentiation
sindy = sindy_cd(xobs,Theta_0,n,lambda,gamma,dt);

% store outputs
nWrongTermsS = sum(sum(abs((true_nz_weights~=0) - (sindy~=0))));
modelErrorS = norm(sindy-true_nz_weights)/norm(true_nz_weights);
successS = norm((true_nz_weights~=0) - (sindy~=0))==0;


%% ENSEMBLES SINDY

%% calculate derivatives
% finite difference differentiation
dxobs_0 = zeros(size(x));
dxobs_0(1,:)=(-11/6*xobs(1,:) + 3*xobs(2,:) -3/2*xobs(3,:) + xobs(4,:)/3)/dt;
dxobs_0(2:size(xobs,1)-1,:) = (xobs(3:end,:)-xobs(1:end-2,:))/(2*dt);
dxobs_0(size(xobs,1),:) = (11/6*xobs(end,:) - 3*xobs(end-1,:) + 3/2*xobs(end-2,:) - xobs(end-3,:)/3)/dt;
            

%% double bagging SINDy
            
%% Bagging SINDy library
% randomly sample library terms without replacement and throw away terms
% with low inclusion probability
nEnsemble1 = round(nEnsemble1P*size(Theta_0,2));
mOutBS = zeros(nEnsemble1,n,nEnsemble2);
libOutBS = zeros(nEnsemble1,nEnsemble2);
for iii = 1:nEnsemble2
    rs = RandStream('mlfg6331_64','Seed',iii); 
    libOutBS(:,iii) = datasample(rs,1:size(Theta_0,2),nEnsemble1,'Replace',false)';
    mOutBS(:,:,iii) = sparsifyDynamics(Theta_0(:,libOutBS(:,iii)),dxobs_0,lambda,n,gamma);
end

inclProbBS = zeros(size(Theta_0,2),n);
for iii = 1:nEnsemble2
    for jjj = 1:n
        for kkk = 1:nEnsemble1
            if mOutBS(kkk,jjj,iii) ~= 0
                inclProbBS(libOutBS(kkk,iii),jjj) = inclProbBS(libOutBS(kkk,iii),jjj) + 1;
            end
        end
    end
end
inclProbBS = inclProbBS/nEnsemble2*size(Theta_0,2)/nEnsemble1;

XiD = zeros(size(Theta_0,2),n);
for iii = 1:n
    libEntry = inclProbBS(:,iii)>ensT;
    XiBias = sparsifyDynamics(Theta_0(:,libEntry),dxobs_0(:,iii),lambda,1,gamma);
    XiD(libEntry,iii) = XiBias;
end

                
%% Double bagging SINDy 
% randomly sample library terms without replacement and throw away terms
% with low inclusion probability
% then on smaller library, do bagging

XiDB = zeros(size(Theta_0,2),n);
XiDBmed = zeros(size(Theta_0,2),n);
XiDBs = zeros(size(Theta_0,2),n);
XiDBeOut = zeros(size(Theta_0,2),n,nEnsemblesDD);
for iii = 1:n
    libEntry = inclProbBS(:,iii)>ensT;

    bootstatDD = bootstrp(nEnsemblesDD,@(Theta,dx)sparsifyDynamics(Theta,dx,lambda,1,gamma),Theta_0(:,libEntry),dxobs_0(:,iii)); 

    XiDBe = [];
    XiDBnz = [];
    for iE = 1:nEnsemblesDD
        XiDBe(:,iE) = reshape(bootstatDD(iE,:),size(Theta_0(:,libEntry),2),1);
        XiDBnz(:,iE) = XiDBe(:,iE)~=0;
        
        XiDBeOut(libEntry,iii,iE) = XiDBe(:,iE);
    end

    % Thresholded bootstrap aggregating (bagging, from bootstrap aggregating)
    XiDBnzM = mean(XiDBnz,2); % mean of non-zero values in ensemble
    XiDBnzM(XiDBnzM<ensembleT) = 0; % threshold: set all parameters that have an inclusion probability below threshold to zero

    XiDBmean = mean(XiDBe,2);
    XiDBmedian = median(XiDBe,2);
    XiDBstd = std(XiDBe')';

    XiDBmean(XiDBnzM==0)=0; 
    XiDBmedian(XiDBnzM==0)=0; 
    XiDBstd(XiDBnzM==0)=0; 
    
    XiDB(libEntry,iii) = XiDBmean;
    XiDBmed(libEntry,iii) = XiDBmedian;
    XiDBs(libEntry,iii) = XiDBstd;
    
end


%% model error and success rates

nWrongTermsDE = sum(sum(abs((true_nz_weights~=0) - (XiD~=0))));
modelErrorDE = norm(XiD-true_nz_weights)/norm(true_nz_weights);
successDE = norm((true_nz_weights~=0) - (XiD~=0))==0;

nWrongTermsDDE = sum(sum(abs((true_nz_weights~=0) - (XiDB~=0))));
nWrongTermsDDE2 = sum(sum(abs((true_nz_weights~=0) - (XiDBmed~=0))));
modelErrorDDE = norm(XiDB-true_nz_weights)/norm(true_nz_weights);
modelErrorDDE2 = norm(XiDBmed-true_nz_weights)/norm(true_nz_weights);
successDDE = norm((true_nz_weights~=0) - (XiDB~=0))==0;
successDDE2 = norm((true_nz_weights~=0) - (XiDBmed~=0))==0;



%% model error and success rate: SINDy, library bagging, double library bagging
modelError = [modelErrorS, modelErrorDE, modelErrorDDE]
successRate = [successS, successDE, successDDE]


%% compare time series for different ensemble methods

polysIN = 1:2; % skip last rows to oncrease speed, doesnt change results
skipLastRows = size(XiDB,1)-9; 

% sindy
Xi = sindy(1:end-skipLastRows,:);
[tspanSINDy,xSINDY]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polysIN),tspan,x0);%,options);  % approximate

% Library bagging 
Xi = XiD(1:end-skipLastRows,:);
[tspanSINDyD,xSINDYXiD]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polysIN),tspan,x0);%,options);  % approximate

% Double bagging 
Xi = XiDB(1:end-skipLastRows,:);
[tspanSINDyDB,xSINDYXiDB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polysIN),tspan,x0);%,options);  % approximate

[tTrue,xTrue]=ode45(@(t,x) lorenz(t,x,Beta),tspan,x0,options);


%% plot
% color propreties
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

lslb = ':'; % line style library bagging

Ct = C9;
Cs = C5;
Clb = C2;
Clb1 = C10;

sizeX = 600;
sizeY = 500;

figure('Position', [10 10 sizeX sizeY])
subplot(3,1,1)
plot(tspan,xobsPlotE(:,1),'Color',Ct,'LineWidth',lw); hold on
plot(tspanSINDy,xSINDY(:,1),'Color',Cs,'LineWidth',lw); hold on
plot(tspanSINDyDB,xSINDYXiDB(:,1),lslb,'Color',Clb,'LineWidth',lwlb); hold on
ylim([-25 45])
legend({'observed dynamics','SINDy','LB-SINDy'},'NumColumns',3)

subplot(3,1,2)
plot(tspan,xobsPlotE(:,2),'Color',Ct,'LineWidth',lw); hold on
plot(tspanSINDy,xSINDY(:,2),'Color',Cs,'LineWidth',lw); hold on
plot(tspanSINDyDB,xSINDYXiDB(:,2),lslb,'Color',Clb,'LineWidth',lwlb); hold on
ylim([-25 25])

subplot(3,1,3)
plot(tspan,xobsPlotE(:,3),'Color',Ct,'LineWidth',lw); hold on
plot(tspanSINDy,xSINDY(:,3),'Color',Cs,'LineWidth',lw); hold on
plot(tspanSINDyDB,xSINDYXiDB(:,3),lslb,'Color',Clb,'LineWidth',lwlb); hold on
ylim([0 50])
sgtitle('Reconstruction: SINDy vs. Library Bagging')



%% plot UQ time series: ensemble forecast (using standard SINDy and double library bagging)
% draw multiple models and average -> bagging forecast (goes to mean forecast for large nE) 

x0 = xobsPlotE(end,:)'; % final point of training data
tspan = 0.01:0.01:5;
nUQ = 1000; % number of samples for UQ
% nE is the number of models use for ensemble forecasting: e.g. nE=10 uses 10 SINDy model and averages the model coefficients
for nE = [1 10 100] % number of ensembles for forecast
    Em = 1; % ensemble forecast method: 1) run models in parallel and take mean of prediction, or 2) take mean of model parameters and run single model
    pct = 95; % plot pct% confidence interval 
    skipLastRows = size(XiDB,1)-9; % speed up sampling by reducing size of library
    polysIN = 1:2; % dont change, sparseGalerkin.m function is optimised for polys = 1:2
    Beta = cell2mat(ode_params);

    plotUQ_Lorenz_timeseries(XiDB(1:end-skipLastRows,:),XiDBs(1:end-skipLastRows,:),x0,tspan,Beta,polysIN,nUQ,pct,nE,sindy(1:end-skipLastRows,:),Em,options)
end



%% plot uncertainty in coefficients

lib = poolDataLIST({'x','y','z'},XiDB,n,polys);    

skipLastRows = size(XiDB,1)-9; % speed up sampling by reducing size of library
plotUQ_Lorenz(XiDBeOut(1:end-skipLastRows,:,:),true_nz_weights(1:end-skipLastRows,:),XiDB,lib)

