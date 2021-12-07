
%% Lotka Volterra lynx hare experimental data: ensemble SINDy 

clear all
close all
clc

% lynx and hare population data from http://www.math.tamu.edu/~phoward/m442/modbasics.pdf
% year, lynx, hare
lhpop = [1900,1901,1902,1903,1904,1905,1906,1907,1908,1909,1910,1911,1912,1913,1914,1915,1916,1917,1918,1919,1920;4,6.10000000000000,9.80000000000000,35.2000000000000,59.4000000000000,41.7000000000000,19,13,8.30000000000000,9.10000000000000,7.40000000000000,8,12.3000000000000,19.5000000000000,45.7000000000000,51.1000000000000,29.7000000000000,15.8000000000000,9.70000000000000,10.1000000000000,8.60000000000000;30,47.2000000000000,70.2000000000000,77.4000000000000,36.3000000000000,20.6000000000000,18.1000000000000,21.4000000000000,22,25.4000000000000,27.1000000000000,40.3000000000000,57,76.6000000000000,52.3000000000000,19.5000000000000,11.2000000000000,7.60000000000000,14.6000000000000,16.2000000000000,24.7000000000000];

tspan = lhpop(1,:)-lhpop(1,1);
xobs = lhpop([3 2],:)'./std(lhpop([3 2],:)');

% true system parameter estimation Seth Hirsh UQ-SINDy paper
a = 0.55;
b = 0.455; 
d = 0.5433;
g = 0.84;  
true_nz_weights = zeros(10,2);
true_nz_weights(2,1) = a;
true_nz_weights(4,1) = -b;
true_nz_weights(3,2) = -g;
true_nz_weights(4,2) = d;


% common parameters
n = 2;
polys = 0:3; % if changed, also change sparseGalerkin.m function, as it is optimised for polys = 0:3
gamma = 0;
common_params = {polys,[]};

tol_ode = 1e-10;         % set tolerance (abs and rel) of ode45
options = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,n));

Theta = build_theta(xobs,common_params);

%% calculate derivatives
dtL = 1;

% Fourth order centered difference with third order forward/backward difference at endpoints.
dxobs(1,:)=(-11/6*xobs(1,:) + 3*xobs(2,:) -3/2*xobs(3,:) + xobs(4,:)/3)/1;
dxobs(2,:)=(-11/6*xobs(2,:) + 3*xobs(3,:) -3/2*xobs(4,:) + xobs(5,:)/3)/1;
dxobs(3:19,:) = (-1/12*xobs(5:end,:) + 2/3*xobs(4:end-1,:) - 2/3*xobs(2:end-3,:) + 1/12*xobs(1:end-4,:));
dxobs(20,:) = (11/6*xobs(end-1,:) - 3*xobs(end-2,:) + 3/2*xobs(end-3,:) - xobs(end-4,:)/3)/1;
dxobs(21,:) = (11/6*xobs(end,:) - 3*xobs(end-1,:) + 3/2*xobs(end-2,:) - xobs(end-3,:)/3)/1;


nEnsemble2 = 150;
ensT = 0.65;
nEnsemble1P = 0.85;
ensembleT = 0.8;
nEnsemblesDD = 1000; % larger ensemble for refined UQ if using plotUQtimeseriesELVbootstrap
lambda = 0.19;


%% Bagging SINDy library
% randomly sample library terms without replacement and throw away terms
% with low inclusion probability
nEnsemble1 = round(nEnsemble1P*size(Theta,2));
mOutBS = zeros(nEnsemble1,n,nEnsemble2);
libOutBS = zeros(nEnsemble1,nEnsemble2);
for iii = 1:nEnsemble2
    rs = RandStream('mlfg6331_64','Seed',iii); 
    libOutBS(:,iii) = datasample(rs,1:size(Theta,2),nEnsemble1,'Replace',false)';
    mOutBS(:,:,iii) = sparsifyDynamics(Theta(:,libOutBS(:,iii)),dxobs,lambda,n,gamma);
end

inclProbBS = zeros(size(Theta,2),n);
for iii = 1:nEnsemble2
    for jjj = 1:n
        for kkk = 1:nEnsemble1
            if mOutBS(kkk,jjj,iii) ~= 0
                inclProbBS(libOutBS(kkk,iii),jjj) = inclProbBS(libOutBS(kkk,iii),jjj) + 1;
            end
        end
    end
end
inclProbBS = inclProbBS/nEnsemble2*size(Theta,2)/nEnsemble1;

XiD = zeros(size(Theta,2),n);
for iii = 1:n
    libEntry = inclProbBS(:,iii)>ensT;
    XiBias = sparsifyDynamics(Theta(:,libEntry),dxobs(:,iii),lambda,1,gamma);
    XiD(libEntry,iii) = XiBias;
end


%% Double bagging SINDy 
% randomly sample library terms without replacement and throw away terms with low inclusion probability
% then, on smaller library, do bagging

XiDB = zeros(size(Theta,2),n);
inclProbDB = zeros(size(Theta,2),n);
XiDBs = zeros(size(Theta,2),n);
XiDBeOut = zeros(size(Theta,2),n,nEnsemblesDD);
XiDBeOut2 = zeros(size(Theta,2),n,nEnsemblesDD);

OOSsmallOut = zeros(nEnsemblesDD,n);

for iii = 1:n
    libEntry = inclProbBS(:,iii)>ensT;
    bootstatDD = bootstrp(nEnsemblesDD,@(Theta,dx)sparsifyDynamics(Theta,dx,lambda,1,gamma),Theta(:,libEntry),dxobs(:,iii)); 

    XiDBe = [];
    XiDBnz = [];
    for iE = 1:nEnsemblesDD
        XiDBe(:,iE) = reshape(bootstatDD(iE,:),size(Theta(:,libEntry),2),1);
        XiDBnz(:,iE) = XiDBe(:,iE)~=0;
        
        XiDBeOut(libEntry,iii,iE) = XiDBe(:,iE);
    end

    % Thresholded bootstrap aggregating (bagging, from bootstrap aggregating)
    XiDBnzM = mean(XiDBnz,2); % mean of non-zero values in ensemble
    inclProbDB(libEntry,iii) = XiDBnzM;
    XiDBnzM(XiDBnzM<ensembleT) = 0; % threshold: set all parameters that have an inclusion probability below threshold to zero

    XiDBmean = mean(XiDBe,2);

    XiDBmean(XiDBnzM==0)=0; 

    XiDB(libEntry,iii) = XiDBmean;
    
    XiDBstd = std(XiDBe')';
    XiDBstd(XiDBnzM==0)=0; 
    XiDBs(libEntry,iii) = XiDBstd;

    for iE = 1:nEnsemblesDD
        XiDBeOutInter = XiDBe(:,iE);
        XiDBeOutInter(XiDBnzM==0)=0; 
        XiDBeOut2(libEntry,iii,iE) = XiDBeOutInter;
    end
end


%% plot UQ time series: ensemble forecast
% draw not from each coefficient probability, but draw 5 different models from the ensemble and average
% we can compare this to the posterior predictive distribution (PPD) from Hirsh 2021.
nUQ = size(XiDBeOut,3);
nE = 5; % number of ensembles for forecast
pct = 95; % plot prctile 
plotUQ_LV_timeseries(XiDB,XiDBeOut2,XiDBs,xobs(1,:),tspan,polys,nUQ,pct,nE,tspan,xobs,options,lhpop)

%% plot uncertainty in coefficients
lib = {'1 ';'u ';'v ';'uv';'vv';'uu'};
XiDBeOutFigure = XiDBeOut(1:6,:,:);
plotUQ_LV(XiDBeOutFigure,true_nz_weights,XiDB,lib)


%% inclusion probability
% multiply inclusion probability of library bagging with second
% inclusion probability of data bagging
inclProb1 = inclProbBS;
inclProb2a = inclProbDB;
inclProb2b = inclProbBS;
inclProb2b(inclProbDB~=0) = inclProb2b(inclProbDB~=0).*inclProbDB(inclProbDB~=0);


