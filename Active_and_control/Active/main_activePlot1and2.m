%%%%%%%%%%%%%%%%%%%
% 
% plots active SINDy figure left and center
%
%

clear all
close all
clc

plotNr = 1; % plot 1 or 2

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
polys = 1:3;
trigs = [];
common_params = {polys,trigs};
gamma = 0;
tol_ode = 1e-10;         % set tolerance (abs and rel) of ode45
options = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));
Beta = cell2mat(ode_params);

% time step
dt = 0.01;

% get true Lorenz system for comparison
true_nz_weights = getTrueWeights(ode_params,common_params,n);


% signal power for noise calculation
[~,x10]=ode45(@(t,x) lorenz(t,x,Beta),dt:dt:10,x0,options);
signal_power = rms(x10(:));


%% general parameters

irng = 2;
% set rnd number for reproduction
rng(irng,'twister')

% rand condition
maxSteps = 1000;
rM = 25; % radius state space
cM = [0 0 25]; % center state space
xRc = rM*(rand(maxSteps,3)-0.5)*2 + cM;

% random initial conditions
% p = sobolset(n);
nS = 20; % number of initial points
xC = xRc(1:nS,:); % random

% Compute Derivative
dx = [];
for i=1:length(xC)
    dx(i,:) = lorenz(0,xC(i,:),Beta);
end


% add noise
eps = 0.05; % noise level
sigma = eps*signal_power;
noise = normrnd(0,sigma,[maxSteps,3]);
xobs = xC + noise(1:size(xC,1),:);

% build library
Theta_0 = build_theta(xobs,common_params);


%% SINDy
% sindy with central difference differentiation
sindy = sindy_cd(xobs,Theta_0,n,lambda,gamma,dt);

% % store outputs
% nWrongTermsS = sum(sum(abs((true_nz_weights~=0) - (sindy~=0))));
% modelErrorS = norm(sindy-true_nz_weights)/norm(true_nz_weights);
% successS = norm((true_nz_weights~=0) - (sindy~=0))==0;


%% ENSEMBLES SINDY

%% calculate derivatives
dxobs_0 = dx;

%% run active loop
iLoop = 1;
nSn = 1;
active = 1;

if plotNr == 1
    % plot 1
    plot1 = 1;
    plot2 = 0;
    nLoop = 5;
else
    % plot 2
    plot1 = 0;
    plot2 = 1;
    nLoop = 81;
end


while true
    
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
    inclProbDB = zeros(size(Theta_0,2),n);
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
        inclProbDB(libEntry,iii) = XiDBnzM;
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


    if active
    % plot variance prediction
    polysIN = 1:2; % skip last rows to oncrease speed, doesnt change results

    if plot1
        tspanE = dt:dt:2;
        x0n = xRc(nS+(iLoop-1)*nSn+(1:nSn),:);
        for iE = 1:nEnsemblesDD
            [~,xSINDYXiDB]=ode45(@(t,x)sparseGalerkin(t,x,XiDBeOut(1:end-(size(XiDB,1)-9),:,iE),polysIN),tspanE,x0n);%,options);  % approximate
            xSINDYXiDBout(:,:,iE) = xSINDYXiDB; 
        end
        maxStd = 0.5;
        for i=1:3
            xStd = reshape(xSINDYXiDBout(:,i,:),size(xSINDYXiDBout,1),size(xSINDYXiDBout,3));
            stdOut(:,i) = std(xStd',1);
            medOut(:,i) = median(xStd',1);
            mm(i)=min([min(find(stdOut(:,i)>maxStd)) max(tspanE)/dt]);
        end
        xCn = medOut(1:min(mm),:);
    end

    if plot2
        % check ambiguity
        for iRt = 1:200
        x0n = rM*(rand(nSn,3)-0.5)*2 + cM; % random
        for iE = 1:nEnsemblesDD
            xSINDYXiDB = sparseGalerkin(0,x0n',XiDBeOut(1:end-(size(XiDB,1)-9),:,iE),polysIN);
            xSINDYXiDBout(:,iE) = xSINDYXiDB; 
        end

        mmOut(iRt) = mean(std(xSINDYXiDBout',1));
        x0nOut(iRt,:) = x0n;
        end
        [~,maxmm] = max(mmOut);
        xCn = x0nOut(maxmm,:);
    end


    if plot1
        lw1 = 1.5;
        lw2 = 0.5;
        fos = 14;
        fosS = 11;
        C1 = [0 128 255]/255;
        blue = [44,127,184]./255;
        % green: y
        green = [49,163,84]./255;
        % orange: z
        orange = [240,59,32]./255;
        colorsNew = [blue; green; orange];
        tEnd = 1;
        greyf = 0.85;
        if iLoop == 1
            sizeX = 400;
            sizeY = 600;
            figure('Position', [10 10 sizeX sizeY])
        end
        subplot(nLoop,3,(iLoop-1)*3+1)
        plot(tspanE,reshape(xSINDYXiDBout(:,1,:),size(xSINDYXiDBout,1),size(xSINDYXiDBout,3)),'Color',colorsNew(1,:),'Linewidth',lw2); hold on
        plot(tspanE,medOut(:,1),'--','Color',[1 1 1]*greyf,'Linewidth',lw1); hold on
        set(gca,'ticklabelinterpreter','latex','FontSize',fosS)
        xlim([0 tEnd])
        yticklabels([])
        if iLoop ~= nLoop
            xticklabels([])
        end
        if iLoop == 1
            title('x','interpreter','latex','FontSize',fos) 
        end
        ylabel(sprintf('IC %d',iLoop),'interpreter','latex','FontSize',fos)
        subplot(nLoop,3,(iLoop-1)*3+2)
        plot(tspanE,reshape(xSINDYXiDBout(:,2,:),size(xSINDYXiDBout,1),size(xSINDYXiDBout,3)),'Color',colorsNew(2,:),'Linewidth',lw2); hold on
        plot(tspanE,medOut(:,2),'--','Color',[1 1 1]*greyf,'Linewidth',lw1); hold on
        set(gca,'ticklabelinterpreter','latex','FontSize',fosS)
        xlim([0 tEnd])
        if iLoop ~= nLoop
            xticklabels([])
        end
        yticklabels([])
        if iLoop == 1
            title('y','interpreter','latex','FontSize',fos) 
        end
        if iLoop == nLoop
            xlabel('Time, s','interpreter','latex','FontSize',fos)
        end
        subplot(nLoop,3,(iLoop-1)*3+3)
        if iLoop == nLoop
            % ghost line for legend
            plot([-10 -9 ],[25 25],'k','Linewidth',0.5); hold on
            plot([-10 -9 ],[25 25],'--','Color',[1 1 1]*greyf,'Linewidth',1); hold on
        end
        plot(tspanE,reshape(xSINDYXiDBout(:,3,:),size(xSINDYXiDBout,1),size(xSINDYXiDBout,3)),'Color',colorsNew(3,:),'Linewidth',lw2); hold on
        plot(tspanE,medOut(:,3),'--','Color',[1 1 1]*greyf,'Linewidth',lw1); hold on
        set(gca,'ticklabelinterpreter','latex','FontSize',fosS)
        xlim([0 tEnd])
        if iLoop ~= nLoop
            xticklabels([])
        end
        yticklabels([])
        if iLoop == 1
            title('z','interpreter','latex','FontSize',fos) 
        end
        if iLoop == nLoop
            legend({'ensemble','median'},'interpreter','latex','FontSize',fos)
        end
    end

    else % random
        xCn = xRc(nS+(iLoop-1)*nSn+(1:nSn),:);
    end

    xnobs = xCn + noise(nS+(iLoop-1)*nSn+(1:nSn),:);
    dxn = [];
    for i=1:size(xCn,1)
        dxn(i,:) = lorenz(0,xCn(i,:),Beta);
    end

    xobs = [xobs; xnobs];
    dx = [dx; dxn];
    dxobs_0 = dx;

    Theta_0 = build_theta(xobs,common_params);

    if plot2
        if iLoop == 1
            lib = poolDataLIST({'x','y','z'},XiDB,n,polys);    
            lib(1) = {'x '};
            lib(2) = {'y '};
            lib(3) = {'z '};
            skipLastRows = size(XiDB,1)-9; % speed up sampling by reducing size of library
            XL = [];
            XL = plotUQ_LorenzActive(XiDBeOut(1:end-skipLastRows,:,:),true_nz_weights(1:end-skipLastRows,:),XiDB,lib,iLoop,XL);
        end
    end

    if iLoop == nLoop
        break
    end

    iLoop = iLoop + 1;

end


if plot2
    lib = poolDataLIST({'x','y','z'},XiDB,n,polys);    
    lib(1) = {'x '};
    lib(2) = {'y '};
    lib(3) = {'z '};
    skipLastRows = size(XiDB,1)-9; % speed up sampling by reducing size of library
    plotUQ_LorenzActive(XiDBeOut(1:end-skipLastRows,:,:),true_nz_weights(1:end-skipLastRows,:),XiDB,lib,iLoop,XL)
end

