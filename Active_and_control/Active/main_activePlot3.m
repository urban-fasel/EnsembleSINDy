%%%%%%%%%%%%%%%%%%%
% 
% plots active SINDy figure right
%
%

clear all
close all
clc

runSweep = 0; % run (1) or load (0) results
active = 1; % run active (1) or passive (0)

if runSweep
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


    %% main loop for different noise realizations
    for irng = 1:1000
        % set rnd number for reproduction
        rng(irng,'twister')

        % rand condition
        maxSteps = 1000;
        rM = 25; % radius state space
        cM = [0 0 25]; % center state space
        xRc = rM*(rand(maxSteps,3)-0.5)*2 + cM;

        % random initial conditions
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

        % data before smoothing for plotting
        xobsPlotE = xobs;

        % build library
        Theta_0 = build_theta(xobs,common_params);


        %% ENSEMBLES SINDY

        %% calculate derivatives
        dxobs_0 = dx;

        %% run active loop
        iLoop = 1;
        nSn = 1;
        nLoop = 81;
        % active = 0;
        % active = 1;

        while true

            if ~active
                %% SINDy
                % sindy with central difference differentiation
                sindy = sindy_cd(xobs,Theta_0,n,lambda,gamma,dt);

                % store outputs
                nWrongTermsS = sum(sum(abs((true_nz_weights~=0) - (sindy~=0))));
                modelErrorS = norm(sindy-true_nz_weights)/norm(true_nz_weights);
                successS = norm((true_nz_weights~=0) - (sindy~=0))==0;

                modelErrorSout(irng,iLoop) = modelErrorS;
                SRSout(irng,iLoop) = successS;
            end

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

            modelErrorDDEout(iLoop) = modelErrorDDE;
            XiDBsOutMax(iLoop) = max(max(XiDBs));
            XiDBsOutMean(iLoop) = mean(mean(XiDBs));
            XiDBsOutMean2(iLoop) = mean(mean(XiDBs(XiDBs~=0)));
            SRout(iLoop) = successDDE;

            tt(iLoop) = nS + (iLoop-1)*nSn;


            if active
            % plot variance prediction
            polysIN = 1:2; % skip last rows to oncrease speed, doesnt change results
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

            if iLoop == nLoop
                break
            end

            iLoop = iLoop + 1;

        end

        XiDBsOutMaxOut(irng,:) = XiDBsOutMax;
        XiDBsOutMeanOut(irng,:) = XiDBsOutMean;
        XiDBsOutMeanOut2(irng,:) = XiDBsOutMean2;
        SRoutOut(irng,:) = SRout;
        modelErrorDDEoutOut(irng,:) = modelErrorDDEout;

    end


    XiDBsOutMaxMean = mean(XiDBsOutMaxOut);
    XiDBsOutMeanMean = mean(XiDBsOutMeanOut);
    XiDBsOutMeanMean2 = mean(XiDBsOutMeanOut2);
    SRoutMean = mean(SRoutOut);
    mEmean = mean(modelErrorDDEoutOut);

    % run: 1000 noise instantiations, 200 points active, also save passive SINDy
    SRSoutMean = mean(SRSout);
    mSmean = mean(modelErrorSout);
    if active
        save('Results/activeOutAll')
        save('Results/activeOut','XiDBsOutMaxMean','XiDBsOutMeanMean','XiDBsOutMeanMean2','SRoutMean','mEmean','tt')
    else
        save('Results/passiveOutAll')
        save('Results/passiveOut','XiDBsOutMaxMean','XiDBsOutMeanMean','XiDBsOutMeanMean2','SRoutMean','mEmean','SRSoutMean','mSmean','tt')
    end

else

    activeOut = load('Results/activeOutAll');
    passiveOut = load('Results/passiveOut');
    tt = activeOut.tt;

end


%% final plot

warning('run active and passive')

fos = 14; % fontsize
fosL = 11; % fontsize legend
ttn = 1;
lw2 = 1.5;
lw1 = 0.75;
C1 = [0 0 0];
C2 = 0.5*[1 1 1];
sizeX = 250;
sizeY = 600;
figure('Position', [10 10 sizeX sizeY])
subplot(3,1,1)
plot(tt,activeOut.XiDBsOutMeanMean2,'Color',C1,'Linewidth',lw2); hold on
plot(tt,passiveOut.XiDBsOutMeanMean2,'--','Color',C2,'Linewidth',lw1); hold on
ylim([0 0.3])
xlim([min(tt)-ttn max(tt)+ttn])
yticks(0:0.05:0.25)
xticklabels([])
set(gca,'ticklabelinterpreter','latex','FontSize',fosL)
ylabel('$\bar{\sigma}$ coeff.','interpreter','latex','FontSize',fos)
legend({'active E-SINDy','E-SINDy'},'interpreter','latex','FontSize',fosL)
leg = legend({'active E-SINDy','E-SINDy'},'interpreter','latex','FontSize',fosL);
leg.ItemTokenSize = [10,12];
subplot(3,1,2)
plot(tt,activeOut.SRoutMean,'Color',C1,'Linewidth',lw2); hold on
plot(tt,passiveOut.SRoutMean,'--','Color',C2,'Linewidth',lw1); hold on
ylim([0 1])
xlim([min(tt)-ttn max(tt)+ttn])
yticks(0:0.25:1)
xticklabels([])
set(gca,'ticklabelinterpreter','latex','FontSize',fosL)
ylabel('success rate','interpreter','latex','FontSize',fos)
subplot(3,1,3)
plot(tt,activeOut.mEmean,'Color',C1,'Linewidth',lw2); hold on
plot(tt,passiveOut.mEmean,'--','Color',C2,'Linewidth',lw1); hold on
ylim([0 0.03])
xlim([min(tt)-ttn max(tt)+ttn])
yticks(0:0.01:0.03)
set(gca,'ticklabelinterpreter','latex','FontSize',fosL)
ylabel('model error','interpreter','latex','FontSize',fos)
xlabel('Nr. data samples','interpreter','latex','FontSize',fos)



