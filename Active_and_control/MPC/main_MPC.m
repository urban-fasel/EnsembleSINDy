%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run E-SINDy MPC
% 
%

clear all
close all
clc

runSweep = 0; % run or load results

if runSweep
    %% Initialize problem
    dt        = 0.01; % Time step
    polyorder = 3;    % Library terms polynomial order
    usesine   = 0;    % Using sine functions in library

    % Lorenz system
    Beta = [10, 8/3, 28];   % system parameters
    x0 = [-8 7 27]';        % initial conditions
    eps = 0.01;             % Noise level
    n = length(x0);         % number of states
    tol_ode = 1e-10;        % set tolerance (abs and rel) of ode45
    options = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));
    [~,x10]=ode45(@(t,x) lorenz(t,x,Beta),dt:dt:10,x0,options); % trajectory for noise signal power calculation
    signal_power = rms(x10(:)); % signal power root mean square for noise 

    % true lorenz system terms for comparison with SINDy
    true_nz_weights = zeros(34,4);
    true_nz_weights(1,1) = -Beta(1);
    true_nz_weights(2,1) = Beta(1);
    true_nz_weights(4,1) = 1;
    true_nz_weights(1,2) = Beta(3);
    true_nz_weights(1,2) = Beta(3);
    true_nz_weights(2,2) = -1;
    true_nz_weights(7,2) = -1;
    true_nz_weights(3,3) = -Beta(2);
    true_nz_weights(6,3) = 1;


    for ii = 1:11

        DD = 10;
        Dt0 = 50;
        T_train = Dt0*dt + (ii-1)*DD*dt;

        tt(ii) = T_train;

        % Time span training
        tspanTrain = 0:dt:T_train; 

        % SPHS training signal
        Pf = 0.5; % Fundamental period
        K = 8; 
        A = 10;
        forcingSPHS = @(x,t) A*sphs(Pf,K,t);
        forcingTrain = forcingSPHS;
        % Define training forcing u(t)
        for i = 1:length(tspanTrain)
            uTrain(i,:) = forcingTrain(0,tspanTrain(i));
        end

        %% Generate Data
        % Integrate forced Lorenz system
        options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n)); % Define ode options
        [~,xTrain]=ode45(@(t,x) lorenzForcing(t,x,forcingTrain(x,t),Beta),tspanTrain,x0,options);
        xTrainNoNoise = xTrain;

        nNoise = 1000;
        nStart = 1;

        for irng = nStart:(nStart+nNoise)

            rng(irng,'twister')
            xTrain = xTrainNoNoise + signal_power*eps*randn(size(xTrainNoNoise));  % Add noise

            % calculate derivatives
            dx = zeros(size(xTrain));
            dx(1,:)=(-11/6*xTrain(1,:) + 3*xTrain(2,:) -3/2*xTrain(3,:) + xTrain(4,:)/3)/dt;
            dx(2:size(xTrain,1)-1,:) = (xTrain(3:end,:)-xTrain(1:end-2,:))/(2*dt);
            dx(size(xTrain,1),:) = (11/6*xTrain(end,:) - 3*xTrain(end-1,:) + 3/2*xTrain(end-2,:) - xTrain(end-3,:)/3)/dt;
            dx(:,n+1) = 0*dx(:,n);          % Add additional column of zeros for derivatives of u, needed in SINDy algorithm

            xTrain = [xTrain uTrain];       % Stack state and input measurements


            %% SINDy model Identification
            Theta = poolData(xTrain,(n+1),polyorder,usesine); % Generate library
            lambda = 0.2;     % lambda is our sparsification hyperparameter.
            Xi = sparsifyDynamics(Theta,dx,lambda,(n+1)); % Run SINDy using sequential treshold least squares
            if nNoise == 1
                poolDataLIST({'x','y','z','u'},Xi,(n+1),polyorder,usesine); % List library terms
            end


            %% ensemble SINDy
            [Xie, Xies, XiE] = ensemble(Theta,dx,lambda);


            %% run MPC
            % standard SINDy
            [xMPC, dxMPC, uMPC, J] = runMPCensemble(dt,Xi,x0,polyorder,usesine,Beta,n);
            % E-SINDy
            [xMPCe, dxMPCe, uMPCe, Je] = runMPCensemble(dt,Xie,x0,polyorder,usesine,Beta,n);

            uMPCout(irng,:) = uMPC;
            Jout(irng) = mean(J);
            xallout(irng,:) = xMPC(1,:);
            yallout(irng,:) = xMPC(2,:);
            zallout(irng,:) = xMPC(3,:);


            uMPCeout(irng,:) = uMPCe;
            Jeout(irng) = mean(Je);
            xeallout(irng,:) = xMPCe(1,:);
            yeallout(irng,:) = xMPCe(2,:);
            zeallout(irng,:) = xMPCe(3,:);

        end

        xallLout(:,:,ii) = xallout;
        yallLout(:,:,ii) = yallout;
        zallLout(:,:,ii) = zallout;
        uallLout(:,:,ii) = uMPCout;

        xeallLout(:,:,ii) = xeallout;
        yeallLout(:,:,ii) = yeallout;
        zeallLout(:,:,ii) = zeallout;
        ueallLout(:,:,ii) = uMPCeout;

        J_MPC(ii) = mean(Jout);
        Je_MPC(ii) = mean(Jeout);

        save('Results/ensembleMPCoutInter');

    end

    save('Results/ensembleMPCout2');

else
    
    % load data
    load('Results/ensembleMPCout');
    
end


%% MPC cost function
fos = 14;
fosS = 12;
sizeX = 400;
sizeY = 350;
figure('Position', [10 10 sizeX sizeY])
plot(tt/dt,Je_MPC,'b','Linewidth',1.2); hold on
plot(tt/dt,J_MPC,'r--','Linewidth',1.2); hold on
xticks(50:25:150)
xlim([48 152])
set(gca,'ticklabelinterpreter','latex','FontSize',fosS)
xlabel('Nr. time steps for training','interpreter','latex','FontSize',fos)
ylabel('$\bar{J}$: MPC cost function','interpreter','latex','FontSize',fos)
legend({'E-SINDy','SINDy'},'interpreter','latex','FontSize',fosS)



%% plot training and MPC trajectories -> long training
xref = [sqrt(72);sqrt(72);27];
xTrainPlot = xTrainNoNoise;
C1 = [1 1 1]*0.7;
C2 = [1 1 1]*0.2;
tTrain = 0:150;
tMPC = tTrain(end):tTrain(end)+400;
lw1 = 1.5;
lw2 = 1;
lw3 = 1.5;
lw4 = 0.5;
lw5 = 1;
fos = 14;
fosS = 10;

nNplot = 3;
sizeX = 400;
sizeY = 350;
figure('Position', [10 10 sizeX sizeY])
iP = 2;
nDLplot = [1 5];
subplot(4,1,1)
plot(tMPC,xeallLout(nNplot,:,nDLplot(iP)),'b','Linewidth',lw1); hold on
plot(tMPC,xallLout(nNplot,:,nDLplot(iP)),'r--','Linewidth',lw2); hold on
plot([tMPC(1) tMPC(end)],[xref(1) xref(1)],':','Color',C2,'Linewidth',lw5)
plot(tTrain,xTrainPlot(:,1)','Color',C1,'Linewidth',lw3); hold on
plot([tTrain(end) tTrain(end)],[-20 20],'k','Linewidth',lw4)
xlim([0 tMPC(end)])
xticks([])
set(gca,'ticklabelinterpreter','latex','FontSize',fosS)
ylabel('x','interpreter','latex','FontSize',fos)
title('150 time steps training','interpreter','latex','FontSize',fos)
leg=legend({'E-SINDy MPC','SINDy MPC','Reference','Training'},'interpreter','latex','FontSize',fosS,'NumColumns',2);%,'Location','Northeast');
leg.ItemTokenSize = [10,12];

subplot(4,1,2)
plot(tMPC,yeallLout(nNplot,:,nDLplot(iP)),'b','Linewidth',lw1); hold on
plot(tMPC,yallLout(nNplot,:,nDLplot(iP)),'r--','Linewidth',lw2); hold on
plot([tMPC(1) tMPC(end)],[xref(2) xref(2)],':','Color',C2,'Linewidth',lw5); hold on
plot(tTrain,xTrainPlot(:,2)','Color',C1,'Linewidth',lw3); hold on
plot([tTrain(end) tTrain(end)],[-20 20],'k','Linewidth',lw4)
xlim([0 tMPC(end)])
xticks([])
set(gca,'ticklabelinterpreter','latex','FontSize',fosS)
ylabel('y','interpreter','latex','FontSize',fos)

subplot(4,1,3)
plot(tMPC,zeallLout(nNplot,:,nDLplot(iP)),'b','Linewidth',lw1); hold on
plot(tMPC,zallLout(nNplot,:,nDLplot(iP)),'r--','Linewidth',lw2); hold on
plot([tMPC(1) tMPC(end)],[xref(3) xref(3)],':','Color',C2,'Linewidth',lw5); hold on
plot(tTrain,xTrainPlot(:,3)','Color',C1,'Linewidth',lw3); hold on
plot([tTrain(end) tTrain(end)],[0 50],'k','Linewidth',lw4)
xlim([0 tMPC(end)])
xticks([])
yticks(0:25:50)
set(gca,'ticklabelinterpreter','latex','FontSize',fosS)
ylabel('z','interpreter','latex','FontSize',fos)

subplot(4,1,4)
plot(tMPC,ueallLout(nNplot,:,nDLplot(iP)),'b','Linewidth',lw1); hold on
plot(tMPC,uallLout(nNplot,:,nDLplot(iP)),'r--','Linewidth',lw2); hold on
plot(tTrain,xTrainPlot(:,4)','Color',C1,'Linewidth',lw3); hold on
plot([tTrain(end) tTrain(end)],[-55 55],'k','Linewidth',lw4)
xlim([0 tMPC(end)])
set(gca,'ticklabelinterpreter','latex','FontSize',fosS)
xlabel('Time steps','interpreter','latex','FontSize',fos)
ylabel('u','interpreter','latex','FontSize',fos)

set(leg,...
    'Position',[0.395462443033854 0.742380952380953 0.482870890299479 0.104285714285714],...
    'NumColumns',2,...
    'Interpreter','latex',...
    'FontSize',10);



%% plot training and MPC trajectories -> short training
xref = [sqrt(72);sqrt(72);27];
xTrainPlot = xTrainNoNoise(1:51,:);
C1 = [1 1 1]*0.7;
C2 = [1 1 1]*0.2;
tTrain = 0:50;
tMPC = tTrain(end):tTrain(end)+400;
lw1 = 1.5;
lw2 = 1;
lw3 = 1.5;
lw4 = 0.5;
lw5 = 1;
fos = 14;
fosS = 10;

nNplot = 3;
sizeX = 400;
sizeY = 350;
figure('Position', [10 10 sizeX sizeY])
iP = 1;
nDLplot = [1 5];
subplot(4,1,1)
plot(tMPC,xeallLout(nNplot,:,nDLplot(iP)),'b','Linewidth',lw1); hold on
plot(tMPC,xallLout(nNplot,:,nDLplot(iP)),'r--','Linewidth',lw2); hold on
plot([tMPC(1) tMPC(end)],[xref(1) xref(1)],':','Color',C2,'Linewidth',lw5)
plot(tTrain,xTrainPlot(1:51,1)','Color',C1,'Linewidth',lw3); hold on
plot([tTrain(end) tTrain(end)],[-20 20],'k','Linewidth',lw4)
xlim([0 tMPC(end)])
xticks([])
set(gca,'ticklabelinterpreter','latex','FontSize',fosS)
ylabel('x','interpreter','latex','FontSize',fos)
title('50 time steps training','interpreter','latex','FontSize',fos)

subplot(4,1,2)
plot(tMPC,yeallLout(nNplot,:,nDLplot(iP)),'b','Linewidth',lw1); hold on
plot(tMPC,yallLout(nNplot,:,nDLplot(iP)),'r--','Linewidth',lw2); hold on
plot([tMPC(1) tMPC(end)],[xref(2) xref(2)],':','Color',C2,'Linewidth',lw5); hold on
plot(tTrain,xTrainPlot(1:51,2)','Color',C1,'Linewidth',lw3); hold on
plot([tTrain(end) tTrain(end)],[-20 20],'k','Linewidth',lw4)
xlim([0 tMPC(end)])
xticks([])
set(gca,'ticklabelinterpreter','latex','FontSize',fosS)
ylabel('y','interpreter','latex','FontSize',fos)

subplot(4,1,3)
plot(tMPC,zeallLout(nNplot,:,nDLplot(iP)),'b','Linewidth',lw1); hold on
plot(tMPC,zallLout(nNplot,:,nDLplot(iP)),'r--','Linewidth',lw2); hold on
plot([tMPC(1) tMPC(end)],[xref(3) xref(3)],':','Color',C2,'Linewidth',lw5); hold on
plot(tTrain,xTrainPlot(1:51,3)','Color',C1,'Linewidth',lw3); hold on
plot([tTrain(end) tTrain(end)],[0 50],'k','Linewidth',lw4)
xlim([0 tMPC(end)])
xticks([])
yticks(0:25:50)
set(gca,'ticklabelinterpreter','latex','FontSize',fosS)
ylabel('z','interpreter','latex','FontSize',fos)

subplot(4,1,4)
plot(tMPC,ueallLout(nNplot,:,nDLplot(iP)),'b','Linewidth',lw1); hold on
plot(tMPC,uallLout(nNplot,:,nDLplot(iP)),'r--','Linewidth',lw2); hold on
plot(tTrain,xTrainPlot(1:51,4)','Color',C1,'Linewidth',lw3); hold on
plot([tTrain(end) tTrain(end)],[-55 55],'k','Linewidth',lw4)
xlim([0 tMPC(end)])
set(gca,'ticklabelinterpreter','latex','FontSize',fosS)
xlabel('Time steps','interpreter','latex','FontSize',fos)
ylabel('u','interpreter','latex','FontSize',fos)





