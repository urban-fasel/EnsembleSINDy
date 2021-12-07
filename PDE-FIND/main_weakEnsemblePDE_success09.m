%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%% Weak-ENSEMBLE PDE-FIND: 
%%%%%%%%%%%%  comparing standard weak PDE-find with weak-ensemble PDE-find
%%%%%%%%%%%% 
%%%%%%%%%%%%
%%%%%%%%%%%% Code by Messenger et al 2020, modified by Urban Fasel 2021
%%%%%%%%%%%%
%%%%%%%%%%%% pde_num selects a PDE system from the list pde_names
%%%%%%%%%%%% noise_ratio sets the signal-to-noise ratio (L2 sense)
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz
%%%%%%%%%%%%
%%%%%%%%%%%% modified by Urban Fasel
%%%%%%%%%%%%
%%%%%%%%%%%% increse noise up to point where success rate drops below 0.9
%%%%%%%%%%%%

clc
clear all
close all

% run weak-ensemble or standard weak pde-find
runEns = 1;

% choose PDE
pde_num = 5;
pde_names = {'burgers.mat','KdV.mat','KS.mat','NLS.mat','rxn_diff'};

% load data 
if pde_num ~= 5
    load(['datasets/',pde_names{pde_num}])
else
    RDparams
end

coarsen_data = [[1 1];[1 1];[1 1]];

% set noise level
% if pde_num == 4 
%     sigma_NR = 0.5; 
% elseif pde_num == 5
%     sigma_NR = 0.2;
% else
%     sigma_NR = 1.0;
% end
noise_dist = 0; 
noise_alg = 0;

nstates = length(U_exact);
dim = length(size(U_exact{1}));
inds = cell(1,dim);
for j=1:dim
    N = length(xs{j});
    inds{j} = 1:coarsen_data(j,1):floor(N/coarsen_data(j,2));
    xs{j} = xs{j}(inds{j});
end
for j=1:nstates
    U_exact{j} = U_exact{j}(inds{:});
end
    
    
% Messenger et al run WSINDy on 200 instantiations of noise and average the results
% of error statistics: 200 runs sufficiently reduces variance in the results
nNoise = 199;
iStart = 1;

errorMax = 0.05;
errorII = 0;
successII = 1;
successMin = 0.90;
sigma_NR0 = 0.3;
sigma_NR = sigma_NR0;
deltaNoise = 0.01;
nS = 5;
iii = 1;


%% results for nNoise=199
% KS 0.90: weak=0.98 (nNoise=199, started 0.9), e-weak=1.59 (nNoise=199, started at 1.4) 

% burgers 0.90: weak=2.09 (nNoise=999, started 1.8), e-weak=3.35 (nNoise=999, started at 3.3)

% KdV 0.90: weak=1.17(nNoise=199, started 1.1), e-weak=2.31 (nNoise=199, started at 1.5) 

% NLS 0.90: weak=0.43(nNoise=199, started 0.4), e-weak=0.93 (nNoise=199, started at 0.9) 

% rxn_diff 0.90: weak=0.08 (was 0.06) (nNoise=199, started 0.06), e-weak=0.32 (nNoise=199, started at 0.3) 

% E-WSINDy, WSINDy, PDE-FIND
sigma_NR_List = [3.35 2.09 0.01;...
    2.31 1.17 0.01;...
    1.59 0.98 0.01;...
    0.93 0.43 0.01;...
    0.32 0.08 0.005];
meanNR = mean(sigma_NR_List(:,1)./sigma_NR_List(:,2));


while successII >= successMin

    
for iN = iStart:iStart+nNoise
    rng(iN,'twister')
    rng_seed = rng().Seed; 

    rng(rng_seed);
    [U_obs,noise,snr,sigma] = gen_noise(U_exact,sigma_NR,noise_dist,noise_alg,rng_seed,0);
    dims = size(U_obs{1});
    dim = length(dims);

    %% Set hyperparameters

    %---------------- weak discretization
    if pde_num ~= 5
        s_x = floor(length(xs{1})/50);
        s_t = floor(length(xs{end})/50);
    else
        s_x = 13;
        s_t = 12;
    end

    phi_class = 1;
    tau = 10^-10;
    tauhat = 0; % p.17: B, KdV, KS: tauhat = 3, NLS, SG, RD: tauhat = 1, KS and NLS choose mx and mt 
    if pde_num == 5
        tauhat = 1;
    end
    toggle_scale = 2;
%     toggle_scale = 0;

    %---------------- regularized least squares solve

    % lambda chosen using information criteria (fixed for weak-ensemble)
    lambda = 10.^(linspace(-4,0,50));
    gamma = 0;

    %---------------- find test function hyperparameters from tau,tauhat: only for rxn_diff
    if tauhat >0
        [m_x,m_t,p_x,p_t,sig_est,corners_all] = findcorners(U_obs,xs,tau,tauhat,max_dx,max_dt,phi_class);
        tols = [-p_x -p_t];
    else
        tols = [tau tau];
    end

    %% ensemble hyperparameters for PDEs
    if pde_num==1
        xIn = [400 0.05 0.6 0.80 400 0.25];
    end

    if pde_num==2
        xIn = [400;0.4;0.95;0.95;100;0.02];
    end

    if pde_num == 3
        xIn = [300 0.05 0.3 0.9 150 0.15];
    end

    if pde_num==4
        xIn = [300 0.05 0.2 0.9 150 0.1];
    end

    if pde_num == 5
        xIn = [200 0.05 0.3 0.95 200 0.05];
    end


    %---------------- library bagging hyperparameters
    if runEns
        LBp.LB = 1; % run library bagging ensemble method
        LBp.DB = 1; % run double bagging -> bootstraping the data after library bagging

        LBp.nE1 = xIn(3); % use 60% of the library terms in each model
        LBp.nE2 = round(xIn(1)); % use 100 models in ensemble 
        LBp.ensT = xIn(2); % threshold parameters below this inclusion probability
        LBp.nE3 = round(xIn(5));
        LBp.ensT2 = xIn(4);
        lambda = xIn(6); % single lambda, manually tuned
    else
        LBp.LB = 0; % run library bagging ensemble method
        LBp.DB = 0; % run double bagging -> bootstraping the data after library bagging
    end


    %% Find dynamics

    if nNoise > 1
        [W,~,~,~,lambda_learned, ~,~, true_nz_weightsOut,resid,dW,~,~,~,~,~,~,~, ~,~,~,~,~,~,~,~,~,~,XiDBs] = ... 
            wsindy_pde_fun(U_obs,xs,lambda,gamma,true_nz_weights,lhs,max_dx,max_dt,polys,trigs,custom_add,custom_remove,use_all_pt,use_cross_dx,...
                            toggle_scale,m_x,m_t,s_x,s_t,tols,phi_class,LBp);
    else
        [W,G,b,M,lambda_hat, tags_pde_G,lib_list_G, true_nz_weightsOut,resid,dW,its_all,thrs_EL,ET_wsindy,tags_pde,lib_list,pdx_list,lhs_ind, Cfs_x,Cfs_t,dx,dt,p_x,p_t,sub_inds,scales,M_full,Theta_pdx,XiDBs] = ... 
            wsindy_pde_fun(U_obs,xs,lambda,gamma,true_nz_weights,lhs,max_dx,max_dt,polys,trigs,custom_add,custom_remove,use_all_pt,use_cross_dx,...
                            toggle_scale,m_x,m_t,s_x,s_t,tols,phi_class,LBp);    
        lambda_learned = lambda_hat;
    end

    %% results

    % model error
    error(iN-iStart+1) = norm(W-true_nz_weightsOut)/norm(true_nz_weightsOut);
    success(iN-iStart+1) = norm((true_nz_weightsOut~=0) - (W~=0))==0;
    
    % error calculated such as in Sam Rudy PDE-FIND 2017
    errorSamRudy(iN-iStart+1) = mean(abs(true_nz_weightsOut(true_nz_weightsOut~=0) - W(true_nz_weightsOut~=0))./abs(true_nz_weightsOut(true_nz_weightsOut~=0)));
    stdSamRudy(iN-iStart+1) = std(abs(true_nz_weightsOut(true_nz_weightsOut~=0) - W(true_nz_weightsOut~=0))./abs(true_nz_weightsOut(true_nz_weightsOut~=0)));

end

errorIIall(iii,:) = error;
successIIall(iii,:) = success;

error(error>1)=1; % max error

errorII(iii) = mean(error);
errorIISR(iii) = mean(errorSamRudy);
successII(iii) = mean(success)
sigma_NRII(iii) = sigma_NR

sigma_NR = sigma_NR + deltaNoise;

iii = iii+1;

end

% model error and success rate for comparison
out = [mean(error) mean(success) mean(errorSamRudy) mean(stdSamRudy)]


figure
plot(errorII)
figure
plot(errorIISR)
% 
figure
plot(successII)

for i = 1:nS
    errorIIalln = errorIIall(i,:);
    errorIIalln(errorIIalln>20*median(errorIIalln)) = [];
    sEE(i) = length(errorIIalln);
    errorIIallOut(i) = mean(errorIIalln);
end

figure
plot(sigma_NR0:deltaNoise:sigma_NR0+(nS-1)*deltaNoise,errorIIallOut)

