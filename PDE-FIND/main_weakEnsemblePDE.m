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

clc
clear all
close all

% run weak-ensemble or standard weak pde-find
runEns = 1;

% choose PDE
pde_num = 2;
pde_names = {'burgers.mat','KdV.mat','KS.mat','NLS.mat','rxn_diff'};

% load data 
if pde_num ~= 5
    load(['datasets/',pde_names{pde_num}])
else
    RDparams
end

coarsen_data = [[1 1];[1 1];[1 1]];

% set noise level
if pde_num == 4 
    sigma_NR = 0.5; 
elseif pde_num == 5
    sigma_NR = 0.2;
else
    sigma_NR = 1.0;
end
noise_dist = 0; 
noise_alg = 0;


% Messenger et al run WSINDy on 200 instantiations of noise and average the results
% of error statistics: 200 runs sufficiently reduces variance in the results
nNoise = 1;%200;

for iN = 1:nNoise
    rng(iN,'twister')
    rng_seed = rng().Seed; 

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
        % Model error and success rate: weak-ensemble
        % nNoise = 200: 0.0252    1.0000
        % Model error and success rate: weak PDE-FIND
        % nNoise = 200: 0.0260    0.9900
    end

    if pde_num==2
        xIn = [400;0.20;0.95;0.95;100;0.02];
        % Model error and success rate: weak-ensemble
        % nNoise = 200: 0.1193    0.9950
        % Model error and success rate: weak PDE-FIND
        % nNoise = 200: 0.2745    0.9350 
    end

    if pde_num == 3
        xIn = [300 0.05 0.3 0.9 150 0.15];
        % Model error and success rate: weak-ensemble
        % nNoise = 200:  0.2473    0.9950
        % Model error and success rate: weak PDE-FIND
        % nNoise = 200: 0.0667    1.0000
    end

    if pde_num==4
        xIn = [300 0.05 0.2 0.9 150 0.1];
        % Model error and success rate: weak-ensemble
        % nNoise = 200: 0.1133    1.0000
        % Model error and success rate: weak PDE-FIND
        % nNoise = 200: 0.1299    0.8200
    end

    if pde_num == 5
        xIn = [200 0.05 0.3 0.95 200 0.05];
        % Model error and success rate: weak-ensemble
        % nNoise = 200: 0.0708    0.9950
        % Model error and success rate: weak PDE-FIND
        % nNoise = 200: 0.7766    0  
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
        [W,~,~,~,lambda_learned, ~,~, true_nz_weightsOut,~,dW,~,~,~,~,~,~,~, ~,~,~,~,~,~,~,~,~,~,XiDBs] = ... 
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
    error(iN) = norm(W-true_nz_weightsOut)/norm(true_nz_weightsOut);
    success(iN) = norm((true_nz_weightsOut~=0) - (W~=0))==0;

    % % error calculated such as in Sam Rudy PDE-FIND 2017
    % errorSamRudy(iN) = mean(abs(true_nz_weightsOut(true_nz_weightsOut~=0) - W(true_nz_weightsOut~=0))./abs(true_nz_weightsOut(true_nz_weightsOut~=0)));
    % stdSamRudy(iN) = std(abs(true_nz_weightsOut(true_nz_weightsOut~=0) - W(true_nz_weightsOut~=0))./abs(true_nz_weightsOut(true_nz_weightsOut~=0)));

end

% model error and success rate for comparison
out = [mean(error) mean(success)] %mean(errorSamRudy) mean(stdSamRudy)]


%% plot 
if nNoise == 1
    print_loc = 1;
    toggle_plot_basis_fcn = 1; 
    toggle_plot_sol = 1; 
    toggle_plot_loss = 1; 

    if print_loc~=0
        print_results(W,G,resid,dW,print_loc,dims,polys,trigs,max_dx,max_dt,lambda_hat,gamma,lhs_ind,tags_pde,m_x,m_t,p_x,p_t,s_x,s_t,scales,ET_wsindy,its_all)
    end

    if toggle_plot_basis_fcn
        figure(1);
        plot_basis_fcn(Cfs_x,Cfs_t,m_x,dx,m_t,dt,max_dx,max_dt,pdx_list,[],scales(end-nstates:end));
    end

    if and(toggle_plot_sol,dim==2)
        figure(2);clf;set(gcf,'units','normalized','outerposition',[0.4 0 0.3 0.45])
        surf(xs{1},xs{2},U_obs{1}', 'EdgeColor','none')
        view([15 55])
        xlabel('$x$','interpreter','latex','fontsize',14)
        ylabel('$t$','interpreter','latex','fontsize',14)
        set(gca, 'TickLabelInterpreter','latex','fontsize',14)
    end

    if length(lambda)>1 && toggle_plot_loss
        figure(3);clf;set(gcf,'units','normalized','outerposition',[0.7 0 0.3 0.45])
        loglog(thrs_EL(2,:),thrs_EL(1,:),'o-')
        xlabel('$\lambda$','interpreter','latex','fontsize',14)
        ylabel('$\mathcal{L}$','interpreter','latex','fontsize',14)
        set(gca, 'TickLabelInterpreter','latex','fontsize',14)
    end
end