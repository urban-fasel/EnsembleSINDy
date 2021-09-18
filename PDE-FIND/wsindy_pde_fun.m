%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: script for recoverying PDE systems
%%%%%%%%%%%% 
%%%%%%%%%%%% pde_num selects a PDE system from the list pde_names
%%%%%%%%%%%% noise_ratio sets the signal-to-noise ratio (L2 sense)
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function [W,G,b,M,lambda_learned,tags_pde_G,lib_list_G,true_nz_weights,resid,dW,its_all,thrs_EL,ET_wsindy,tags_pde,lib_list,pdx_list,lhs_ind, Cfs_x,Cfs_t,dx,dt,p_x,p_t,sub_inds,scales,M_full,Theta_pdx,XiDBs] = wsindy_pde_fun(U_obs,xs,lambda,gamma,true_nz_weight_tags,...
    lhs,max_dx,max_dt,polys,trigs,custom_add,custom_remove,use_all_pt,use_cross_dx,...
    toggle_scale,supp_phi_x,supp_phi_t,s_x,s_t,tols,phi_class,LBp)

    
n = length(U_obs);
dims = size(squeeze(U_obs{1}));
dim = length(dims);

if toggle_scale > 0
    scale_u = zeros(1,n);
    if toggle_scale == 1
        scale_u_fcn = @(v) norm(v(:)/norm(v(:),1)^(1/max(polys)),max(polys))^(max(polys)/(max(polys)-1));
    elseif toggle_scale == 2
        scale_u_fcn = @(v) norm(v(:)/norm(v(:),2)^(1/max(polys)),2*max(polys))^(max(polys)/(max(polys)-1));
    elseif toggle_scale == Inf
        scale_u_fcn = @(v) norm(v,inf)/(10^(1/max(polys)));
    end
    scale_u(1) = scale_u_fcn(U_obs{1}(:));
    for k=2:n 
        scale_u(k) = scale_u_fcn(U_obs{k}(:));
    end
else
    scale_u = [];
end

tic,
[tags_pde,lib_list,pdx_list,lhs_ind,true_nz_weights] = get_lib_tags(n,dim,lhs,max_dx,max_dt,polys,trigs,use_cross_dx,use_all_pt,custom_remove,custom_add,true_nz_weight_tags);
[Cfs_x,Cfs_t,dx,dt,p_x,p_t,sub_inds] = get_testfcn_weights(dims,xs,max_dx,max_dt,supp_phi_x,supp_phi_t,s_x,s_t,tols,phi_class);
if isempty(scale_u)
    scales = [];
    M_full = [];
else
    [scales,M_full] = get_scales(scale_u,p_x,supp_phi_x,dx,max_dx,p_t,supp_phi_t,dt,max_dt,lib_list,dim,phi_class,Cfs_x,Cfs_t);
end
Theta_pdx = get_lib_columns(n,lib_list,U_obs,Cfs_x,Cfs_t,supp_phi_x,supp_phi_t,dx,dt,sub_inds,dim,scales);
if length(lambda)>1
    [W,G,b,resid,dW,its_all,thrs_EL,M] = wsindy_pde_RGLS_seq(lambda,gamma,Theta_pdx,lhs_ind,true_nz_weights,M_full,LBp);
    XiDBs = [];
else
    [W,G,b,resid,dW,its_all,thrs_EL,M,XiDBs] = wsindy_pde_RGLS(lambda,gamma,Theta_pdx,lhs_ind,true_nz_weights,M_full,LBp);
end
ET_wsindy = toc;

tags_pde_G = tags_pde{~ismember(1:size(W,1),lhs_ind)};
lib_list_G = lib_list(~ismember(1:size(W,1),lhs_ind),:);
lambda_learned = thrs_EL(min(end,4),thrs_EL(min(end,4),:)>0);
lambda_learned = lambda_learned(end);

end