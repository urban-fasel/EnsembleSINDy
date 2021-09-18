%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: get test function values on scaled reference grid,
%%%%%%%%%%%% generate indices for query points (sub_inds)
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function [Cfs_x,Cfs_t,dx,dt,p_x,p_t,sub_inds] = get_testfcn_weights(dims,xs,max_dx,max_dt,supp_phi_x,supp_phi_t,s_x,s_t,tols,phi_class)

dim = length(dims);

if or(2*supp_phi_x+1>min(dims(1:end-1)),2*supp_phi_t+1>min(dims(end)))
    disp('ERROR: Test function not compactly supported')
    return
end

sub_inds = cell(1,dim);
for j=1:dim-1
    sub_inds{j} = 1:s_x:dims(j)-2*supp_phi_x;
end
sub_inds{dim} = 1:s_t:dims(dim)-2*supp_phi_t;

dx = xs{1}(2)-xs{1}(1);
dt = xs{end}(2)-xs{end}(1);
[Cfs_x,p_x] = phi_int_weights(supp_phi_x,max_dx,tols(1),phi_class);
[Cfs_t,p_t] = phi_int_weights(supp_phi_t,max_dt,tols(2),phi_class);

end
