%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: solve regularized LSP with sequential thresh.
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function [W,G,b,resid,dW,its_all,thrs_EL,M,XiDBs] = wsindy_pde_RGLS(lambda,gamma,Theta_pdx,lhs_ind,axi,M_scale,LBp)

num_eq = length(lhs_ind);
[K,m] = size(Theta_pdx);
    
b = zeros(K,num_eq);
W = zeros(m-num_eq,num_eq);
if ~isempty(axi)
    dW = cell(num_eq+1,1);
else
    dW = [];
end
G = Theta_pdx(:,~ismember(1:m,lhs_ind));
its_all = zeros(num_eq,1);


M = [];
for k=1:num_eq
    b(:,k) = Theta_pdx(:,lhs_ind(k));
    if isempty(M_scale)
%         [W(:,k),its,thrs_EL] = sparsifyDynamics(G, b(:,k), lambda, gamma, M);
        if LBp.LB
            [W(:,k),its,thrs_EL,XiDBs(:,k)] = sparsifyDynamicsLB(G, b(:,k), lambda, gamma, M, LBp);    
        else
            [W(:,k),its,thrs_EL] = sparsifyDynamics(G, b(:,k), lambda, gamma, M); 
            XiDBs = [];
        end
    else
        M = [M M_scale(~ismember(1:m,lhs_ind))/M_scale(lhs_ind(k))];
%         [W(:,k),its,thrs_EL] = sparsifyDynamics(G, b(:,k), lambda, gamma, M(:,end));
        if LBp.LB
            [W(:,k),its,thrs_EL,XiDBs(:,k)] = sparsifyDynamicsLB(G, b(:,k), lambda, gamma, M(:,end), LBp);    
        else
            [W(:,k),its,thrs_EL] = sparsifyDynamics(G, b(:,k), lambda, gamma, M(:,end)); 
            XiDBs = [];
        end
    end
    if ~isempty(axi)
        dW{k+1} = W(:,k)-axi(:,k);
        dW{1}(k) = norm(dW{k+1},inf)/norm(axi(:,k),inf);
    end
    its_all(k) = its;
end

if ~isempty(M_scale)
    resid = ((b./M_scale(lhs_ind)') - (G./M_scale(~ismember(1:m,lhs_ind))')*W)/norm(b./M_scale(lhs_ind)');
else
    resid = (b - G*W)/norm(b);
end


end
