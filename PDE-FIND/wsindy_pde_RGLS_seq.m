%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: solve regularized LSP with sequential thresh.
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function [W,G,b,resid,dW,its_all,rnorms,M] = wsindy_pde_RGLS_seq(lambdas,gamma,Theta_pdx,lhs_ind,axi,M_scale,LBp)

num_eq = length(lhs_ind);
[K,m] = size(Theta_pdx);
    
G = Theta_pdx(:,~ismember(1:m,lhs_ind));
b = zeros(K,num_eq);
for k=1:num_eq
    b(:,k) = Theta_pdx(:,lhs_ind(k));
end
%--------------------------------------------------------------
W_ls = G \ b;
GW_ls = norm(G*W_ls);

proj_cost = [];
overfit_cost = [];
rnorms = [];

alpha = 1/2; 
W = zeros(m-num_eq,num_eq);
for l=1:length(lambdas)
    lambda = lambdas(l);
    M = [];
    for k=1:num_eq
        if isempty(M_scale)
            M = [M M_scale(~ismember(1:m,lhs_ind))/M_scale(lhs_ind(k))];
%             [W(:,k),~,~] = sparsifyDynamics(G, b(:,k), lambda, gamma, M(:,end));
            if LBp.LB
                [W(:,k),~,~] = sparsifyDynamicsLB(G, b(:,k), lambda, gamma, M(:,end), LBp);    
            else
                [W(:,k),~,~] = sparsifyDynamics(G, b(:,k), lambda, gamma, M(:,end)); 
            end
        else
            M = M_scale(~ismember(1:m,lhs_ind))/M_scale(lhs_ind(k));
%             [W(:,k),~,~] = sparsifyDynamics(G, b(:,k), lambda, gamma, M);
            if LBp.LB
                [W(:,k),~,~] = sparsifyDynamicsLB(G, b(:,k), lambda, gamma, M, LBp);    
            else
                [W(:,k),~,~] = sparsifyDynamics(G, b(:,k), lambda, gamma, M); 
            end
            W(:,k) = W(:,k)./M;
        end
    end
    
    proj_cost = [proj_cost 2*alpha*norm(G*(W-W_ls))/GW_ls];
    overfit_cost = [overfit_cost 2*(1-alpha)*length(find(W~=0))/length(W(:))];
    rnorms = [rnorms proj_cost(end) + overfit_cost(end)];
end

l = find(rnorms == min(rnorms),1);

lambda = lambdas(l);
% lambda = 0.223;
its_all = zeros(num_eq,1);
if ~isempty(axi)
    dW = cell(num_eq+1,1);
else
    dW = [];
end

M = [];
for k=1:num_eq
    if ~isempty(M_scale)
        M = [M M_scale(~ismember(1:m,lhs_ind))/M_scale(lhs_ind(k))];
        if LBp.LB
            [W(:,k),its,~] = sparsifyDynamicsLB(G, b(:,k), lambda, gamma, M(:,end), LBp);    
        else
            [W(:,k),its,~] = sparsifyDynamics(G, b(:,k), lambda, gamma, M(:,end)); 
        end
    else
        if LBp.LB
            [W(:,k),its,~] = sparsifyDynamicsLB(G, b(:,k), lambda, gamma, M, LBp);
        else
            [W(:,k),its,~] = sparsifyDynamics(G, b(:,k), lambda, gamma, M);
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
rnorms = [rnorms;lambdas; [[rnorms(1:l);lambdas(1:l)] zeros(2,length(lambdas)-l)]; proj_cost; overfit_cost];
end
