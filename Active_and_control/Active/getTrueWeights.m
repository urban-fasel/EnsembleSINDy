%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy: function for building data matrix Theta_0
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy: Galerkin-based Data-Driven Model
%%%%%%%%%%%% Selection"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz
%%%%%%%%%%%%
%%%%%%%%%%%% Modified by Urban Fasel, 2021 


function true_nz_weights = getTrueWeights(ode_params,common_params,n)

polys = common_params{1}; 
trigs = common_params{2}; 

sigma = ode_params{1}; 
beta = ode_params{2}; 
rho = ode_params{3};
weights = {[[0 1 0 sigma];[1 0 0 -sigma]],...
                      [[1 0 0 rho];[1 0 1 -1];[0 1 0 -1]],...
                      [[1 1 0 1];[0 0 1 -beta]]}; 
                  
tags = getTags(n,polys,trigs);
true_nz_weights = get_true_weights(weights,tags,n);

end

function tags = getTags(n,polys,trigs)

% n = size(xobs,2);

tags = [];
for p = 1:length(polys)
    monom_powers = partitionNk(polys(p),n);
    tags = [tags;monom_powers];
end
for k=1:length(trigs)
    trig_inds = [-trigs(k)*1i*eye(n);trigs(k)*1i*eye(n)];
    tags = [tags; trig_inds];
end
end


function total_parts = partitionNk(N,k)
    if N>1
        A = partitions(N,0:N,N,k);
        num_unique_parts = size(A,1);
        unique_parts = zeros(num_unique_parts,k);
        for j=1:num_unique_parts
            temp_part = [];
            for m= 1:N+1
                if A(j,m)>0
                    for l = 1:A(j,m)
                        temp_part = [temp_part m-1];
                    end
                end
            end
            unique_parts(j,:) = temp_part;
        end
        total_parts = [];
        for j=1:num_unique_parts
            total_parts = [total_parts;unique(perms(unique_parts(j,:)),'rows')];
        end
    elseif N==1
        total_parts = eye(k);
    elseif N==0
        total_parts = zeros(1,k);
    end
end

function true_nz_weights = get_true_weights(weights,tags,n)
    true_nz_weights = zeros(size(tags,1),n);
    for i = 1:length(weights)
        weights_i = weights{i};
        [l1,l2] = size(weights_i);
        for j = 1:l1
            true_nz_weights(all(weights_i(j,1:l2-1) == tags,2),i) = weights_i(j,l2);
        end
    end
end