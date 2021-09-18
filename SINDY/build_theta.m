%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy: function for building data matrix Theta_0
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy: Galerkin-based Data-Driven Model
%%%%%%%%%%%% Selection"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function Theta_0 = build_theta(xobs,common_params)

%%%%%%%%%%%%%%%%% get params

polys = common_params{1}; 
trigs = common_params{2}; 

%%%%%%%%%%%%%%%%% build Theta_0 & true_nz_weights

[Theta_0,~] = poolDataGen(xobs,polys,trigs);

end

function [Theta_0,tags] = poolDataGen(xobs,polys,trigs)

n = size(xobs,2);
ind = 0;

tags = [];
for p = 1:length(polys)
    monom_powers = partitionNk(polys(p),n);
    num_monoms = size(monom_powers,1);
    for j = 1:num_monoms
        Theta_0(:,ind+1) = prod(xobs.^(monom_powers(j,:)),2);
        ind = ind+1;
    end
    tags = [tags;monom_powers];
end
for k=1:length(trigs)
    trig_inds = [-trigs(k)*1i*eye(n);trigs(k)*1i*eye(n)];
    Theta_0 = [Theta_0 sin(trigs(k)*xobs) cos(trigs(k)*xobs)];
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
