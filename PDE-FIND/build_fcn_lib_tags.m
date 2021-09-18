%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: generate 'tags' for terms to include in model
%%%%%%%%%%%% library. For instance, for solution data (u,v) over (x,y,t), 
%%%%%%%%%%%% the term (u^p*v^q)_{xxt} is denoted [p q 2 0 1] 
%%%%%%%%%%%% the first two entries give the powers of u and v 
%%%%%%%%%%%% the last entries give the partial derivative orders w/r/t (x,y,t) 
%%%%%%%%%%%% For trig terms such as (sin(ku))_{xyt}, use 
%%%%%%%%%%%% [-1ik 0  1 1 1], and for (cos(ku))_{xyt}, use
%%%%%%%%%%%% [  0 1ik 1 1 1], and for (cos(kv))_{xyt}, use
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function [tags,lib_list,pdx_list] = build_fcn_lib_tags(n,dim,max_deriv_x,max_deriv_t,polys,trigs,use_cross_derivs,use_all_dt)
tags = [];
for p = 1:length(polys)
    tags = [tags;partitionNk(polys(p),n)];
end
for k=1:length(trigs)
    trig_inds = [-trigs(k)*1i*eye(n);trigs(k)*1i*eye(n)];
    tags = [tags; trig_inds];
end

if and(use_cross_derivs,use_all_dt)
    inds_grid = cell(dim,1);
    for d=1:dim-1
        inds_grid{d} = 1:max_deriv_x+1;
    end
    inds_grid{dim} = 1:max_deriv_t+1;
    [inds_grid{:}] = ndgrid(inds_grid{:});

    pdxs = (0:max(max_deriv_x,max_deriv_t))';
    pdx_list = [];
    for d=1:dim-1
        pdx_list = [pdx_list pdxs(inds_grid{d},:)];
    end
    pdx_list = [pdx_list pdxs(inds_grid{end},:)];
elseif and(use_cross_derivs,~use_all_dt)
    inds_grid = cell(dim-1,1);
    for d=1:dim-1
        inds_grid{d} = 1:max_deriv_x;
    end
    [inds_grid{:}] = ndgrid(inds_grid{:});

    pdxs = (0:max_deriv_x)';
    pdx_list = [];
    for d=1:dim-1
        pdx_list = [pdx_list pdxs(inds_grid{d},:)];
    end
    pdx_list = [zeros(1,dim);[pdx_list [max_deriv_t;zeros(size(pdx_list,1)-1,1)]]];
elseif and(~use_cross_derivs,use_all_dt)    
    pdx_list = zeros(1,dim);
    for k=1:max_deriv_t                                                        
        pdx_list = [pdx_list;[zeros(1,dim-1) k]];
    end
    for k=1:dim-1
        for j=1:max_deriv_x
            pdx_list = [pdx_list;[zeros(1,k-1) j zeros(1,dim-k)]];
        end
    end
elseif and(~use_cross_derivs,~use_all_dt)
    pdx_list = zeros(1,dim);
    pdx_list = [pdx_list;[zeros(1,dim-1) max_deriv_t]];
    for k=1:dim-1
        for j=1:max_deriv_x
            pdx_list = [pdx_list;[zeros(1,k-1) j zeros(1,dim-k)]];
        end
    end
end

if any(ismember(polys,0))
    lib_list = zeros((size(tags,1)-1)*size(pdx_list,1),size(tags,2)+size(pdx_list,2));    
    for i=2:size(tags,1)
        for j=1:size(pdx_list,1)
            lib_list(1+(i-2)*size(pdx_list,1)+j,:) = [tags(i,:) pdx_list(j,:)];
        end
    end
else
    [a,b] = ndgrid(1:size(tags,1),1:size(pdx_list,1));
    lib_list = [tags(a,:) pdx_list(b,:)];
end
end
