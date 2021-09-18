%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: generate 'tags' for terms to include in model,
%%%%%%%%%%%% indices for true model terms, true model weights 
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function [tags_pde,lib_list,pdx_list,lhs_ind,true_nz_weights] = get_lib_tags(n,dim,lhs,max_dx,max_dt,polys,trigs,use_cross_dx,use_all_dt,custom_remove,custom_add,true_nz_weight_tags)

if isempty(custom_remove)
    custom_remove = ones(1,n+dim)*Inf;
end

[~,lib_list,pdx_list] = build_fcn_lib_tags(n,dim,max_dx,max_dt,polys,trigs,use_cross_dx,use_all_dt);

if ~use_all_dt
    lib_list = lib_list(or(~lib_list(:,end)>0,ismember(lib_list,lhs,'rows')),:);
end
lib_list = unique([lib_list;custom_add],'rows');
lib_list = lib_list(~ismember(lib_list,custom_remove,'rows'),:);

[tags_pde,lib_list] = build_str_tags(lib_list,dim,n);

lhs_ind=zeros(size(lhs,1),1);
try
    for k=1:size(lhs,1)
        lhs_ind(k) = find(ismember(lib_list,lhs(k,:),'rows'));
    end
catch
    disp('ERROR: LHS not computed')
    return
end
if ~isempty(true_nz_weight_tags{1})
    num_eqs = size(lhs,1);
    try 
        true_nz_weights = zeros(size(lib_list,1)-num_eqs,num_eqs);
        for k=1:num_eqs
            axi_temp = tags2axi(true_nz_weight_tags{k},lib_list);
            true_nz_weights(:,k) = axi_temp(~ismember(1:size(lib_list,1),lhs_ind));
        end
    catch
        u_input = input('True terms missing from library, proceed anyway?'); 
        if u_input~=1 
            return;
        else
            true_nz_weights = [];
        end
    end
else
    true_nz_weights = [];
end
end
