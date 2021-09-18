%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: generate string for identified PDE
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function str = print_pde(w,tags_pde,lhs)
    nz_inds = find(w);
    if ~isempty(nz_inds)
        str = [lhs,' = ',num2str(w(nz_inds(1)),6), tags_pde{nz_inds(1)}];
        for k=2:length(nz_inds)
            str = [str,' + ',num2str(w(nz_inds(k)),6),tags_pde{nz_inds(k)}];
        end
    else
        str = '';
    end
end
