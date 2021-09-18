%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: generate true model weights 
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function true_nz_weights = tags2axi(true_nz_weight_tags,lib_list)
  m = size(lib_list,1);
  true_nz_weights = zeros(m,1);
  [~,loc] = ismember(true_nz_weight_tags(:,1:end-1),lib_list,'rows');
  true_nz_weights(loc) = true_nz_weight_tags(:,end);
end