
%% parameters reaction diffusion

% if exist('datasets/rxn_diff.mat')
%     load('datasets/rxn_diff.mat')
if exist('C:\Users\ufase\OneDrive - UW\Documents\GitHub\LargeData\rxn_diff.mat')
    load('C:\Users\ufase\OneDrive - UW\Documents\GitHub\LargeData\rxn_diff.mat')
else
    warning('Generate reaction diffusion data, might take a while...')
    reaction_diffusion
end

U_exact{1} = u;
U_exact{2} = v;

clear u v

xs{1} = x';
xs{2} = y';
xs{3} = t';

true_nz_weights{1} = [1 0 2 0 0 1/10;...
                      1 0 0 2 0 1/10;...
                      1 2 0 0 0 -1;...
                      3 0 0 0 0 -1;...
                      0 3 0 0 0 1;...
                      2 1 0 0 0 1;...
                      1 0 0 0 0 1];
                  
true_nz_weights{2} = [0 1 2 0 0 1/10;...
                      0 1 0 2 0 1/10;...
                      0 1 0 0 0 1;...
                      1 2 0 0 0 -1;...
                      3 0 0 0 0 -1;...
                      0 3 0 0 0 -1;...
                      2 1 0 0 0 -1];
                  
lhs = [1 0 0 0 1;...
       0 1 0 0 1;];
   
max_dx = 5;
max_dt = 1;
polys = 0:4;
trigs = [];
use_all_pt = 0;
use_cross_dx = 0;
custom_add = [];
custom_remove = [];

m_t = 14;
m_x = 13;

