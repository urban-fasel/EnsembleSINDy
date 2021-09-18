%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: print results of WSINDy_PDE to 'filename', with
%%%%%%%%%%%% filename = [] leading to results printed in command window
%%%%%%%%%%%%  
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function print_results(W,G,resid,dW,filename,dims,polys,trigs,max_dx,max_dt,lambda_learned,gamma,lhs_ind,tags_pde,supp_phi_x,supp_phi_t,p_x,p_t,s_x,s_t,scales,ET_wsindy,its_all)

[m,n] = size(W);

if ~isequal(filename,1)    
    filename = fopen(filename,'a');
end

for k=1:n
    tags_pde_rdx = tags_pde(~ismember(1:m+n,lhs_ind));
    str_wsindy = print_pde(W(:,k),tags_pde_rdx,tags_pde{lhs_ind(k)});
    fprintf(filename,['\nRecovered PDE: ',str_wsindy]);
    fprintf(filename,'\nRelative Res: ||b-G*W||_2/||b||_2 = %.2e',norm(resid(:,k)));
    if ~isempty(dW)
        fprintf(filename,'\nMax Weight Error: max|W-W_{true}| = %.2e\n', dW{1}(k));
    end
end

% fprintf(filename,'      \n');
% fprintf(filename,'polys = ');
% fprintf(filename,'%u ',polys);
% fprintf(filename,'\ntrigs = ');
% fprintf(filename,'%u ',trigs);
% fprintf(filename,'\nMax derivs [t x] = ');
% fprintf(filename,'%u ',[max_dt max_dx]);
% fprintf(filename,'\n[m_x m_t] = ');
% fprintf(filename,'%u ',[supp_phi_x supp_phi_t]);
% fprintf(filename,'\n[s_x s_t] = ');
% fprintf(filename,'%u ',[s_x s_t]);
% fprintf(filename,'\n[p_x p_t] = ');
% fprintf(filename,'%u ',[p_x p_t]) ;
% fprintf(filename,'\n scales = ');
% fprintf(filename,'%u ',scales) ;
% 
% fprintf(filename,'\n      \n');
% fprintf(filename,'Size of dataset = ');
% fprintf(filename,'%u ',dims);
% fprintf(filename,'\nSize G = ');
% fprintf(filename,'%u ',size(G));
% if gamma >0
%     fprintf(filename,'\nCond G = %.2e',cond([G;gamma*eye(m)]));
% else
%     fprintf(filename,'\nCond G = %.2e',cond(G));
% end
% fprintf(filename,'\n[lambda_hat gamma] = ');
% fprintf(filename,'%.3e ',[lambda_learned gamma]);
% 
% 
% fprintf(filename,'\n      \n');
% fprintf(filename,'Elapsed time WSINDy = %4.4f \n',ET_wsindy);
% fprintf(filename,'STLS its = ');
% fprintf(filename,'%u ',its_all);
% fprintf(filename,'\n ');
% 
% if ~all(filename==1)
%     fclose(filename);
% end

end