%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: generate strings for library terms
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function [tags_pde,lib_list] = build_str_tags(lib_list,dim,n)
    tags_pde = {};
    ind = 1;
    remove_inds = zeros(size(lib_list,1),1);
    for j=1:size(lib_list)
        tags = lib_list(j,:);
            if dim == 2
                str_pdx = [repelem('x',tags(end-1)),repelem('t',tags(end))];
            elseif dim == 3
                str_pdx = [repelem('x',tags(end-2)),repelem('y',tags(end-1)),repelem('t',tags(end))];
            elseif dim == 4
                str_pdx = [repelem('x',tags(end-3)),repelem('y',tags(end-2)),repelem('z',tags(end-1)),repelem('t',tags(end))];
            end

            if isreal(tags)
                tags_pde{ind} = ['u^{', strrep(num2str(tags(1:n)),'  ',','),'}_{',str_pdx,'}'];
                ind = ind+1;
            else
                ind_symb = imag(tags(1:n));
                symb = sum(ind_symb);
                if symb>0
                    tags_pde{ind} = ['cos(',num2str(abs(symb)),'u^{',strrep(num2str(ind_symb),'  ',','),'})_{',str_pdx,'}'];
                    ind = ind+1;
                else
                    tags_pde{ind} = ['sin(',num2str(abs(symb)),'u^{',strrep(num2str(abs(ind_symb)),'  ',','),'})_{',str_pdx,'}'];
                    ind = ind+1;
                end
            end
%        end
    end
    lib_list = lib_list(~remove_inds,:);
end    
