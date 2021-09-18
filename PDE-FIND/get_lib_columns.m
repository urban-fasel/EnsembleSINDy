%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: compute G and b, listed together in Theta_pdx, 
%%%%%%%%%%%% using convolutions over separable test functions via convNDfft
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function Theta_pdx = get_lib_columns(n,lib_list,U_obs,Cfs_x,Cfs_t,supp_phi_x,supp_phi_t,dx,dt,sub_inds,dim,scales)

Ns = size(U_obs{1});

if isempty(scales)
    scales = ones(n+dim,1);
end

Cfs_ffts = cell(dim,1);
[mm,nn] = size(Cfs_x);
for k=1:dim-1
    Cfs_ffts{k} = [zeros(mm,Ns(k)-nn) (supp_phi_x*dx*scales(n+k)).^(-(0:mm-1)').*Cfs_x/nn];
    Cfs_ffts{k} = fft(Cfs_ffts{k},[],2);
end
[mm,nn] = size(Cfs_t);
Cfs_ffts{dim} = [zeros(mm,Ns(dim)-nn) (supp_phi_t*dt*scales(n+dim)).^(-(0:mm-1)').*Cfs_t/nn];
Cfs_ffts{dim} = fft(Cfs_ffts{dim},[],2);
    
Theta_pdx = [];
ind = 1;
    
while ind<size(lib_list,1)+1
    tags = lib_list(ind,1:n);
    if isreal(tags)
        fcn = (U_obs{1} / scales(1)).^tags(1);
        for k=2:n
            fcn = fcn.*((U_obs{k} / scales(k)).^tags(k));
        end
    else
        ind_freq = find(imag(tags(1:n)));
        freq = sum(imag(tags(1:n)));        
        if freq<0
            fcn = sin(abs(freq)*U_obs{ind_freq});
        else
            fcn = cos(abs(freq)*U_obs{ind_freq});
        end
    end

    while all(lib_list(ind,1:n) == tags)
        test_conv_cell = {};
        for k=1:dim-1
            test_conv_cell{k} = Cfs_ffts{k}(lib_list(ind,n+k)+1,:);
        end
        test_conv_cell{dim} = Cfs_ffts{dim}(lib_list(ind,n+dim)+1,:);
        fcn_conv = convNDfft(fcn,test_conv_cell,sub_inds,2);
        Theta_pdx(:,ind) = fcn_conv(:);
        ind = ind+1;
        if ind > size(lib_list,1)
            break
        end
    end
end
end



% function Theta_pdx = get_lib_columns(n,lib_list,U_obs,Cfs_x,Cfs_t,supp_phi_x,supp_phi_t,dx,dt,sub_inds,dim,scales)
% 
% Theta_pdx = [];
% ind = 1;
% 
% if isempty(scales)
%     while ind<size(lib_list,1)+1
%         tags = lib_list(ind,1:n);
%         if isreal(tags)
%             fcn = U_obs{1}.^tags(1);
%             for k=2:n
%                 fcn = fcn.*(U_obs{k}.^tags(k));
%             end
%         else
%             ind_freq = find(imag(tags(1:n)));
%             freq = sum(imag(tags(1:n)));        
%             if freq<0
%                 fcn = sin(abs(freq)*U_obs{ind_freq});
%             else
%                 fcn = cos(abs(freq)*U_obs{ind_freq});
%             end
%         end
% 
%         while all(lib_list(ind,1:n) == tags)
%             test_conv_cell = {};
%             for k=1:dim-1
%                 test_conv_cell{k} = Cfs_x(lib_list(ind,n+k)+1,:)' * ((supp_phi_x*dx)^(-lib_list(ind,n+k))) * 1/(2*supp_phi_x+1);
%             end
%             test_conv_cell{dim} = Cfs_t(lib_list(ind,end)+1,:)' * ((supp_phi_t*dt)^(-lib_list(ind,end))*dt) * 1/(2*supp_phi_x+1);
%             fcn_conv = convNDfft(fcn,test_conv_cell,sub_inds,1);
%             Theta_pdx(:,ind) = fcn_conv(:);
%             ind = ind+1;
%             if ind > size(lib_list,1)
%                 break
%             end
%         end
%     end
% else
%     while ind<size(lib_list,1)+1
%         tags = lib_list(ind,1:n);
%         if isreal(tags)
%             fcn = (U_obs{1} / scales(1)).^tags(1);
%             for k=2:n
%                 fcn = fcn.*((U_obs{k} / scales(k)).^tags(k));
%             end
%         else
%             ind_freq = find(imag(tags(1:n)));
%             freq = sum(imag(tags(1:n)));        
%             if freq<0
%                 fcn = sin(abs(freq)*U_obs{ind_freq});
%             else
%                 fcn = cos(abs(freq)*U_obs{ind_freq});
%             end
%         end
% 
%         while all(lib_list(ind,1:n) == tags)
%             test_conv_cell = {};
%             for k=1:dim-1
%                 test_conv_cell{k} = Cfs_x(lib_list(ind,n+k)+1,:)' * ((supp_phi_x*scales(n+k)*dx)^(-lib_list(ind,n+k)) * 1/(2*supp_phi_x+1));
%             end
%             test_conv_cell{dim} = Cfs_t(lib_list(ind,end)+1,:)' * ((supp_phi_t*scales(end)*dt)^(-lib_list(ind,end)) * 1/(2*supp_phi_t+1));
%             fcn_conv = convNDfft(fcn,test_conv_cell,sub_inds,1);
%             Theta_pdx(:,ind) = fcn_conv(:);
%             ind = ind+1;
%             if ind > size(lib_list,1)
%                 break
%             end
%         end
%     end
% end