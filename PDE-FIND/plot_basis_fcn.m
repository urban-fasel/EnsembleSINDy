%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: plot test functions used to identify PDE
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function plot_basis_fcn(Cfs_x,Cfs_t,supp_phi_x,dx,supp_phi_t,dt,max_dx,max_dt,pdx_list,pos,scales)
    if isempty(pos)
        pos = [0.4 0.5 0.6 0.45];
    end
    if isempty(scales)
        scales = [1 1];
    end
    clf;set(gcf,'units','normalized','outerposition',pos)
    nn= max_dx+max_dt+1;
    
    pdx_list = [pdx_list(2,:);pdx_list(1,:);pdx_list(3:end,:)];
    
    for j=1:nn
        subplot(2,ceil(nn/2),j)
        phi_x = Cfs_x(pdx_list(j,1)+1,:)' * (supp_phi_x*scales(1)*dx)^(-pdx_list(j,1));
        phi_t = Cfs_t(pdx_list(j,end)+1,:) * (supp_phi_t*scales(2)*dt)^(-pdx_list(j,end));
        surf(-scales(2)*supp_phi_t*dt:scales(2)*dt:scales(2)*supp_phi_t*dt,-scales(1)*supp_phi_x*dx:scales(1)*dx:scales(1)*supp_phi_x*dx,phi_x*phi_t, 'EdgeColor','none') 
        axis tight 
        title(['$\alpha^{',num2str(j-1),'} = (', num2str(pdx_list(j,1)),',',num2str(pdx_list(j,end)),')$'],'interpreter','latex','fontsize',14)
        view([55 15])
        xlabel('$t$','interpreter','latex','fontsize',14)
        ylabel('$x$','interpreter','latex','fontsize',14)
        set(gca, 'TickLabelInterpreter','latex','fontsize',14)
    end
end
