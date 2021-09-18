
function [scales,M_scale] = get_scales(scale_u,p_x,supp_phi_x,dx,max_dx,p_t,supp_phi_t,dt,max_dt,lib_list,dim,phi_class,Cfs_x,Cfs_t)

    if phi_class==1
        scale_x = (prod(p_x-(0:floor(max_dx/2)-1))/prod(1:ceil(max_dx/2))*prod(1:max_dx))^(1/max_dx) / (supp_phi_x*dx);
        scale_t = (prod(p_t-(0:floor(max_dt/2)-1))/prod(1:ceil(max_dt/2))*prod(1:max_dt))^(1/max_dt) / (supp_phi_t*dt);   % enforce unit inf norm                
    elseif phi_class==2
        pow_x = ceil(max_dx); pow_t = ceil(max_dt);
        scale_x = norm(Cfs_x(end,:),1)^(1/pow_x)/supp_phi_x/dx;
        scale_t = norm(Cfs_t(end,:),1)^(1/pow_t)/supp_phi_t/dt;
    end
    
    scales = [scale_u repmat(scale_x,1,dim-1) scale_t];
    M_scale = scales.^(-lib_list);
    M_scale(imag(M_scale)~=0)=1;
    M_scale = prod(M_scale,2);
    
end
