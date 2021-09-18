%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: get test function values on scaled reference grid,
%%%%%%%%%%%% computed from max derivative (d) and half-support size m, such
%%%%%%%%%%%% that 1D test function phi is supported on 2m+1 points
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function [Cfs,p] = phi_int_weights(m,d,tol,phi_class)

    if phi_class == 1
        if tol<0
            p = -tol;
        else
            p = ceil(max(log(tol)/log((2*m-1)/m^2),d+1));     % choose p so that the penultimate grid-point has value eps (for decay)
        end
        t = (0:m)/m;
        t_L = zeros(d+1,m+1);                             % store (1+t)^q, (1-t)^q
        t_R = zeros(d+1,m+1);                                  
        for j=1:m
            t_L(:,j)  = (1+t(j)).^(fliplr(p-d:p))';         
            t_R(:,j)  = (1-t(j)).^(fliplr(p-d:p))';
        end

        ps = ones(d+1,1);                                  % derivative coefficients
        for q=1:d
            ps(q+1) = (p-q+1)*ps(q);
        end
        t_L = ps.*t_L;
        t_R = ((-1).^(0:d)'.*ps).*t_R;

        Cfs = zeros(d+1,2*m+1);                            % Values of derivatives at grid points
        Cfs(1,:) = [fliplr(t_L(1,:).*t_R(1,:)) t_L(1,2:end).*t_R(1,2:end)];
        P = fliplr(pascal(d+1));    
        for k=1:d
            binoms = diag(P,d-k);
            Cfs_temp = zeros(1,m+1);
            for j=1:k+1
                Cfs_temp = Cfs_temp + binoms(j)*t_L(k+2-j,:).*t_R(j,:);
            end
            Cfs(k+1,:) = [(-1)^k*fliplr(Cfs_temp) Cfs_temp(2:end)];
        end
    elseif phi_class == 2
        if and(tol < 0,m>0)
            p = -tol;
            a = m*p;
        elseif and(tol > 0, m>0)
            a = sqrt(-2*log(tol));
            p = (1-1/m)/a;
        elseif and(tol > 0, m <= 0)
            p = -m;
            a = sqrt(-2*log(tol));
            m = ceil(1+a*p);
            p = p / m;
        end
        x = linspace(-a,a,2*m-1);
        dx = x(2)-x(1);
        x=[x(1)-dx x x(end)+dx];

        Cfs = ones(d+1,2*m+1);
        Cfs(2,:) = x;

        for k=3:d+1
            Hnp = (k-2)*Cfs(k-2,:);
            Cfs(k,:) = x.*Cfs(k-1,:)-Hnp;
        end

        e = exp(-x.^2/2);
        s = (-p*m).^((0:d)');

        Cfs = Cfs.*(s*e);
    end
    
end    
    
    
    
    