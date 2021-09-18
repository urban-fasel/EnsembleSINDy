function [Xi,its,thrs_EL] = sparsifyDynamics(Theta,dXdt,lambda,gamma,M)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz
%
% modified by Daniel A. Messenger, 2020 to prevent return of zero vector
% and include regularization
%
% compute Sparse regression: sequential least squares

n = min(size(dXdt));
nn = size(Theta,2);

if  gamma ~= 0
    Theta = [Theta;gamma*eye(nn)];
    dXdt = [dXdt;zeros(nn,n)];
end

Xi = Theta \ dXdt;  % initial guess: Least-squares
if ~isempty(M)
    Xi = M.*Xi;
end

if isempty(M)
    thrs_EL = lambda;
else
    bnds = norm(dXdt)./vecnorm(Theta)'.*M; 
    LBs = lambda*max(1,bnds);
    UBs = 1/lambda*min(1,bnds);
    thrs_EL = [LBs bnds UBs];
end

smallinds = 0*Xi;
for j=1:nn
    if ~isempty(M)
        smallinds_new = or(abs(Xi)<LBs,abs(Xi)>UBs);
        if all(smallinds_new(:)==smallinds(:))
            its = j;
            return
        else
            smallinds = smallinds_new;
            Xi(smallinds)=0;    
            for ind=1:n
                Xi(~smallinds,ind) = M(~smallinds).*(Theta(:,~smallinds)\dXdt(:,ind));
            end
        end
    else
        smallinds_new = (abs(Xi)<lambda);
        if all(smallinds_new(:)==smallinds(:))
            its = j;
            return
        else
            smallinds = smallinds_new;
            Xi(smallinds)=0;
            for ind = 1:n        
                biginds = ~smallinds(:,ind);
                Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind);
            end
        end
    end
end
its = nn;
end

%    bnds = norm(dXdt_reg).^2./abs(dXdt_reg'*Theta_reg)'.*M; 


%             while size(find(smallinds)) == size(Xi(:))
%                 lambda = lambda/2;
%                 smallinds = (abs(Xi)<lambda);
%             end
