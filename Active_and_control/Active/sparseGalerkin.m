function dy = sparseGalerkin(t,y,ahat,polyorder)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz
%
% Modified Urban Fasel, 2021

if max(polyorder) == 2
    yPool = poolDataFastLorenz(y');
elseif max(polyorder) == 3
    yPool = poolDataFastLV(y');
end

dy = (yPool*ahat)';