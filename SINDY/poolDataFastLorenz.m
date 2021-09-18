%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy: function for building data matrix Theta_0
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy: Galerkin-based Data-Driven Model
%%%%%%%%%%%% Selection"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz
%%%%%%%%%%%%
%%%%%%%%%%%% Modified by Urban Fasel, 2021 

function Theta_0 = poolDataFastLorenz(xobs)

ind = 0;

Theta_0 = zeros(1,9);

monom_powers = [1,0,0;0,1,0;0,0,1];
num_monoms = 3; 
for j = 1:num_monoms
    Theta_0(:,ind+1) = prod(xobs.^(monom_powers(j,:)),2);
    ind = ind+1;
end
monom_powers = [0,1,1;1,0,1;1,1,0;0,0,2;0,2,0;2,0,0];
num_monoms = 6;
for j = 1:num_monoms
    Theta_0(:,ind+1) = prod(xobs.^(monom_powers(j,:)),2);
    ind = ind+1;
end

end
