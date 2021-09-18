%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% SINDy central difference
%%%%%%%%%%%% 
%%%%%%%%%%%% 

function sindy = sindy_cd(xobs,Theta_0,n,lambda,gamma,dtL)

dxobs_0(1,:)=(-11/6*xobs(1,:) + 3*xobs(2,:) -3/2*xobs(3,:) + xobs(4,:)/3)/dtL;
dxobs_0(2:size(xobs,1)-1,:) = (xobs(3:end,:)-xobs(1:end-2,:))/(2*dtL);
dxobs_0(size(xobs,1),:) = (11/6*xobs(end,:) - 3*xobs(end-1,:) + 3/2*xobs(end-2,:) - xobs(end-3,:)/3)/dtL;

sindy = sparsifyDynamics(Theta_0,dxobs_0,lambda,n,gamma);

end
