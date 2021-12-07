function J = objectiveFCNlorenz(u,x,N,xref,u0,p,Q,R,Ru)
%% Cost function of nonlinear MPC
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1
%   x:      current state at time k
%   N:      prediction horizon
%   xref:   state references, constant from time k+1 to k+N
%   u0:     previous controller output at time k-1
%   p:      model parameters
%   Q,R,Ru: MPC cost function weights      
%
% Output:
%   J:      objective function cost
%

%% Integrate system
Ns = size(x,1);
xk = zeros(Ns,N+1); xk(:,1) = x;
for ct=1:N
    % Obtain plant state at next prediction step.
    xk(:,ct+1) = rk4u(@sparseGalerkinControl,xk(:,ct),u(ct),p.Ts,1,[],p);
end
xk = xk(:,2:N+1);


%% Cost Calculation
% Set initial plant states, controller output and cost
uk = u(1);
J = 0;
% Loop through each prediction step
for ct=1:N
    % Obtain plant state at next prediction step
    xk1 = xk(:,ct);
    % Accumulate state tracking cost from x(k+1) to x(k+N)
    J = J + (xk1-xref)'*Q*(xk1-xref);
    % Accumulate MV rate of change cost from u(k) to u(k+N-1)
    if ct==1
        J = J + (uk-u0)'*R*(uk-u0) + (p.beta0-uk)'*Ru*(p.beta0-uk);
    else
        J = J + (uk-u(ct-1))'*R*(uk-u(ct-1)) + (p.beta0-uk)'*Ru*(p.beta0-uk);
    end
    % Update uk for the next prediction step
    if ct<N
        uk = u(ct+1);
    end
end

