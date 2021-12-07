function [xHistory, dxMPC, uHistory, J] = runMPCensemble(dt,Xi,x0,polyorder,usesine,Beta,n)

%% Initialize MPC
Nvec        = 5;        % Choose prediction horizon over which the optimization is performed
Imax        = 0.05;     % Constraint: infectious population <= zMax  
Ts          = 0.01;     % Sampling time
Duration    = 4;        % Run control for 100 time units
Q           = [1 1 1]; 	% State weights
R           = 0.001;    % Control variation du weights
Ru          = 0.001;    % Control weights
x0n         = x0;       % Initial condition
uopt0       = 0;        % Set initial control input to zero

% Constraints on control optimization
LB          = -50*ones(Nvec,1);%0.3*p.beta0*ones(Nvec,1);	% Lower bound of control input
UB          = 50*ones(Nvec,1);%p.beta0*ones(Nvec,1);     % Upper bound of control input

% Reference state, which shall be achieved: only deviation from xref1(3) (infectious cases) is penalized
% xref1 = [-sqrt(72);-sqrt(72);27];
xref1 = [sqrt(72);sqrt(72);27];

% MPC parameters used in objective and constraint function
pMPC.ahat = Xi(:,1:n);
pMPC.polyorder = polyorder;
pMPC.usesine = usesine;
pMPC.Ts = Ts;
pMPC.beta0 = 0;
pMPC.Imax = Imax;   
    

%% run MPC    
% Options for optimization routine
options = optimoptions('fmincon','Display','none');%,'MaxIterations',100);

% Start simulation
x        = x0n;     % Initial state
uopt     = uopt0.*ones(Nvec,1); % Initial control input
xHistory = x;       % Stores state history
uHistory = uopt(1); % Stores control history
tHistory = 0;       % Stores time history
rHistory = xref1;   % Stores reference (could be trajectory and vary with time)


% For each iteration: take measurements & optimize control input & apply control input
for ct = 1:(Duration/Ts)   

    xref = xref1;

    % NMPC with full-state feedback
    COSTFUN = @(u) objectiveFCNlorenz(u,x,Nvec,xref,uopt(1),pMPC,diag(Q),R,Ru); % Cost function
    try
        uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,[],options); % Optimization without constraint
    catch
        uopt = 0;
    end

    % Integrate system: Apply control & Step one timestep forward
    x = rk4u(@lorenzForcing,x,uopt(1),Ts/1,1,[],Beta);      % Runge-Kutta scheme of order 4 for control system
    
    xHistory = [xHistory x];                    % Store current state
    uHistory = [uHistory uopt(1)];              % Store current control
    tHistory = [tHistory tHistory(end)+Ts/1];   % Store current time
    rHistory = [rHistory xref];                 % Store current reference
end

dxMPC = xHistory - xref;

% Evaluate objective function
J = evaluateObjectiveFCN(uHistory,xHistory,rHistory,diag(Q),R,Ru,pMPC);
