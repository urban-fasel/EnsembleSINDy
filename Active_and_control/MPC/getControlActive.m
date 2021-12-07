function uopt = getControlActive(DD,dt,Xi,x0,uopt0,select_model,polyorder,usesine)

Nvec        = DD;%14;           % Choose prediction horizon over which the optimization is performed
% Imax        = 0.05;         % Constraint: infectious population <= zMax  
Ts          = dt;%1;            % Sampling time
n = 3;
% Ton         = 7;            % Time when control starts
% Tcontrol    = 7;            % Time when control updates: once a week
% Duration    = T_opt;        % Run control for 100 time units
Q           = 1000*[1 1 1];	%[0 0 1 0];	% State weights
R           = 0;%0.01;%0.1;          % Control variation du weights
Ru          = 0;%0.01;%0.1;          % Control weights
x0n         = x0';          % Initial condition
% uopt0       = p.beta0;      % Set initial control input to zero

% Constraints on control optimization
LB          = -15*ones(Nvec,1);%0.3*p.beta0*ones(Nvec,1);	% Lower bound of control input
UB          = 15*ones(Nvec,1);%p.beta0*ones(Nvec,1);     % Upper bound of control input

% Reference state, which shall be achieved: only deviation from xref1(3) (infectious cases) is penalized
xref1 = [0; 0; 0; 0];

% MPC parameters used in objective and constraint function
if strcmp(select_model,'DMD')
    pMPC.sys = sysmodel_DMDc;
elseif strcmp(select_model,'SINDy')
    pMPC.ahat = Xi(:,1:n,:);
    pMPC.polyorder = polyorder;
    pMPC.usesine = usesine;
end
pMPC.Ts = Ts;
% pMPC.beta0 = p.beta0;
% pMPC.Imax = Imax;   
    

%% run MPC    
% Options for optimization routine
options = optimoptions('fmincon','Display','none');%,'MaxIterations',100);

% Start simulation
x        = x0n;     % Initial state
uopt     = uopt0.*ones(Nvec,1); % Initial control input
% xHistory = x;       % Stores state history
% uHistory = uopt(1); % Stores control history
% tHistory = 0;       % Stores time history
% rHistory = xref1;   % Stores reference (could be trajectory and vary with time)


% For each iteration: take measurements & optimize control input & apply control input
% for ct = 1:(Duration/Ts)   
%     if ct*Ts>Ton            % Turn control on
%         if mod(ct*Ts,Tcontrol) == 0  % Update (Optimize) control input: once a week
            % Set references
            xref = xref1;

            % NMPC with full-state feedback
            COSTFUN = @(u) objectiveFCNactive(u,x,Nvec,xref,uopt(1),pMPC,diag(Q),R,Ru,select_model); % Cost function
%             if constrON
%                 CONSFUN = @(u) constraintFCN(u,x,Nvec,pMPC,select_model); % Constraint function
%                 uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,options); % Optimization with constraint
%             else
                uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,[],options); % Optimization without constraint
%                 uopt = fmincon(COSTFUN,uopt,[],[],[],[],[],[],[],options); % Optimization without constraint
%             end
%         end
%     else    % If control is off
%         uopt = uopt0.*ones(Nvec,1);
%         xref = [nan; nan; nan; nan];
%     end
% 
%     % Integrate system: Apply control & Step one timestep forward
%     x = rk4u(@SEIR,x,uopt(1),Ts/1,1,[],p);      % Runge-Kutta scheme of order 4 for control system
%     xHistory = [xHistory x];                    % Store current state
%     uHistory = [uHistory uopt(1)];              % Store current control
%     tHistory = [tHistory tHistory(end)+Ts/1];   % Store current time
%     rHistory = [rHistory xref];                 % Store current reference
% end