function [ID, IDn] = ID_name(input)
%% ID name generator 
% Generates ID and IDn for an specific simulation using the controling 
% parameters.

%% Read Parameters
    name    = input.name;           % Simulation Name
    f_case  = input.f_case;         % {1}Relaxation model, {2}Euler Limit
    theta   = input.theta;          % {-1}BE, {0}MB, {1}FD.
    feq_model = input.feq_model;    % {1}UU.  {2}ES.
    tau     = input.r_time;         % Relaxation time (if f_case = 1)
    method  = input.method;         % {1}CPR, {2}CPR, {3}CPR,
    IC_case = input.IC;             % IC:{1}~{14}. See SBBGK_IC1d.m
    P_deg   = input.P;              % Polinomial Degree
    nx      = input.K;              % Number of elements
    RK_degree = input.RK_degree;    % Number of RK stages
    RK_stages = input.RK_stages;    % Number of RK stages

% Name
name1 = name(1:2);
name2 = name(3:end-1);
name3 = name(end-1:end);

% Statistic Used
switch theta
    case {-1}
        statistic = 'BE';
    case {0}
        statistic = 'MB';
    case {1}
        statistic = 'FD';
    otherwise
        error('not a valid theta value')
end
% Model of the Equilibrium Distribution
switch feq_model
    case {1} %UU
        feq = '-UU';
    case {2} %ES
        feq = '-ES';
    otherwise
        error('Feq case not available')
end
% Method
switch method
    case {1} % NDG
        advec = 'NDG';
        % Polinomial Degree
        P_degree = num2str(P_deg);
    case {2} % MDG
        advec = 'MDG';
        % Polinomial Degree
        P_degree = num2str(P_deg);
    case {3} % CPR/FR
        advec = 'CPR';
        % Polinomial Degree
        P_degree = num2str(P_deg);
    otherwise
        error('Advection method not available')
end

% Elements used
elements  = ['X',num2str(nx)];

% RK stages
RKs = num2str(RK_stages);
RKd = num2str(RK_degree);

% Initial Condition Number
ic = num2str(IC_case);

% Relaxation model or Euler limit
switch f_case
    case{1} % Relaxation time
        % Compute Relaxation frequency: 'omega',
        omega = 1/tau;
        omega = ['w',num2str(omega)];
    case{2} % Euler Limit
        omega = 'EL';
    otherwise
        error('Feq case not available')
end

% Tecplot format
f = '.plt';

%% Generate ID
IDn = [name1,feq,name2,' ',statistic,advec,name3,' ',...
    elements,'P',P_degree,'RK',RKd,RKs,' ',omega,'IC',ic,f];
ID = [name1,feq,name2,'-',statistic,advec,name3,'-',...
    elements,'P',P_degree,'RK',RKd,RKs,'-',omega,'IC',ic];
