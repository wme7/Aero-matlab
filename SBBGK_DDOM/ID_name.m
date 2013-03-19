function [ID, IDn] = ID_name(name,theta,nx,P_deg,RK_stages,tau,IC_case,fmodel,f_case,method)
%% ID name generator 
% Generates ID and IDn for an specific simulation using the controling 
% parameters.

%% Read Parameters
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
switch fmodel
    case{1} %UU
        feq = '-UU';
    case{2} %ES
        feq = '-ES';
    otherwise
        error('Feq case not available')
end
% Method
switch method
    case{1} % Upwind
        advec = 'Upwind';
        P_degree = num2str(1);
    case{2} % TVD
        advec = 'TVD';
        P_degree = num2str(1);
    case{3} % WENO3
        advec = 'WENO3';
        P_degree = num2str(1);
    case{4} % WENO5
        advec = 'WENO5';
        P_degree = num2str(1);
    case{5} % DG
        advec = 'DGM';
        % Polinomial Degree
        P_degree = num2str(P_deg);
    otherwise
        error('Advection method not available')
end

% Elements used
elements  = ['X',num2str(nx)];

% RK stages
RKs = num2str(RK_stages);

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
    elements,'P',P_degree,'RK',RKs,' ',omega,'IC',ic,f];
ID = [name1,feq,name2,'-',statistic,advec,name3,'-',...
    elements,'P',P_degree,'RK',RKs,'-',omega,'IC',ic];
