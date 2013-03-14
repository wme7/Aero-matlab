function [ID, IDn] = ID_name(name,theta,nx,P_deg,RK_stages,tau,IC_case,fmodel,f_case)
%% ID name generator 
% Generates ID and IDn for an specific simulation using the controling 
% parameters.

%% Read Parameters
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
        feq = 'UU';
    case{2} %ES
        feq = 'ES';
    otherwise
        error('Feq case not available')
end
% Elements used
elements  = num2str(nx);
% Polinomial Degree
P_degree = num2str(P_deg);
% RK stages
RKs = num2str(RK_stages);
% Relaxation time
omega = 1/tau;
omega = num2str(omega);
% Initial Condition Number
ic = num2str(IC_case);
% Relaxation model or Euler limit
switch f_case
    case{1} % Relaxation Model
        % Do nothing
    case{2} % Euler Limit
        omega = 'EL';
    otherwise
        error('Feq case not available')
end
% Tecplot format
f = '.plt';

%% Generate ID
IDn = [name,statistic,feq,'X',elements,' ','P',P_degree,'RK',RKs,'w',omega,'IC',ic,f];
ID = [name,statistic,feq,'X',elements,'-','P',P_degree,'RK',RKs,'w',omega,'IC',ic];
