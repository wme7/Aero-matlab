%% Run Cases for ploting Data
clc; clear all; close all; 

if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 4
end

%% Number of cases
cases = 1:4;                    % Number of Cases to evaluate

%% Simulation Parameters
name        ='SBBGK1d'; % Simulation Name
CFL         = 0.05;     % CFL condition
%r_time      = 0.0001;   % Relaxation time  <-- Parameter To be studied
tEnd        = 0.10;     % End time
theta       = 0;        % {-1} BE, {0} MB, {1} FD.
quad        = 2;        % for NC = 1 , GH = 2
method      = 1;        % for TVD = 1, WENO3 = 2, WENO5 = 3
IC_case     = 1;        % Reimann IC cases available: {1,2,3,4,5,6,7}
plot_figs   = 0;        % 0: no, 1: yes please!
write_ans   = 1;        % 0: no, 1: yes please!
% Using DG
P_deg       = 0;        % Polinomial Degree
Pp          = P_deg+1;  % Polinomials Points
% Using RK integration time step
RK_stages   = 4;        % Number of RK stages
% Grid size and Elements
nx  = 100;              % Desided number of points in our domain 

%% Particular Parameters for Runing Cases
tau = [1/10 1/100 1/1000 1/10000];  % Relaxation time per case

%% Excecute SBBGK in Parallel using 4 processors
parfor i = cases
    SBBGK_1d_func(name,CFL,tau(i),tEnd,theta,quad,method,IC_case, ...
        plot_figs,write_ans,P_deg,Pp,RK_stages,nx)
end
matlabpool close

%% if everything is ok then, tell me:
fprintf('All Results have been succesfully saved!\n')
