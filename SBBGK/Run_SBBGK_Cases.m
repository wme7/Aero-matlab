%% Run Cases for ploting Data
clc; clear all; close all; 

if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 4
end

%% Number of cases
cases = 1:4;                    % Number of Cases to evaluate

%% Simulation Parameters
name        ='SBBGK1d'; % Simulation Name
CFL         = 0.30;     % CFL condition
%r_time      = 0.0001;   % Relaxation time  <-- Parameter To be studied
tEnd        = 0.2;      % End time
theta       = 1;        % {-1} BE, {0} MB, {1} FD.
quad        = 1;        % for NC = 1 , GH = 2
method      = 1;        % for TVD = 1, WENO3 = 2, WENO5 = 3
IC_case     = 7;        % Reimann IC cases available: {1,2,3,4,5,6,7}
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
tau = [0.1 0.01 0.001 0.0001];  % Relaxation time per case

%% Excecute SBBGK in Parallel using 4 processors
parfor i = cases
    SBBGK_1d_func(name,CFL,tau(i),tEnd,theta,quad,method,IC_case, ...
        plot_figs,write_ans,P_deg,Pp,RK_stages,nx)
end
matlabpool close

%% if everything is ok then, tell me:
fprintf('All Results have been succesfully saved!\n')

%% Compute Exact solution for Comparison
% [xt,u,rho,p,e] = Exact_Riemann(IC_case);
%     % plot exact solution
%     if plot_fig == 1
%         % Plotting instructions
%         subplot (2,2,1)
%         plot(xt,u);
%         title('Plot of U v/s x/t');
%         xlabel ('x/t');
%         ylabel ('u');
%         axis tight;
%         subplot (2,2,2)
%         plot(xt,rho);
%         title('Plot of Density v/s x/t');
%         xlabel ('x/t');
%         ylabel ('rho');
%         axis tight;
%         subplot (2,2,3)
%         plot(xt,p);
%         title('Plot of Pressure v/s x/t');
%         xlabel ('x/t');
%         ylabel ('P');
%         axis tight;
%         subplot (2,2,4)
%         plot(xt,e);
%         title('Plot of Internal Energy v/s x/t');
%         xlabel ('x/t');
%         ylabel ('E');
%         axis tight;
%     end