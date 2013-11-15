%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Test 1-D wave equation with CPR/FR implementation
%
%               du/dt + df/dx = 0,  for x \in [a,b]
%                 where f = f(u): linear/nonlinear
%
%              coded by Manuel Diaz, NTU, 2013.11.29
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: A flux reconstruction approach to high-order schemes including
% Discontinuous Galerkin methods. H.T. Huynh, AIAA 2007.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

% Fixed Parameters
tEnd = 2*pi; % One cycle for every test
IC = 2; % sinusoidal function

% Parameters
K   = [3,4,5];
cfl = [0.01,0.001,0.0001]; % to ensure stability
nE  = [20,40,80,160,320];
%quad = {'Legendre','LGL'};

% Number of parameters
p1 = length(K);
p2 = length(nE);

% Allocate space for results
table = zeros(p2,2,p1);
Norm = zeros(size(table));
Deg = zeros(size(table));

%% Run Tests
for l = 1:p1
    for n = 1:p2
        tic;
        [Norm(n,1,l),Deg(n,1,l),Norm(n,2,l),Deg(n,2,l)] = ...
            TestCPRfun('linear',cfl(l),tEnd,IC,K(l),nE(n),'LGL');
        toc;        
    end
end

%% Display Result
for l = 1:p1
    fprintf('***************************************************************\n')
    fprintf(' Degree %d\n',K(l));
    fprintf('***************************************************************\n')
    fprintf(' nE \t L1-Norm \t Degree \t Linf-Norm \t Degree\n');
    for n = 1:p2
        fprintf('%3.0f \t %1.2e \t %2.2f \t\t %1.2e \t %2.2f \n',...
            nE(n),Norm(n,1,l),Deg(n,1,l),Norm(n,2,l),Deg(n,2,l));
    end
end
fprintf('\n');
% By observing the degree of accuaracy with respect to the CFL number, it
% is suggested that we should tune the CFL condition for each case in order
% to get the fastest, stable and most accuarate solution.

% Manuel Diaz, NTU, 2013
% End of Test