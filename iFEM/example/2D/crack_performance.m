function [N,errH1,cost] = crack_performance(maxN,theta,estimateoption)
%% CRACK Problem
%
% crack_performance solves Poisson equation in a crack domain with AFEM.
% It test the performance of iFEM for large number of unknowns. The default
% setting is: 
%     maxN = 5e5;     theta = 0.5;    estimateoption = 'recovery';
%
% [N,errH1,cost] = crack_performance(1e5,0.75,'residual') use the setting
%     maxN = 1e5;     theta = 0.75;    estimateoption = 'residual';
% and returns the number of elements, H1 error and computational cost.
%
% One can vary the input arguments to test different parameters and
% estimators.
%
% Example
%
%     [N1,err1H1] = crack_performance(3e4,0.5,'recovery');
%     [N2,err2H1] = crack_performance(3e4,0.5,'residual');
%     figure(1); clf;
%     r1 = showrate(N1,err1H1,10,'-*');
%     figure(1); hold on
%     r2 = showrate(N2,err2H1,10,'-k*');
%     legend('||Du-Du_h|| recovery type error estimator',['cN^{' num2str(r1) '}'],...
%     '||Du-Du_h|| residual type error estimator',['cN^{' num2str(r2) '}'], ...
%     'LOCATION','Best');
% 
% Example    
%     [N1,err1H1] = crack_performance(3e4,0.25);
%     [N2,err2H1] = crack_performance(3e4,0.5);
%     [N3,err3H1] = crack_performance(3e4,0.75);
%     figure(2); clf;
%     r1 = showrate(N1,err1H1,10,'-*');
%     figure(2); hold on
%     r2 = showrate(N2,err2H1,10,'-k*');
%     r3 = showrate(N3,err3H1,10,'-r*');
%     legend('||Du-Du_h|| \theta = 0.25',['cN^{' num2str(r1) '}'],...
%     '||Du-Du_h|| \theta = 0.50',['cN^{' num2str(r2) '}'], ...
%     '||Du-Du_h|| \theta = 0.75',['cN^{' num2str(r3) '}'], ...
%        'LOCATION','Best');     
%
% See also  crack, Lshape
%
% Copyright (C) Long Chen.

profile on
close all
%% Parameters
switch nargin
    case 0
        maxN = 5e5;     theta = 0.5;    estimateoption = 'recovery';
    case 1
        theta = 0.5;    estimateoption = 'recovery';
    case 2
        estimateoption = 'recovery';
end
maxIt = 50; 
N = zeros(maxIt,1);     cost = zeros(maxIt,1);  errH1 = zeros(maxIt,1);

%%  Generate an initial mesh
node = [1,0; 0,1; -1,0; 0,-1; 0,0; 1,0];    % nodes
elem = [5,1,2; 5,2,3; 5,3,4; 5,4,6];        % elements
elem = label(node,elem);                    % label the mesh
bdEdge = setboundary(node,elem,'Dirichlet');    % Dirichlet boundary condition

%%  Get a fine mesh by uniform bisection
for k = 1:2
    [node,elem,bdEdge] = uniformbisect(node,elem,bdEdge);
end

%% Set up PDE data
pde.f = @f;
pde.g_D = @exactu;

%%  Adaptive Finite Element Method
for k = 1:maxIt
    t = cputime;
    % Step 1: SOLVE
    u = Poisson(node,elem,pde,bdEdge);
    % Step 2: ESTIMATE
    if strcmp(estimateoption,'recovery')
        eta = estimaterecovery(node,elem,u);       % recovery type
    else
        eta = estimateresidual(node,elem,u,pde);    % residual type
    end
    cost(k) = cputime-t; 
    errH1(k) = getH1error(node,elem,@Du,u);
    N(k) = size(elem,1);
    if (N(k)>maxN), break; end        
    % Step 3: MARK
    markedElem = mark(elem,eta,theta);
    % Step 4: REFINE
   [node,elem,bdEdge] = bisect(node,elem,markedElem,bdEdge);
end

%% Plot computational cost
N = N(1:k);	cost = cost(1:k);   errH1 = errH1(1:k);
figure(1); clf; 
r1 = showrate(N,errH1,5,'-*');
legend('||Du-Du_h||',['cN^{' num2str(r1) '}'],'LOCATION','Best');
figure(2); clf;
r2 = showrate(N(10:end),cost(10:end),5,'k-+');
title('Computational cost vs Number of elements', 'FontSize', 14);
legend('Cost for AFEM loop',['cN^{' num2str(r2) '}'],'LOCATION','Best')
profile viewer
%%
% The step SOLVE dominates the whole simulation. 
end  %End of function CRACK_PERFORMANCE

%%  Data of CRACK
function z = f(p) % load data (right hand side function)
z = ones(size(p,1),1);
end

function z = exactu(p) % exact solution
r = sqrt(sum(p.^2,2));
z = sqrt(0.5*(r-p(:,1)))-0.25*r.^2;
end

function z = Du(p) % derivative of the exact solution
r = sqrt(sum(p.^2,2));
z(:,1) = (p(:,1)./r-1)./sqrt(8*(r-p(:,1)))-0.5*p(:,1); % u_x
z(:,2) = p(:,2)./r./sqrt(8*(r-p(:,1)))-0.5*p(:,2);     % u_y
end