%% AMG TEST II: Different boundary conditions
%
% We consider the linear finite element discretization of Poisson equation
% on the unstructured mesh with Dirichlet or Neumann boundary conditions. 

%%
load lakemesh
%% Dirichlet problem
close all;
[N,itStep,time,err] = amgtest(node,elem,1);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
r = showrate(N,time,2);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r) '}'] ,'Fontsize', 14);

%% Neumann problem
close all;
[N,itStep,time,err] = amgtest(node,elem,0);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
r = showrate(N,time,2);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r) '}'] ,'Fontsize', 14);