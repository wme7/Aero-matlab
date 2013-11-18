%% AMG TEST II: Different boundary conditions
%
% We consider the linear finite element discretization of Poisson equation
% on the unstructured mesh with Dirichlet or Neumann boundary conditions. 

%%
load lakemesh
%% Dirichlet problem
[N,itStep,time,err] = amgtest(node,elem,1);
%% 
colHeaders = {'Unknowns','Iterations','Time(sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
r = showrate(N.*log(N),time,3);
xlabel('NlogN'); ylabel('Time');
title(['Complexity is (NlogN)^{' num2str(r) '}'] );

%% Neumann problem
[N,itStep,time,err] = amgtest(node,elem,0);
%% 
colHeaders = {'Unknowns','Iterations','Time(sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
r = showrate(N.*log(N),time,3);
xlabel('NlogN'); ylabel('Time');
title(['Complexity is (NlogN)^{' num2str(r) '}'] );