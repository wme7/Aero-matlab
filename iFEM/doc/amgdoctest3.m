%% AMG TEST III: Robustness to time discretization
%
% We consider the linear finite element discretization of heat equation on
% the unstructured mesh with Neumann boundary conditions. We test the
% implicit time discretization with various time stepsizes dt = 1/h^4,
% 1/h^2, 1/h, 1.

%%
load lakemesh
%% Time step: h^4. Mass matrix will dominate.
close all;
dt = 'h^4';
[N,itStep,time,err] = amgtest(node,elem,2,[],dt);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
r = showrate(N,time,3);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r,3) '}'] ,'Fontsize', 14);

%% Time step: h^2
close all;
dt = 'h^2';
[N,itStep,time,err] = amgtest(node,elem,2,[],dt);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
r = showrate(N,time,3);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r,3) '}'] ,'Fontsize', 14);

%% Time step:  h
close all;
dt = 'h';
[N,itStep,time,err] = amgtest(node,elem,2,[],dt);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
r = showrate(N,time,3);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r,3) '}'] ,'Fontsize', 14);

%% Time step:  1. Like regularized Neumann problem.
close all;
dt = '1';
[N,itStep,time,err] = amgtest(node,elem,2,[],dt);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
r = showrate(N,time,3);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r,3) '}'] ,'Fontsize', 14);