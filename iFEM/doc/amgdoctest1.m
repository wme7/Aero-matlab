%% AMG TEST I: DIFFERENT MESHES
% 
% We consider linear finite element discretization of the Poisson equation
% with homongenous Dirichlet boundary condition on different meshes. 

%%
clear all; close all;
%% Uniform mesh
[node,elem] = squaremesh([0,1,0,1],0.1);
[node,elem] = uniformrefine(node,elem);
[node,elem] = uniformrefine(node,elem);
showmesh(node,elem);
snapnow;
[N,itStep,time,err,errHist] = amgtest(node,elem);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
colHeaders = {'Iterations','Approximate Error', 'Residual'};
makeHtmlTable([(0:itStep(end))' errHist],[],[],colHeaders,[],6);
%%
r = showrate(N,time,2);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r) '}'] ,'Fontsize', 14);

%% Circle mesh
close all;
[node,elem] = circlemesh(0,0,1,0.2);
[node,elem] = uniformrefine(node,elem);
[node,elem] = uniformrefine(node,elem);
showmesh(node,elem);
snapnow;
[N,itStep,time,err] = amgtest(node,elem);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
r = showrate(N,time,2);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r) '}'] ,'Fontsize', 14);

%% Unstructured mesh
close all;
load lakemesh
showmesh(node,elem);
snapnow;
[N,itStep,time,err] = amgtest(node,elem);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
r = showrate(N,time,2);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r) '}'],'Fontsize', 14);