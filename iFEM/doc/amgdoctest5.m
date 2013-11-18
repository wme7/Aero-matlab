%% AMG TEST V: DIFFERENT INTERPOLATION OPERATORS
% 
% We consider the effect of using different interpolation operators to
% interpolate fine grid values by coarse grid values. It is tested through
% the following choices in option.interpolation
%
% * 's' standard interpolation. Use the matrix A_fc as a weighted average
% of all connected coarse nodes.
% * 't' two-points interpolation. Use at most two connected coarse
% nodes.
% * 'a' aggegration (one-point) interpolation. Use the strongest connected
% coarse node.

%%
clear all; close all;
%% Unstructured mesh in 2-D
load lakemesh
showmesh(node,elem);
%%
% Standard interpolation
option.interpolation = 's';
[N,itStep,time,err] = amgtest(node,elem,[],option);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
clf;
r = showrate(N,time,1);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r) '}'],'Fontsize', 14);
%%
% Two-points interpolation
load lakemesh
option.interpolation = 't';
[N,itStep,time,err] = amgtest(node,elem,[],option);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
clf;
r = showrate(N,time,1);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r) '}'],'Fontsize', 14);
%%
% One point interpolation
load lakemesh
option.interpolation = 'a';
[N,itStep,time,err] = amgtest(node,elem,[],option);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
clf;
r = showrate(N,time,1);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r) '}'],'Fontsize', 14);

%% Unstructured mesh in 3-D
clear all; close all
load bunny;
showboundary3(node,elem);
view([-179 74]);
%%
% Standard interpolation
option.interpolation = 's';
[N,itStep,time,err] = amgtest3(node,elem,1,option);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
clf;
r = showrate(N,time,1);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r) '}'],'Fontsize', 14);
%%
% Two-points interpolation
load bunny
option.interpolation = 't';
[N,itStep,time,err] = amgtest3(node,elem,1,option);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
clf;
r = showrate(N,time,1);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r) '}'],'Fontsize', 14);
%%
% One point interpolation
load bunny
option.interpolation = 'a';
[N,itStep,time,err] = amgtest3(node,elem,1,option);
%% 
colHeaders = {'Unknowns','Iterations','Time (sec)','Error'};
makeHtmlTable([N itStep time err],[],[],colHeaders,[],6);
%%
clf;
r = showrate(N,time,1);
xlabel('N'); ylabel('Time');
title(['Complexity is N^{' num2str(r) '}'],'Fontsize', 14);

%% Conclusion
% In general, if we use more coarse grids, we get more accurate
% interpolation and require less iteration steps. But the matrix on coarse
% level becomes denser, see nnz/Nc, which requires, indeed, more
% computational time in both smoother and coarse grid solver. In the
% extremly case, the one-point interpolation almost keeps the sparsity but
% iteration steps increase as level increases in the speed of J. The
% average using all neighboring points keeps the iteration steps but the
% sparsity increase a lot, especially in 3-D. The two-points interpolation
% seems a good balance.
