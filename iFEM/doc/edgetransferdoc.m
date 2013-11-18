%% Check Transfer operator for edge
% The script file to check 2D transfer operator for edge.
% Face element and edge element share the same transfer in 2D.
%
% One can refer to 
% <../edgetransfer.pdf edge transfer> for document of derivations.
%
%  Created by Ming Wang at May, 2012. Revised at July, 2012. 

%% *Transfer operators for lowest order edge element*
%
%% Step 1: A simple uniform mesh and its uniform refinement
% A coarse grid
node = [0,0; 1,0; 1,1; 0,1];   
elem = [2,3,1; 4,1,3];         
% get edge
[elem2edge,edge] = dofedge(elem); %#ok<*ASGLU>
% plot coarse mesh
subplot(1,2,1); showmesh(node,elem); findedge(node,edge);
set(gcf,'Units','normal'); set(gcf,'Position',[0.25,0.25,0.5,0.5]);
% A fine grid obtained by uniform refine 
[nodef,elemf] = uniformrefine(node,elem);
% get edge
[elem2edgef,edgef] = dofedge(elemf);
% plot fine mesh
subplot(1,2,2); showmesh(nodef,elemf); findedge(nodef,edgef);

%% Step 2: Check the prolongation 
%  In one words, we intend to numerically verify that
%
% $$ ND_1(T_{c}) \subset ND_1(T_{f}) $$.
%
%  The function used is: u = a + b \cross x,  with a = [1 2] and b = 1.
%
u = inline('[1 - x(:,2), 2 + x(:,1)]','x');
% Approach 1: Direct interpolation to the fine grid
uI_f = edgeinterpolate(u,nodef,edgef);

% Approach 2: interpolation to coarse grid and prolongate to fine grid
uI_c = edgeinterpolate(u,node,edge);
pro = edgetransfer(elem,elemf);
u_c2f = pro*uI_c;

% Compare two interpolations
disp(' Edges     uI_f        u_c2f      uI_f-u_c2f');
disp(num2str([(1:size(edgef,1))' uI_f u_c2f uI_f-u_c2f]))
clear all;

%% *Transfer operators for lowest order face element*
%
%% Step 1: A simple uniform mesh and its uniform refinement
% A coarse grid
node = [0,0; 1,0; 1,1; 0,1];   
elem = [2,3,1; 4,1,3];         
% get edge
[elem2edge,edge] = dofedge(elem);
% plot coarse mesh
subplot(1,2,1); showmesh(node,elem); findedge(node,edge);
set(gcf,'Units','normal'); set(gcf,'Position',[0.25,0.25,0.5,0.5]);
% A fine grid obtained by uniform refine 
[nodef,elemf] = uniformrefine(node,elem);
% get edge
[elem2edgef,edgef] = dofedge(elemf);
% plot fine mesh
subplot(1,2,2); showmesh(nodef,elemf); findedge(nodef,edgef);


%% Step 2: Check the prolongation 
%  In one words, we intend to numerically show that:
%
% $$ RT_0(T_{c}) \subset RT_0(T_{f}) $$.
%
%  The function used is: u = a + bx,  with a = [1 2] and b = 1.
u = inline('[1 + x(:,1), 2 + x(:,2)]','x');
% Approach 1: Direct interpolation to the fine grid
uI_f = faceinterpolate(u,nodef,elemf,'RT0');

% Approach 2: interpolation to coarse grid and prolongate to fine grid
uI_c = faceinterpolate(u,node,elem,'RT0');
pro = edgetransfer(elem,elemf); % Transfer operator
u_c2f = pro*uI_c;

% Compare the solution:
disp(' Edges     uI_f        u_c2f      uI_f-u_c2f');
disp(num2str([(1:size(edgef,1))' uI_f u_c2f uI_f-u_c2f]))