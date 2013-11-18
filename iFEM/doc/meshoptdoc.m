%% Mesh Smoothing and Optimization
% 
%% Improve geometric mesh quality
% 
% The function [node,elem] = optmesh(node,elem) will optimize the shape
% regularity of triangles in the input mesh (node,elem) and outputs a
% better mesh (node,elem). 

clear all; close all;
load airfoilperturb
figure(1); subplot(1,2,1); 
showmesh(node,elem); title('original mesh');
figure(2); subplot(1,2,1); 
showmeshquality(node,elem); axis([0 1 0 2700]);
[node,elem] = optmesh(node,elem);
figure(1); subplot(1,2,2); 
showmesh(node,elem); title('smoothed mesh');
figure(2); subplot(1,2,2); 
showmeshquality(node,elem); axis([0 1 0 2700]);
%%
% We explain algorithms implemented in optimesh.m in the following.

%% ODT-based mesh smoothing
%
% In the function <../../../mesh/html/meshsmoothing.html meshsmoothing>, we move one node at a time inside its
% patch, which consists of all triangles surrounding this node, such that
% the interpolation error to a quadratic function is minimized. The
% function |meshsmoothing| will keep the topology of the input mesh, i.e.,
% the node index and connectivity of nodes are unchanged.
%
% In the simplest case, the scheme is to move the node to the average of
% circumenters of triangles in the local patch. Details can be found in the
% paper <http://math.uci.edu/~chenlong/CH2008.html ODTmesh>. 
theta = [-2*pi/3 -pi/3 0 pi/3 2*pi/3 pi]';
node = [cos(theta), sin(theta)];
node(end+1,:) = 0;
elem = delaunayn(node);
node(end,:) = rand(1,2)*0.4;
figure(1); subplot(1,2,1);
showmesh(node,elem); findnode(node,'all','noindex');
c = circumcenter(node,elem);
hold on; plot(c(:,1),c(:,2),'r.','MarkerSize',16)
node(end,:) = mean(c);
plot(node(end,1),node(end,2),'b.','MarkerSize',16)
figure(1); subplot(1,2,2);
showmesh(node,elem); findnode(node,'all','noindex');

