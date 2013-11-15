function [x,c] = transform_coord(xi,a,b)
%% Transform Normalized Local Coordinates to an expecific [a,b] range.
%**************************************************************************
% Evaluate:
%
%  $ x = (b+a)/2 + (b-a)/2 \xi $
% 
% Notation
%   x:  Global Coordinates system
%   xi: Local Coordiante system
%
% Where    a <=  x  <=  b  
%   and   -1 <=  xi <= +1
%
% Based on 
%   Parviz Moin; Fundamental of Engineering Numerical Analysis
%   2nd Edition, Cambridge Edit. Chapter 3.6 PP. 42
%
% Coded by Manuel Diaz 2012.12.05
%
%**************************************************************************
%% Compute global coordinates
x = (b+a)/2 + (b-a)/2*(xi);

%% Compute quadrature constant
c = (b-a)/2;
