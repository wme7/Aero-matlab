%% Generate Table of ENO/WENO Stencils coeficients
%**************************************************************************
% Based on:
% Chi-Wang Shu's Lectures notes on: 'ENO and WENO schemes for Hyperbolic
% Conservation Laws' 
%
% coded by Manuel Diaz, 02.09.2012, NTU Taiwan.
% Compare with Eqs. 2.20 & 2.21 for uniform grids!
%**************************************************************************
clear all; close all; clc;

%% Build Table 2.21.
c = zeros(8,7,7);
for k = 1:7
    for i = 1:k+1 % dummy index 
        r = i-2; % for range -1 to k-1
        for l = 1:7 % dummy index 
            j = l-1; % for range 0 to 6
            c(i,l,k) = c_rj(r,j,k);
        end
    end
end
