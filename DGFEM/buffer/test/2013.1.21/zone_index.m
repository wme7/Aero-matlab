function [euler_index,sbbkg_index] = zone_index(x,xa,xb,dx)
% Compute operation boolean zones
euler_zone = x>(xa-dx);
sbbkg_zone = x<(xb+dx);
% Compute indexes
euler_index = find(euler_zone);
sbbkg_index = find(sbbkg_zone);
