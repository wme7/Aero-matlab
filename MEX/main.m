% Main
clc; clear all; close all;
%**************************************************************************
disp('Bisection Method starts...');
tic;
% Call bisection
root = bisection2(-10,10,1e-6);
disp(['Root = ' num2str(root)]);
a = toc;
% Time to compute
disp(['Elapsed time ' num2str(a)]);
%**************************************************************************
disp('Bisection Method starts...');
tic;
% Call bisection
root = bisection2_mex(-10,10,1e-6);
disp(['Root = ' num2str(root)]);
b = toc;
% Time to compute
disp(['Elapsed time ' num2str(b)]);
%**************************************************************************
% Compare
disp(['Mex is  ' num2str(a/b) 'times faster' ]);