% My first program
% other info

clear;

A = magic(4)

% sum row elements
sum(A,2)

% sum columns elements
sum(A,1)

% square matrix's elements
A.^2

% compute natural logarithm of 2
% according with Wiki: ln(x) = (x-1)-(x-1)^2/2 + (x-1)^3/3 -(x-1)^4/4 ...
% so
x = 2

ln_approx = (x-1) 
ln_approx = (x-1) - (x-1)^2/2 
ln_approx = (x-1) - (x-1)^2/2 + (x-1)^3/3 
ln_approx = (x-1) - (x-1)^2/2 + (x-1)^3/3 -(x-1)^4/4
ln_approx = (x-1) - (x-1)^2/2 + (x-1)^3/3 -(x-1)^4/4 + (x-1)^5/5 
ln_approx = (x-1) - (x-1)^2/2 + (x-1)^3/3 -(x-1)^4/4 + (x-1)^5/5 -(x-1)^6/6

ln_exact = log(2)

% hey not bad! ; )