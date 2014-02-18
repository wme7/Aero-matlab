% This is something I always wondered...
% Manuel Diaz, NTU, 2014.01.21

%
A = magic(6);

% 
idx = 1:6; 
idx_new = repmat(idx,[6,1]);

% 
x = linspace(0,1,6);
v = sin(2*pi*x);

% operation
b = A.*v(idx_new); 
disp(b);

% conclusion:
% matlab indexes = c++ pointers ;)
