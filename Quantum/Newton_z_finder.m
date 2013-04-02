%% Newton Method for finding z,
% Implement root finding method.

% Fugacity's initial guess,
z0 = 0.1;

% Spatial parameters,
E = 0.5; j = 1.0; n = 1.0;

% statistic
theta = 1;

% Evaluate
X_1d(z0,n,j,E,theta,1)
dX_1d(z0,n,theta,1)