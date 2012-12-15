function peakVals = springDampSolver23(m,k,b,totalTime)

%% Parameter Sweep of ODEs
% This is a parameter sweep study of a 2nd order ODE system.
%
%   m*X'' + b*X' + k*X = 0
%
% We solve the ODE for a time span of 0 to 25 seconds, with initial
% conditions X(0) = 0 and X'(0) = 1. We sweep the parameters "b" and "k"
% and record the peak values of X for each condition. At the end, we plot a
% surface of the results.

% Copyright 2009-2011 The MathWorks, Inc.

%% Solve ODE

[kGrid, bGrid] = meshgrid(b, k);
peakVals = nan(size(kGrid));

for idx = 1:numel(peakVals)
   [T, Y] = ode23(@(t,y) odesystem(t, y, m, bGrid(idx), kGrid(idx)), ...
      [0, totalTime], ...  % simulate for totalTime seconds
      [0, 1]) ;            % initial conditions
   
   peakVals(idx) = max(Y(:, 1));
end


%#ok<*ASGLU>