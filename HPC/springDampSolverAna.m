function peakVals = springDampSolverAna(m,k,b,totalTime)

% Copyright 2011 The MathWorks, Inc.

t = (0:0.1:totalTime)';  % Define a vector to evaulate on

%% Solving Analytically
% The following is the analytical form of the solution

[kGrid, bGrid] = meshgrid(b, k);
peakVals = nan(size(kGrid));

for idx = 1:numel(peakVals)
   t2 = 1./m;
   t3 = bGrid(idx).^2;
   t6 = 4.*kGrid(idx).*m;
   t4 = t3 - t6;
   t9 = (bGrid(idx).*t.*t2)./2;
   t5 = 1./exp(t9);
   t7 = t4.^(1./2);
   t8 = (t.*t2.*t7)./2;
   t10 = sinh(t8);
   t11 = 1./t4.^(1./2);
   y(:, 1) = 2.*m.*t10.*t11.*t5;
   y(:, 2) = t5.*cosh(t8) - bGrid(idx).*t10.*t11.*t5;
   
   % Get real solution
   y = real(y);
   
   peakVals(idx) = max(y(:, 1));   
end