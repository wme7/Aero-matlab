function compareResults(bVals, kVals, P1, P2, P3, P4, createplot)

% Copyright 2011 The MathWorks, Inc.

if nargin == 6
   createplot = false;
end

% Calculate differences
d12 = P1-P2;
d13 = P1-P3;
d14 = P1-P4;
d23 = P2-P3;
d24 = P2-P4;
d34 = P3-P4;

fprintf('RMSE between ODE15s and ODE23     : %0.4f\n', sqrt(mean(d12(:).^2)));
fprintf('RMSE between ODE15s and ODE45     : %0.4f\n', sqrt(mean(d13(:).^2)));
fprintf('RMSE between ODE15s and Analytical: %0.4f\n', sqrt(mean(d14(:).^2)));
fprintf('RMSE between ODE23 and ODE45      : %0.4f\n', sqrt(mean(d23(:).^2)));
fprintf('RMSE between ODE23 and Analytical : %0.4f\n', sqrt(mean(d24(:).^2)));
fprintf('RMSE between ODE45 and Analytical : %0.4f\n', sqrt(mean(d34(:).^2)));

if createplot
   fh = figure('Units', 'normalized', 'Position', [.1 .1 .8 .8]);
   ax(1) = subplot(2,3,1);
   surf(bVals, kVals, d12);
   xlabel('b'); ylabel('k'); zlabel('P1 - P2');
   
   ax(2) = subplot(2,3,2);
   surf(bVals, kVals, d13);
   xlabel('b'); ylabel('k'); zlabel('P1 - P3');
   
   ax(3) = subplot(2,3,3);
   surf(bVals, kVals, d14);
   xlabel('b'); ylabel('k'); zlabel('P1 - P4');
   
   ax(4) = subplot(2,3,4);
   surf(bVals, kVals, d23);
   xlabel('b'); ylabel('k'); zlabel('P2 - P3');
   
   ax(5) = subplot(2,3,5);
   surf(bVals, kVals, d24);
   xlabel('b'); ylabel('k'); zlabel('P2 - P4');
   
   ax(6) = subplot(2,3,6);
   surf(bVals, kVals, d34);
   xlabel('b'); ylabel('k'); zlabel('P3 - P4');
   
   lh = linkprop(ax, 'View');
   
   setappdata(fh, 'LinkPropHandle', lh);
   
   view(ax(1), 50, 30);
end

