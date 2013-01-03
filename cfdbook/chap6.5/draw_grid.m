%DRAW_GRID
% Screen plot of grid

tic

[X,Y] = meshgrid([0,cumsum(dx)],[0,cumsum(dy)]);
figure(1), clf, hold on
title('Grid','FontSize',16)
plot(X,Y,'k')	% Vertical grid lines
axis('equal')
plot(X',Y','k') % Horizontal grid lines

tijd = toc; disp(['draw_grid time = ',num2str(tijd)])
