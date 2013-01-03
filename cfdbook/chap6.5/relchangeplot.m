%RELCHANGEPLOT
% Plot of relative change per time step

figure(5), clf
yy = linspace(yv(1),yv(end),30);
semilogy(1:length(relchange),relchange,'k-')
hold on
title('Relative change per time step','FontSize',16)
hold off
