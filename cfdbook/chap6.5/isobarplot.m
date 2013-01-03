%ISOBARPLOT
% Screen plot of isobars

figure(2), clf, hold on
title('Isobars','FontSize',16)
pp = reshape(p0,size(XP')); pp = pp';
cvals = linspace(min(min(pp)),max(max(pp)),10);
contour(xv,yu,pp,cvals,'k')
hold off
