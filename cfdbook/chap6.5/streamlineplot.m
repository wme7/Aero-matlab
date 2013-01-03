%STREAMLINEPLOT
% Screen plot of streamlines

figure(3), clf
title('Streamlines','FontSize',16)
hold on
uq = reshape(u1,size(XU')); uq = uq';
vq = reshape(v1,size(XV')); vq = vq';
sf = zeros(size(X));		% Streamfunction
sf(1,:) = [0, - cumsum(DXV(1,:).*vq(1,:))];
for k = 2:K+1
  sf(k,:) = sf(k-1,:) + uq(k-1,:).*DYU(k-1,:);
end
cvals = [linspace(min(min(sf)),max(max(sf)),30),0,0.995*max(max(sf))];
contour(xu,yv,sf,cvals,'k')
uq = (uq(:,1:J) + uq(:,2:J+1))/2;
vq = (vq(1:K,:) + vq(2:K+1,:))/2;
%quiver(xv,yu,uq,vq,0.9,'k')
hold off
