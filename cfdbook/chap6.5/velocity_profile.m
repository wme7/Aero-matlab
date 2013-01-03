%VELOCITY PROFILE

% Function called:	exact_solution

if geval == 1|geval == 7  % Horizontal Poiseuille flow; plot of outflow profile
  figure(4), clf
  yy = linspace(yv(1),yv(end),30);
  plot(exact_solution(0, 0, yy, 'right'),yy,'k-')
  hold on
  title('Velocity profiles at inflow and outflow boundaries','FontSize',16)
  uq = reshape(u1,size(XU')); uq = uq';
  plot(uq(:,J+1),yu,'o')
  plot(uq(:,1),yu,'*')
  hold off
elseif geval == 2	% Vertical Poiseuille flow;   plot of outflow profile
  figure(4), clf
  xx = linspace(xu(1),xu(end),30);
  plot(xx,exact_solution(0, xx, 0, 'upper'),'k-')
  hold on
  title('Velocity profiles at inflow and outflow boundaries','FontSize',16)
  vq = reshape(v1,size(XV')); vq = vq';
  plot(xv,vq(K+1,:),'o')
  plot(xv,vq(1,:),'*')
  hold off
elseif geval == 3
  figure(4), clf
  uq = reshape(u1,size(XU')); uq = uq';
  plot(uq(:,J+1),yu,'o')
  hold on
  plot(uq(:,1),yu,'*')
  hold off
else
end
