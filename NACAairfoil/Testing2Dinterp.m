%% 2D interpolation

[X,Y] = meshgrid(-0.5:0.01:1.5);
[vx,vy] = NacaAirfoil('0012');
dd = sqrt((X-0.5).^2+(Y-0.5).^2)-2;
%Vq = interp2(X,Y,dd,vx,vy);
%plot3(vx,vy,Vq);
%axis('equal')

[in,on] = inpolygon(X,Y,vx,vy);
airfoil = in+on;
pcolor(airfoil)