lx = 10;
ly = 10;

[x,y]=meshgrid(1:ly,1:lx);

obst_x = lx/5+1;
obst_y = ly/2+1;
obst_r = ly/10+1;

obst = (x-obst_x).^2 + (y-obst_y).^2 <= obst_r.^2;
obst(:,[1,ly]) = 1;
bbRegion = find(obst);

count = 0;
for nn = 1:10
    if (mod(nn,2)==1);
        count = count + 1;
        plot(x,rho);
        M(count)=getframe;
    end
end