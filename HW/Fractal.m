col=30;
m=400;
cx=0;
cy=0;
l=1.5;
x=linspace(cx-l,cx+l,m);
y=linspace(cy-l,cy+l,m);
[X,Y]=meshgrid(x,y);
c= -.745429;
Z=X+1i*Y;
for k=1:col;
Z=Z.^2+c;
W=exp(-abs(Z));
end
colormap prism(256)
pcolor(W);
shading flat;
axis('square','equal','off');