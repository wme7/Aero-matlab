% NACA 6414%
clc
clear
format compact

m1=input('Please enter one digit number for maxi camber line:')
p1=input('Please enter one digit number for distance of m and chordline:')
t1=input('Please enter two digit number for thickness of airfoil:')

m=m1*0.01;
p=p1*0.1;
t=t1*0.01;

x1=[0:0.01:p];
x2=[p:0.01:1];
x=[x1,x2]


yc1=(m*(2*p*x1-x1.^2))/(p^2);
yc2=(m/(1-p)^2)*((1-2*p)+2*p*x2-x2.^2);
yt=(t/0.2)*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4)

yc=[yc1,yc2]

ang1=(atan((m/p^2)*(2*p-2*x1)))
ang2=(atan((m/(1-p)^2)*(2*p-2*x2)))
ang=[ang1,ang2]

xu=x-yt.*sin(ang)
yu=yc+yt.*cos(ang)
xl=x+yt.*sin(ang)
yl=yc-yt.*cos(ang)

xplot=[xu,xl]
yplot=[yu,yl]

plot(xplot,yplot)