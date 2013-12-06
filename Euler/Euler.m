%----------------------_-----------------------
% Riemann Solver for solving shock-tube problem
% by Manuel Diaz
%-----------------------_----------------------
% Institute of Applied Mechanics
% Aerodynamic Design and Analsis Lab
% Laboratory 007
% Date 03/09/2011
%------------------------_---------------------
% NOTE:
% 6 Possible cases of formation of shock waves 
% and expansion fan have been considered, in-
% cluding the case of cavitation.
% A Cavitation Check is the is incorporated in
% the code. It further prevents plotting for
% possible but physically unlikely case of ex-
% Expansion shocks.

clear;
clc;
close all;

% Ratio of specific heats for air
gamma = 1.4;
% Problem definition: Conditions at time t=0
% Case 1 : Sod's Problem
ch=0;
while(ch==0)
fprintf ('Choose one of the following cases :- \n');
fprintf ('\n \t Case 1: Sods problem \n');
fprintf ('\t Case 2: Left running expansion and rightrunning "STRONG" shock \n');
fprintf ('\t Case 3: Left running shock and right runningexpansion \n');
fprintf ('\t Case 4: Double shock \n');
fprintf ('\t Case 5: Double expansion \n');
fprintf ('\t Case 6: Cavitation \n');
cas=input ('\nEnter a case no. <1-6>: ');
if cas==1
% Case 1:Left Expansion & right Shock
fprintf('Case 1:Sods problem \n');
rho1=1;
rho4=0.125;
u1=0;
u4=0;
p1=1;
p4=0.1;
fprintf ('P1 = %f \n',p1);
fprintf ('P4 = %f \n',p4);
fprintf ('U1 = %f \n',u1);
fprintf ('U4 = %f \n',u4);
fprintf ('rho1 = %f \n',rho1);
fprintf ('rho4 = %f \n',rho4);
ch=1;
elseif cas==2
% Case 2:Strong Expansion & Shock
fprintf('Case 2:Strong ...Expansion & Shock \n');
rho1=3;
rho4=2;
u1=0;
u4=0;
p1=1000;
p4=0.01;
fprintf ('P1 = %f \n',p1);
fprintf ('P4 = %f \n',p4);
fprintf ('U1 = %f \n',u1);
fprintf ('U4 = %f \n',u4);
fprintf ('rho1 = %f \n',rho1);
fprintf ('rho4 = %f \n',rho4);
ch=1;
elseif cas==3
% Case 3:Shock & Expansion
fprintf('Case 3:Shock & Expansion \n');
rho1=1;
rho4=1;
u1=0;
u4=0;
p1=7;
p4=10;
fprintf ('P1 = %f \n',p1);
fprintf ('P4 = %f \n',p4);
fprintf ('U1 = %f \n',u1);
fprintf ('U4 = %f \n',u4);
fprintf ('rho1 = %f \n',rho1);
fprintf ('rho4 = %f \n',rho4);
ch=1;
elseif cas==4
% Case 4:Double Shock
fprintf('Case 4:Double Shock \n');
rho1=6;
rho4=6;
u1=20;
u4=-6;
p1=450;
p4=45;
fprintf ('P1 = %f \n',p1);
fprintf ('P4 = %f \n',p4);
fprintf ('U1 = %f \n',u1);
fprintf ('U4 = %f \n',u4);
fprintf ('rho1 = %f \n',rho1);
fprintf ('rho4 = %f \n',rho4);
ch=1;
elseif cas==5
% Case 5:Double Expansion
fprintf('Case 5:Double Expansion \n');
rho1=1;
rho4=2.5;
u1=-2;
u4=2;
p1=40;
p4=40;
fprintf ('P1 = %f \n',p1);
fprintf ('P4 = %f \n',p4);
fprintf ('U1 = %f \n',u1);
fprintf ('U4 = %f \n',u4);
fprintf ('rho1 = %f \n',rho1);
fprintf ('rho4 = %f \n',rho4);
ch=1;
elseif cas==6
% Case 6:Cavitation
fprintf('Case 6:Cavitation \n');
rho1=1;
rho4=1;
u1=-20;
u4=20;
p1=0.40;
p4=0.40;
fprintf ('P1 = %f \n',p1);
fprintf ('P4 = %f \n',p4);
fprintf ('U1 = %f \n',u1);
fprintf ('U4 = %f \n',u4);
fprintf ('rho1 = %f \n',rho1);
fprintf ('rho4 = %f \n',rho4);
ch=1;
else
fprintf ('Please enter an appropriate choice \n');
end % for case selection if-else loop
end % for case selection while loop
% Calculation of flow parameters at initial condition
a1=sqrt(gamma*p1/rho1);
a4=sqrt(gamma*p4/rho4);
M1=u1/a1;
M4=u4/a4;
if u1<0 && u4>0 && (u1+(2/(gamma-1))*a1)<=(u4-(2/(gamma-1))*a4)
%chk for cavitation
fprintf('\n Cavitation is observed in this case of double expansion. No solution possible \n');
break;
else
% Secant Method for getting pressure P*
p23up=(((gamma-1)/2*(u1-u4)+a1+a4)/((a1*(p1)^((2*gamma)/(gamma-1)))+(a4*(p4)^((2*gamma)/(gamma-1)))))^((2*gamma)/(gamma-1)); % upper limit
p23down=(rho1*a1*p4+rho4*a4*p1-(rho1*a1*rho4*a4*(u4-u1)))/(rho1*a1+rho4*a4); %lower limit by linear theory
s=0;
if p23down>=p1
s=1;
end
ss=0;
if p23down>=p4
ss=1;
end
if s==1
m1=rho1*a1*sqrt(1+((gamma+1)*(p23up-p1)/(2*gamma*p1)));
m1d=rho1*a1*sqrt(1+((gamma+1)*(p23downp1)/(2*gamma*p1)));
else
m1=rho1*a1*(gamma-1)/(2*gamma)*(1-p23up/p1)/(1-(p23up/p1)^((gamma-1)/(2*gamma)));
m1d=rho1*a1*(gamma-1)/(2*gamma)*(1-p23down/p1)/(1-(p23down/p1)^((gamma-1)/(2*gamma)));
end
if ss==1
m4=rho4*a4*sqrt(1+((gamma+1)*(p23up-p4)/(2*gamma*p4)));
m4d=rho4*a4*sqrt(1+((gamma+1)*(p23down-p4)/(2*gamma*p4)));
else
m4=rho4*a4*(gamma-1)/(2*gamma)*(1-p23up/p4)/(1-(p23up/p4)^((gamma-1)/(2*gamma)));
m4d=rho4*a4*(gamma-1)/(2*gamma)*(1-p23down/p4)/(1-(p23down/p4)^((gamma-1)/(2*gamma)));
end
p23=(m1*p4+m4*p1-m1*m4*(u4-u1))/(m1+m4);
f=p23up-p23;
p23d=(m1d*p4+m4d*p1-m1d*m4d*(u4-u1))/(m1d+m4d);
ff=p23d-p23;
j=0;
% iteration procedure starts from here
while abs(f)>0.000001
p23up=p23up-(f*(p23up-p23d)/(f-ff));
if s==1
m1=rho1*a1*sqrt(1+((gamma+1)*(p23upp1)/(2*gamma*p1)));
m1d=rho1*a1*sqrt(1+((gamma+1)*(p23downp1)/(2*gamma*p1)));
else
m1=rho1*a1*(gamma-1)/(2*gamma)*(1-p23up/p1)/(1-(p23up/p1)^((gamma-1)/(2*gamma)));
m1d=rho1*a1*(gamma-1)/(2*gamma)*(1-p23down/p1)/(1-(p23down/p1)^((gamma-1)/(2*gamma)));
end
if ss==1
m4=rho4*a4*sqrt(1+((gamma+1)*(p23up-p4)/(2*gamma*p4)));
m4d=rho4*a4*sqrt(1+((gamma+1)*(p23down-p4)/(2*gamma*p4)));
else
m4=rho4*a4*(gamma-1)/(2*gamma)*(1-p23up/p4)/(1-(p23up/p4)^((gamma-1)/(2*gamma)));
m4d=rho4*a4*(gamma-1)/(2*gamma)*(1-p23down/p4)/(1-(p23down/p4)^((gamma-1)/(2*gamma)));
end
p23=(m1*p4+m4*p1-m1*m4*(u4-u1))/(m1+m4);
f=p23up-p23;
p23d=(m1d*p4+m4d*p1-m1d*m4d*(u4-u1))/(m1d+m4d);
ff=p23d-p23;
j=j+1;
if j>450000;
fprintf ('No convergance \n');
break;
end
end % for while loop of secant method
% Root finder ends
%Calculation of flow parameters depending whether shock or
% an
%expansion is observed
u23=(m1*u1+m4*u4-(p4-p1))/(m1+m4);
if s==1
rho2=rho1*(1+(((gamma+1)/(gamma-1))*p23/p1))/(((gamma+1)/(gamma-1))+p23/p1);
if u23>u1
fprintf('Expansion shock-not physically possible\n');
return;
end
else
rho2=rho1*(p23/p1)^(1/gamma);
end
if ss==1
rho3=rho4*(1+(((gamma+1)/(gamma-1))*p23/p4))/(((gamma+1)/(gamma-1))+p23/p4);
if u23<u4
fprintf('Expansion shock-not physically possible\n');
return;
end
else
rho3=rho4*(p23/p4)^(1/gamma);
end
a2=sqrt(gamma*p23/rho2);
a3=sqrt(gamma*p23/rho3);
end % for cavitation test
% Print calculated flow quantities
fprintf ('\n Solution of Riemann problem :- \n');
fprintf ('P* = %f \n',p23);
fprintf ('U* = %f \n',u23);
fprintf ('rho2 = %f \n',rho2);
fprintf ('rho3 = %f \n',rho3);
% Variable initialization
cs12=0;
cs34=0;
expc121=0;
expc122=0;
expc341=0;
expc342=0;
% Calculation of shock/expansion speeds
if s==1
cs12l=a1*sqrt(1+(gamma+1)/(2*gamma)*(p23-p1)/p1);
cs12=u1-abs(cs12l);
else
expc121=u1-a1;
expc122=u23-a2;
end
if ss==1
cs34r=a4*sqrt(1+(gamma+1)/(2*gamma)*(p23-p4)/p4);
cs34=u4+abs(cs34r);
else
expc341=u4+a4;
expc342=u23+a3;
end
%Array construction
maxxt=max([cs12 cs34 expc121 expc122 expc341 expc342]);
minxt=min([cs12 cs34 expc121 expc122 expc341 expc342]);
offsetxt=0.1*(maxxt-minxt);
if s==1
xt(1)=cs12-offsetxt;
incr=abs(offsetxt)/1500;
for i=1:1500
xt(i+1)=xt(i)+incr;
u(i)=u1;
rho(i)=rho1;
p(i)=p1;
e(i)=p(i)/(gamma-1)/rho(i);
end
xt(1500)=cs12;
incr=abs(u23-cs12)/1500;
for i=1501:3000
xt(i+1)=xt(i)+incr;
u(i)=u23;
rho(i)=rho2;
p(i)=p23;
e(i)=p(i)/(gamma-1)/rho(i);
end
xt(3000)=u23;
else
xt(1)=expc121-offsetxt;
incr=abs(offsetxt)/1000;
for i=1:1000
xt(i+1)=xt(i)+incr;
u(i)=u1;
rho(i)=rho1;
p(i)=p1;
e(i)=p(i)/(gamma-1)/rho(i);
end
xt(1000)=expc121;
incr=abs(expc122-expc121)/1000;
for i=1001:2000
xt(i+1)=xt(i)+incr;
if expc122>=0
u(i)=2/(gamma+1)*(xt(i)-a1)+(gamma-1)/(gamma+1)*u1;
a=((gamma-1)/(gamma+1)*(xt(i)-u1))+(2/(gamma+1)*a1);
else
u(i)=2/(gamma+1)*(xt(i)+a1)+(gamma-1)/(gamma+1)*u1;
a=(-(gamma-1)/(gamma+1)*(xt(i)-u1))+(2/(gamma+1)*a1);
end
rho(i)=rho1*(a/a1)^(2/(gamma-1));
p(i)=p1*(a/a1)^(2*gamma/(gamma-1));
e(i)=p(i)/(gamma-1)/rho(i);
end
xt(2000)=expc122;
incr=abs(expc122-u23)/1000;
for i=2001:3000
xt(i+1)=xt(i)+incr;
u(i)=u23;
rho(i)=rho2;
p(i)=p23;
e(i)=p(i)/(gamma-1)/rho(i);
end
xt(3000)=u23;
end
if ss==1
incr=abs(u23-cs34)/1500;
for i=3001:4500
xt(i)=xt(i-1)+incr;
u(i)=u23;
rho(i)=rho3;
p(i)=p23;
e(i)=p(i)/(gamma-1)/rho(i);
end
xt(4500)=cs34;
incr=abs(offsetxt)/1500;
for i=4501:6000
xt(i)=xt(i-1)+incr;
u(i)=u4;
rho(i)=rho4;
p(i)=p4;
e(i)=p(i)/(gamma-1)/rho(i);
end
else
incr=abs(expc342-u23)/1000;
for i=3001:4000
xt(i)=xt(i-1)+incr;
u(i)=u23;
rho(i)=rho3;
p(i)=p23;
e(i)=p(i)/(gamma-1)/rho(i);
end
xt(4000)=expc342;
incr=abs(expc342-expc341)/1000;
for i=4001:5000
xt(i)=xt(i-1)+incr;
if expc341>=0
u(i)=2/(gamma+1)*(xt(i)-a4)+(gamma-1)/(gamma+1)*u4;
a=((gamma-1)/(gamma+1)*(xt(i)-u4))+(2/(gamma+1)*a4);
else
u(i)=2/(gamma+1)*(xt(i)+a4)+(gamma-1)/(gamma+1)*u4;
a=(-(gamma-1)/(gamma+1)*(xt(i)-u4))+(2/(gamma+1)*a4);
end
rho(i)=rho4*(a/a4)^(2/(gamma-1));
p(i)=p4*(a/a4)^(2*gamma/(gamma-1));
e(i)=p(i)/(gamma-1)/rho(i);
end
xt(5000)=expc341;
incr=abs(offsetxt)/1000;
for i=5001:6000
xt(i)=xt(i-1)+incr;
u(i)=u4;
rho(i)=rho4;
p(i)=p4;
e(i)=p(i)/(gamma-1)/rho(i);
end
end
% Plotting instructions
subplot (2,2,1)
plot(xt,u);
title('Plot of U v/s x/t');
xlabel ('x/t');
ylabel ('u');
axis tight;
subplot (2,2,2)
plot(xt,rho);
title('Plot of Density v/s x/t');
xlabel ('x/t');
ylabel ('rho');
axis tight;
subplot (2,2,3)
plot(xt,p);
title('Plot of Pressure v/s x/t');
xlabel ('x/t');
ylabel ('P');
axis tight;
subplot (2,2,4)
plot(xt,e);
title('Plot of Internal Energy v/s x/t');
xlabel ('x/t');
ylabel ('E');
axis tight;
% ---------END OF CODE----------- %