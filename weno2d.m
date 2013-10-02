%---------------------------------------
function weno2d()
%%
% Advection of a tracer by a constant streamfunction
% in a biperiodic square domain
% 
% Creation date: 10/13/2011
% email: roullet@univ-brest.fr
% 
%%
% Copyright 2011 Guillaume Roullet
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%%
param.scheme='weno5';

param.integrator='rk3';
param.flux='weno5';
param.splitting='LxF'; 
param.tracer='tiles';
%param.tracer='twopatches';
param.psi='rotation'; 
%param.psi='translation'; 

if strcmp(param.scheme,'roms')
  param.integrator='lf-am3';
  param.flux='up3';
  param.splitting='max'; 
end

if strcmp(param.scheme,'weno5')
  param.integrator='rk3'
  param.flux='weno5';
  param.splitting='LxF'; 
end

if strcmp(param.scheme,'perso')
  param.integrator='lf-am3'
  param.flux='up1';
  param.splitting='max'; 
end



t=0;
tend=1; % final time
cfl=0.8;


n=64; % resolution in x and y (square domain)
dx=1/n;

x=(1:n)*dx;
y=(1:n)*dx;

[xx,yy]=meshgrid(x,y); % T cells
[xv,yv]=meshgrid(x+dx/2,y+dx/2); % PSI points



% several streamfunctions
switch param.psi
 case 'cos'
  psi=cos(xv*2*pi).*sin(yv*2*pi)/(2*pi);
 case 'eye'
  psi=5*(yv.*(1-yv)).*(1+0.4*sin(xv*pi).^2)/(2);
 case 'line'
  psi=xv+yv;
 case 'bounded'
  psi=(xv.*(1-xv).*yv.*(1-yv));
 case 'rotation'
  psi=(xv-0.5).^2+(yv-0.5).^2;
 case 'translation'  
  psi=-yv;
end
% velocity 
u=-(psi-circshift(psi,[1 0]))/dx;
v=+(psi-circshift(psi,[0 1]))/dx;


%several tracers
switch param.tracer
 case 'twopatches'
  ra=sqrt((xx-0.38).^2+(yy-0.5).^2);
  rb=sqrt((xx-0.62).^2+(yy-0.5).^2);
  q=xx*0;
  q(ra<.1|rb<.1)=1;
 case 'tiles'
  q=mod(floor(mod(xx*6,2))+floor(mod(yy*6,2)),2);
end


% set the timestep
maxu=max(abs(u(:)));
maxv=max(abs(v(:)));
dt=dx/max(maxu,maxv)*cfl;
dt=tend/floor(tend/dt);
nt=tend/dt;

%figure(1)

param.dx=dx;

u0=u;
v0=v;

qb=q;
qbb=q;

integ=param.integrator;

%% main time loop
for kt=1:nt
%     param.integrator=integ;
%     if(kt<=2)
%         param.integrator='rk3';
%     end
    
    switch param.integrator
   case 'rk3'
    %RK3
 %   u=u0*cos(kt*pi/nt);   
 %   v=v0*cos(kt*pi/nt);
    q1 = q + dt*adv2d(q,u,v,param);
 %   u=u0*cos((kt+1)*pi/nt);   
 %   v=v0*cos((kt+1)*pi/nt);
    q2 = 0.75*q + 0.25*q1 + 0.25*dt*adv2d(q1,u,v,param);
 %   u=u0*cos((kt+1/3)*pi/nt);   
 %   v=v0*cos((kt+1/3)*pi/nt);     
    q3 = (q+2*q2+2*dt*adv2d(q2,u,v,param))/3;
    
    %swap
    qbb=qb;
    qb=q;
    q=q3;
    
   case 'ab3_explose'
    coef=5./12.;
%    u=u0*cos(kt*pi/nt);   
%    v=v0*cos(kt*pi/nt);    
    dq=adv2d(q,u,v,param)*dt;
    if(kt==1)
      dq1=dq;
      dq0=dq;
    end
    q=q+(1.5+coef)*dq-(0.5+2*coef)*dq1+coef*dq0;
    %swap
    dq0=dq1;
    dq1=dq;
    
   case 'ab3'
    coef=5./12.;
%    u=u0*cos((kt+0.5)*pi/nt);   
%    v=v0*cos((kt+0.5)*pi/nt);    
    qh=(1.5+coef)*q-(0.5+2*coef)*qb+coef*qbb;
    %swap
    qbb=qb;
    qb=q;
    q=q+adv2d(qh,u,v,param)*dt;

   case 'lf-am3'   
       gamma=1./12;
% predictor
% LF : qq is at t=n+1
%    u=u0*cos(kt*pi/nt);   
%    v=v0*cos(kt*pi/nt);    
    qq=qb+adv2d(q,u,v,param)*dt*2;
% AM3 : interpolate back to t=n+1/2   
    qh=(0.5+2*gamma)*q+(0.5-gamma)*qq-gamma*qb;
% corrector    
%    u=u0*cos((kt+0.5)*pi/nt);   
%    v=v0*cos((kt+0.5)*pi/nt);    
    %swap
    qb=q;
    q=q+adv2d(qh,u,v,param)*dt;
   
  end

  
  t  = t +dt;
 
  if(mod(kt,50)==0 || kt==nt)
    imagesc(q);axis xy;caxis([-.1 1.1])
    title(sprintf('time=%f',t))
    drawnow
  end


  
end

  disp('fini')

%---------------------------------------
function dqdt=adv2d(q,u,v,param)
% compute the RHS of dq/dt
%
% Along x
[fop,fom]=fluxsplitting(u,param);
fom=circshift(fom,[0 1]);
switch param.flux
 case 'weno5'
  dqdt =      weno5(fop,q,[0 -1]);
  dqdt = dqdt+weno5(fom,q,[0 +1]);
 case 'up3'
  dqdt =      up3(fop,q,[0 -1]);
  dqdt = dqdt+up3(fom,q,[0 +1]);
 case 'up1'
  dqdt =      up1(fop,q,[0 -1]);
  dqdt = dqdt+up1(fom,q,[0 +1]);
end
% Along y
[fop,fom]=fluxsplitting(v,param);
fom=circshift(fom,[1 0]);
switch param.flux
 case 'weno5'
  dqdt = dqdt+weno5(fop,q,[-1 0]);
  dqdt = dqdt+weno5(fom,q,[+1 0]);
 case 'up3'
  dqdt = dqdt+up3(fop,q,[-1 0]);
  dqdt = dqdt+up3(fom,q,[+1 0]);
 case 'up1'
  dqdt = dqdt+up1(fop,q,[-1 0]);
  dqdt = dqdt+up1(fom,q,[+1 0]);
end
dqdt=dqdt/param.dx;

%---------------------------------------
function [fop,fom]=fluxsplitting(u,param)
% split the flux into right-going and left-going
switch param.splitting
 case 'max'
  fop= max(u,0);
  fom=-min(u,0);
 case 'LxF'
  au=max(abs(u(:)));
  fop =  0.5*(u+au);  
  fom = -0.5*(u-au);
end

  
%---------------------------------------
function z=up1(vel,fo,ind)
% upwind 1st order
  Fp = fo.*vel;
  z = -(Fp-circshift(Fp,ind));
    
%---------------------------------------
function z=up3(vel,fo,ind)
% upwind 3rd order
  fm =circshift(fo,+ind);
  fp =circshift(fo,-ind);

  Fp = 0.5*(fp+fo)-(fm-2*fo+fp)/6;
  Fp = Fp.*vel;
  z = -(Fp-circshift(Fp,ind));
  
%---------------------------------------
function z=weno5(vel,fo,ind)
% weno 5th order
  e=1e-6;

  g1=1/10;
  g2=3/5;
  g3=3/10;
    
  fmm=circshift(fo,ind*2);
  fm =circshift(fo,ind);
  fp =circshift(fo,-ind);
  fpp=circshift(fo,-ind*2);

  fp1 = (2*fmm - 7*fm + 11*fo )/6.;
  fp2 = ( -fm  + 5*fo + 2*fp  )/6.;
  fp3 = (2*fo  + 5*fp - fpp)/6.;

  b1 = 13/12*(fmm-2*fm+fo ).^2 + 1/4*(fmm-4*fm+3*fo).^2; 
  b2 = 13/12*(fm -2*fo+fp ).^2 + 1/4*(fm-fp).^2;
  b3 = 13/12*(fo -2*fp+fpp).^2 + 1/4*(3*fo-4*fp+fpp).^2;

  wt1 = g1./(e+b1).^2;
  wt2 = g2./(e+b2).^2;
  wt3 = g3./(e+b3).^2;

  sum_wt=wt1+wt2+wt3;

  w1 = wt1./sum_wt;
  w2 = wt2./sum_wt;
  w3 = wt3./sum_wt;

  Fp = w1 .* fp1 + w2 .* fp2 + w3 .* fp3;
  Fp = Fp.*vel;
  z = -(Fp-circshift(Fp,ind));