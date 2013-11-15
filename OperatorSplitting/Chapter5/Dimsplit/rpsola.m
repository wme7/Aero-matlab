function z=rpsola(x,y,t)
  % The solution of a Riemann problem (continuous)  
  z=zeros(size(x));
  lne=0.25*t;
  lm=-0.5*t;
  lM=0.5*t; 
  
  ind=find((x>=lM)&(y>=lM));
  z(ind)=0.5;
  
  ind=find((x<lM)&(x>=lne)&(y>=lne)&(y>=x));
  z(ind)=x(ind)/t;
  
  ind=find((x>=lne)&(y<lM)&(y>=lne)&(y<x));
  z(ind)=y(ind)/t;
  
  ind=find(((x>=lne)&(y<lne))|((x<lne)&(y>=lne)));
  z(ind)=0.25;
  
  ind=find((x<lne)&(x>=lm)&(x>=y));
  z(ind)=x(ind)/t;
  
  ind=find((y<lne)&(y>=lm)&(x<y));
  z(ind)=y(ind)/t;

  ind=find((x<=lm)&(y<=lm));
  z(ind)=-0.5;