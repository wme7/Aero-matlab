function z=rpsolb(x,y,t)
% the solution of a Riemann problem (discontinuous)
  z=0.5*ones(size(x));
  ln1=3*t/8;
  ln2=-t/8;
    
  ind=find((x>ln2)&(y>ln2)&(y>-x/3)&(x>-y/3));
  z(ind)=-0.5;
  
  ind=find((y>ln1)&(x<=ln2));
  z(ind)=0.25;
  
  ind=find((x>ln1)&(y<=ln2));
  z(ind)=0.25;