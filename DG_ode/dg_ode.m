function  [err_max,err_node] = dg_ode(n,p)
pl=p+1;
u0=1;%initial condition
h=2./n;
A=zeros(p+1,p+1);
uh = zeros(n,p+1);
[xl,w]=gauleg(pl);
[P]=legtable(xl,p);
t=linspace(0,2,n+1);
for i=0:p
 for j=0:p
      if j>i& rem(j-i,2)==1
          A(i+1,j+1)=2;
      end    
 end
end  
C=ones(p+1,p+1);
for i=0:p
 for j=0:p
   sum_loc=0.;
  for k=1:pl
      tl=h/2*(xl(k)+1)+t(1);
   sum_loc=sum_loc+w(k)*fun(tl)*P(i+1,k)*P(j+1,k);
  end
  B(i+1,j+1)=(h./2.)*sum_loc;
 % f=inline('exp(((x+1).*h./2)).*leg_fun(x,i).*leg_fun(x,j)');
 % B(i+1,j+1)=(h./2).*int(f,x,-1,1);
%  b_0(j+1,1)=(-1).*(-1)^j;
%  b(j+1,1)=(-1).*series(i).*(-1)^j;
 end
 rhs(i+1)=-u0*(-1)^i;
end 
G=A+B-C;
U=G'\rhs';
sum_err=0;
for k=1:pl
tl=h/2*(xl(k)+1)+t(1);
 %U_exac(k)=u0*exp(tl);
% U_exac(k)=u0*exp(tl^2/2);
 % U_exac(k)=u0*exp(tl^3/3);
 U_exac(k)=u0*exp(exp(tl)-1);
sum_loc=0;
 for i=0:p
 sum_loc=sum_loc+U(i+1)*P(i+1,k);    
 end
 U_num(k)=sum_loc;
 sum_err=sum_err+w(k)*(U_exac(k)-U_num(k))^2;
end
err(1)=max(U_exac-U_num);
err2(1)=sum_err*h/2;
sum_loc=0;
 for i=0:p
 sum_loc=sum_loc+U(i+1);    
 end
 Un_exac=u0*exp(exp(t(1)+h)-1);
 Un_num=sum(U);
errn(1)=abs(Un_exac-Un_num);
%sum_loc);
for m=2:n
sum_loc=0;
for i=0:p
 sum_loc=sum_loc+U(i+1);    
 end
u_old=sum_loc;
 
for i=0:p
 for j=0:p
   sum_loc=0.;
  for k=1:pl
      tl=h/2*(xl(k)+1)+t(m);
   sum_loc=sum_loc+w(k)*fun(tl)*P(i+1,k)*P(j+1,k);
  end
  B(i+1,j+1)=(h./2.)*sum_loc;
end
 rhs(i+1)=-u_old*(-1)^i;
end 
G=A+B-C;
U=G'\rhs';
sum_err=0;
for k=1:pl
tl=h/2*(xl(k)+1)+t(m);
 %U_exac(k)=u0*exp(tl);
% U_exac(k)=u0*exp(tl^2/2);
  %U_exac(k)=u0*exp(tl^3/3);
 U_exac(k)=u0*exp(exp(tl)-1);
sum_loc=0;
 for i=0:p
 sum_loc=sum_loc+U(i+1)*P(i+1,k);    
 end
 U_num(k)=sum_loc;
 sum_err=sum_err+w(k)*(U_exac(k)-U_num(k))^2;
end
err(m)=max(abs(U_exac-U_num));
err2(m)=sum_err*h/2;
sum_loc=0;
 for i=0:p
 sum_loc=sum_loc+U(i+1);    
 end
 Un_exac=u0*exp(exp(t(m)+h)-1);
 Un_num=sum(U);
errn(m)=abs(Un_exac-Un_num);
%errn(m)=abs(u0*exp(exp(t(m)+h)-1)-sum(U));
%sum_loc);
end
err_max=max(err);
err2f=sqrt(sum(err2));
err_max=err2f;
err_node=max(errn);
