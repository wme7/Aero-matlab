function h = cutfunc(x,a,b)
c1 = (x <= a);
c2 = (x>a & x<b);
c3 = (x >= b);
h = 0*c1+ (x-a)/(b-a).*c2 + 1*c3;
% c1 = (x <= a);
% c2 = (x>a & x<b);
% c3 = (x >= b);
% h = 0*c1+ 0.5*c2 + 1*c3;




%  for i=2:nx
% % 
%      h(1,i)= (x(1,i)-1)/(-1);
% % 
% %     h(1,1)=1;
%  end