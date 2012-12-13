function func=func0(x)



for i=1:size(x,2)

%  if (x(i)>0.4)&(x(i)<0.6)
%    func(i)=(x(i)-0.4)^10*(x(i)-0.6)^10*10^20;
%  else
%    func(i)=0;
%  end

  if (x(i)>0.4) & (x(i)<0.6)  
    func(i)=1;
  else 
    func(i)=0;
  end


% func(i)=sin(2*pi*(x(i)));

end

end