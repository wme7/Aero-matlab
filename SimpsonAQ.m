function I = SimpsonAQ(a, b, epsilon)
%% Adaptative quadrature algorithm using simpson's quadrature

I1=Simpson(a,b);
I2=Simpson2(a,b);

m=a+(b-a)/2;    %midpoint
hi=b-a;         %interval step size
if m<=a & m>=b
    fprintf('input error: a>b \n');
end
if (I2-I1)/15<=hi*epsilon
    I=I2;
else
    I=SimpsonAQ(a,m,epsilon)+SimpsonAQ(m,b,epsilon);
end