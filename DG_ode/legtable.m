function [P] = legtable(x,m)
%
% input:  row vector x \subset (-1,1) of points
%         m the maximal order of the Legendre polynomials
% output: a matrix P=P(m+1,length(x)) containing the values of the
%         Legendre polynomials at the points x
%
  l=length(x);
  P=ones(m+1,l);
  if m>0
    P(2,:)=x;
    for i=2:m
      P(i+1,:)=((2*i-1)*x.*P(i,:)-(i-1)*P(i-1,:))/i;
    end;
  end
