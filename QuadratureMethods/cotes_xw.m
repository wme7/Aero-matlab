function [x,w,k] = cotes_xw(a,b,N,nodes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% cotes.m
%
% 2 to 11 point Newton-Cotes nodes and weight values.
%
% Inputs: 
%
% a - lower limit of integration
% b - upper limit of integration
% N - Requested number of grid points
% nodes - number of nodes in Newton-Cotes formula
%
% Sample Usage:
%
% >>[x,w,k]cotes_wx(0,pi/2,20,5)
%
% x : Nodes in x
% w : Weight constants
% k : Quadrature constant
%
% Written by: Manuel Diaz
% Based on http://mathworld.wolfram.com/Newton-CotesFormulas.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = (nodes-1)*ceil(N/(nodes-1));  N1=N+1;
x = linspace(a,b,N1)';
h = x(2)-x(1); 

switch nodes

    case{2} % Trapezoidal Rule
    w = kron(ones(1,N/1),[2])'; w(1)=1; w(N1)=1;
    k = 1/2*h;
    
    case{3} % Simpson's Rule
    w = kron(ones(1,N/2),[2 4])'; w(1)=1; w(N1)=1;
    k = 1/3*h;
    
    case{4} % Simpson's 3/8 Rule
    w = kron(ones(1,N/3),[2 3 3])'; w(1)=1; w(N1)=1;
    k = 3/8*h;

    case{5} % Boole's 4/90 Rule
    w = kron(ones(1,N/4),[14 32 12 32])'; w(1)=7; w(N1)=7;
    k = 2/45*h;

    case{6}         
    w = kron(ones(1,N/5),[38 75 50 50 75])'; w(1)=19; w(N1)=19;
    k = 5/288*h;

    case{7} 
    w = kron(ones(1,N/6),[82 216 27 272 27 261])'; w(1)=41; w(N1)=41;
    k = 1/140*h;

    case{8} 
    w = kron(ones(1,N/7),[1502 3577 1323 2989 2989 1323 3577])'; 
    w(1)=751; w(N1)=751;
    k = 7/17280*h;

    case{9}
    w = kron(ones(1,N/8),[1978 5888 -928 10496 -4540 10496 -928 5888])'; 
    w(1)=989; w(N1)=989;
    k = 4/14175*h;

    case{10}
    w = kron(ones(1,N/9),[5714 15741 1080 19344 5778 ...
        5778 19344 1080 15741])'; 
    w(1)=2857; w(N1)=2857;
    k = 9/89600*h;

    case{11}
    w = kron(ones(1,N/10),[32134 106300 48525 272400 260550 427368 ...
        260550 272400 48525 106300])'; 
    w(1)=16067; w(N1)=16067;
    k = 5/299367*h;
        
    otherwise
        error('Order must be between 2 and 11'); 
end
return