function Q=cotes(f,a,b,N,nodes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% cotes.m
%
% 2 to 11 point Summed Newton-Cotes numerical integration formulas.
%
% Inputs: 
%
% f - function to be integrated
% a - lower limit of integration
% b - upper limit of integration
% N - Requested number of grid points
% nodes - number of nodes in Newton-Cotes formula
%
% Sample Usage:
%
% >>cotes(@sin,0,pi/2,20,5)
%
% ans =
%         0.999999999501637 
%
% Written by: Greg von Winckel
% Contact: gregvw(at)math(dot)unm(dot)edu
% URL: http://math.unm.edu/~gregvw
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N=(nodes-1)*ceil(N/(nodes-1));  N1=N+1;
x=linspace(a,b,N1)';
h=x(2)-x(1);    g=f(x);

endpts=g(1)+g(N1);

switch nodes

    case{2} % Trapezoidal Rule
    Q=(h/2)*(endpts+2*sum(g(2:N)));

    case{3} % Simpson's Rule
    Q=(h/3)*(endpts+4*sum(g(2:2:N))+2*sum(g(3:2:N)));

    case{4} % Simpson's 3/8 Rule
    Q=(3*h/8)*(endpts+3*sum(g(2:3:N)+g(3:3:N))+2*sum(g(4:3:N)));

    case{5} % Boole's 4/90 Rule
    Q=(4*h/90)*(7*endpts+32*sum(g(2:4:N))+...
        12*sum(g(3:4:N))+32*sum(g(4:4:N))+14*sum(g(5:4:N)));

    case{6}         
    Q=(5*h/288)*(19*endpts+75*sum(g(2:5:N)+g(5:5:N))+...
        50*sum(g(3:5:N)+g(4:5:N))+38*sum(g(6:5:N)));

    case{7} 
    Q=(6*h/840)*(41*endpts+216*sum(g(2:6:N)+g(6:6:N))+...
        27*sum(g(3:6:N)+g(5:6:N))+272*sum(g(4:6:N))+...
        +82*sum(g(7:6:N)));

    case{8} 
    Q=(7*h/17280)*(751*endpts+3577*sum(g(2:7:N)+g(7:7:N))+...
        1323*sum(g(3:7:N)+g(6:7:N))+2989*sum(g(4:7:N)+g(5:7:N))+...
        +1502*sum(g(8:7:N)));

    case{9}
    Q=(4*h/14175)*(989*endpts+5888*sum(g(2:8:N)+g(8:8:N))-...
       928*sum(g(3:8:N)+g(7:8:N))+10496*sum(g(4:8:N)+g(6:8:N))-...
       4540*sum(g(5:8:N))+1978*sum(g(9:8:N)));

    case{10}
    Q=(9*h/89600)*(2857*endpts+15741*sum(g(2:9:N)+g(9:9:N))+...
       1080*sum(g(3:9:N)+g(8:9:N))+19344*sum(g(4:9:N)+g(7:9:N))+...
       5778*sum(g(5:9:N)+g(6:9:N))+5714*sum(g(10:9:N)));

    case{11}
    Q=(5*h/299376)*(16067*endpts+106300*sum(g(2:10:N)+g(10:10:N))-...
       48525*sum(g(3:10:N)+g(9:10:N))+272400*sum(g(4:10:N)+g(8:10:N))-...
       260550*sum(g(5:10:N)+g(7:10:N))+427368*sum(g(6:10:N))+...
       32134*sum(g(11:10:N)));
    otherwise
        error('Order must be between 2 and 11'); 
end

