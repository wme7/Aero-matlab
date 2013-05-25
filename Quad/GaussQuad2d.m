function [weights,location] = GaussQuad2d(option)
% Guass weights and abscissas.
% Taken from http://en.wikipedia.org/wiki/Gaussian_quadrature

%% Data Table
switch option
    case '1x1'
        x = 0;
        w = 2;
    case '2x2'
        x = [-sqrt(3)/3,sqrt(3)/3]';
        w = [1,1]';
    case '3x3'
        x = [-sqrt(3/5),0,sqrt(3/5)]';
        w = [5/9,8/9,5/9]';
    case '4x4'
        x1 = sqrt((3-2*sqrt(6/5))/7);
        x2 = sqrt((3+2*sqrt(6/5))/7);
        x = [-x1,-x2,x2,x1]';
        w1 = (18+sqrt(30))/36;
        w2 = (18-sqrt(30))/36;
        w = [w1,w2,w2,w1]';
    case '5x5'
        x1 = (1/3)*sqrt(5-2*sqrt(10/7));
        x2 = (1/3)*sqrt(5+2*sqrt(10/7));
        x3 = 0;
        x = [-x1,-x2,x3,x2,x1]';
        w1 = (322 + 13*sqrt(70))/900;
        w2 = (322 - 13*sqrt(70))/900;
        w3 = 128/225;
        w = [w1,w2,w3,w2,w1]';
    otherwise
        error('ngp max value = 5')
end

%% Compute Arrays,
[xi1,xi2] = meshgrid(x,x);
[wi1,wi2] = meshgrid(w,w);
n = length(w)^2;
location = [reshape(xi1,n,1),reshape(xi2,n,1)];
weights  = reshape(wi1.*wi2,n,1);
            