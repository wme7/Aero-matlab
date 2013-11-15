function [weights,location] = GaussQuad1d(ngp)
% Guass weights and abscissas.
% Taken from http://en.wikipedia.org/wiki/Gaussian_quadrature

% Data Table
switch ngp
    case{1} % ngp =1
        location = 0;
        weights = 2;
    case{2} % ngp =2
        location = [-sqrt(3)/3,sqrt(3)/3]';
        weights = [1,1]';
    case{3} % ngp =3
        location = [-sqrt(3/5),0,sqrt(3/5)]';
        weights = [5/9,8/9,5/9]';
    case{4} % ngp =4
        x1 = sqrt((3-2*sqrt(6/5))/7);
        x2 = sqrt((3+2*sqrt(6/5))/7);
        location = [-x1,-x2,x2,x1]';
        w1 = (18+sqrt(30))/36;
        w2 = (18-sqrt(30))/36;
        weights = [w1,w2,w2,w1]';
    case{5} % ngp =5
        x1 = (1/3)*sqrt(5-2*sqrt(10/7));
        x2 = (1/3)*sqrt(5+2*sqrt(10/7));
        x3 = 0;
        location = [-x1,-x2,x3,x2,x1]';
        w1 = (322 + 13*sqrt(70))/900;
        w2 = (322 - 13*sqrt(70))/900;
        w3 = 128/225;
        weights = [w1,w2,w3,w2,w1]';
    otherwise
        error('ngp max value = 5')
end
            