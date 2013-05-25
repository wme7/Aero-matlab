% ............................................................. 
function [weights,locations]=gauss2d(option)
    % Gauss quadrature in 2D
    % option '5x5'
    % ...
    % option '2x2'
    % option '1x1' 
    % locations: Gauss point locations
    % weights: Gauss point weights
        
%% Data Table
switch option
    case '1x1'
        x = 0;
        w = 2;
        pattern = [1];
    case '2x2'
        x = [-sqrt(3)/3,sqrt(3)/3]';
        w = [1,1]';
        pattern = [1 2 4 3];
    case '3x3'
        x = [-sqrt(3/5),0,sqrt(3/5)]';
        w = [5/9,8/9,5/9]';
        pattern = [1 3 9 7 2 6 8 4 5];
    case '4x4'
        x1 = sqrt((3-2*sqrt(6/5))/7);
        x2 = sqrt((3+2*sqrt(6/5))/7);
        x = [-x1,-x2,x2,x1]';
        w1 = (18+sqrt(30))/36;
        w2 = (18-sqrt(30))/36;
        w = [w1,w2,w2,w1]';
        pattern = [1 4 16 13 2 3 8 12 15 14 9 5 6 7 11 10];
    case '5x5'
        x1 = (1/3)*sqrt(5-2*sqrt(10/7));
        x2 = (1/3)*sqrt(5+2*sqrt(10/7));
        x3 = 0;
        x = [-x1,-x2,x3,x2,x1]';
        w1 = (322 + 13*sqrt(70))/900;
        w2 = (322 - 13*sqrt(70))/900;
        w3 = 128/225;
        w = [w1,w2,w3,w2,w1]';
        pattern = [1 5 25 21 2 3 4 10 15 20 24 23 22 16 11 6 7 8 9 14 19 18 17 12 13];
    otherwise
        error('ngp max value = 5')
end

%% Compute Arrays,
[xi1,xi2] = meshgrid(x,x);
[wi1,wi2] = meshgrid(w,w);
n = length(w)^2;
locations = [reshape(xi2,n,1),reshape(xi1,n,1)]; %locations=[x,y]
weights  = reshape(wi1.*wi2,n,1);

%% Arrange gauss points 
locations = locations(pattern,:);
weights = weights(pattern,:);

end  % end function gauss2d