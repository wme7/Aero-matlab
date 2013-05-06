function [weights,locations] = gauss1d(ngp)
%GAUSS1D    Gauss quadrature in one dimension.
%         Get gauss point locations in the parent element domain [-1, 1] and the corresponding weights
%         [WEIGHTS, LOCATIONS] = GAUSS1D(NGP)
% 
switch ngp
    case 1
        locations = 0;
        weights  = 2;
    case 2
        locations = [-0.57735027, 0.57735027];
        weights  = [1,           1];
    case 3
        locations = [-0.7745966692,  0.0,           0.7745966692];
        weights = [0.5555555556,     0.8888888889,  0.5555555556];
    case 4
        locations = [-0.861136311594052, -0.339981043584856, 0.339981043584856, 0.861136311594052];
        weights = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454];
    case 5
        locations = [-0.906179845938664, -0.538469310105683, 0.0, 0.538469310105683, 0.906179845938664];
        weights = [0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189];
    otherwise
        disp('do not support more than five Gauss points')
end
