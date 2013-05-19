% ............................................................. 

    function [weights,locations]=gauss2d(option)
    % Gauss quadrature in 2D
    % option '2x2'
    % option '1x1' 
    % locations: Gauss point locations
    % weights: Gauss point weights
        
    switch option
        case '2x2'
    
        locations=...
          [ -0.577350269189626 -0.577350269189626;
             0.577350269189626 -0.577350269189626;
             0.577350269189626  0.577350269189626;
            -0.577350269189626  0.577350269189626];
        weights=[ 1;1;1;1]; 
    
        case '1x1'
        
        locations=[0 0];
        weights=[4];
    end

    end  % end function gauss2d
    