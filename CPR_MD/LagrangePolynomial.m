classdef LagrangePolynomial
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LAGRANGEPOLYNOMIAL class
    %   Build Lagrange interpolation polynomials
    %
    %   by Manuel Diaz, NTU, 2013.10.12
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        x0
        Pdeg
        nPoints
    end
    
    properties(Dependent = true, SetAccess = private)
       lagrangePolynomial
       dlagrangePolynomial
    end
    
    methods
        function obj = LagrangePolynomial(x0)
           obj.x0 = x0;
           obj.Pdeg = length(x0)-1;
           obj.nPoints = length(x0);
        end
        
        function l = get.lagrangePolynomial(obj)
            % Build Lagrange base functions
            syms x;
            for i=1:obj.nPoints
                l(i)=x/x;
                for j=1:obj.nPoints
                    if(i ~= j)
                        l(i)=l(i)*(x-obj.x0(j))/(obj.x0(i)-obj.x0(j));
                    end
                end
            end
        end
        
        function D = get.dlagrangePolynomial(obj)
            % Build derivate of Lagrange base functions
            syms x;
            for i=1:obj.nPoints
                l(i)=x/x;
                for j=1:obj.nPoints
                    if(i ~= j)
                        l(i)=l(i)*(x-obj.x0(j))/(obj.x0(i)-obj.x0(j));
                    end
                end
                D(i)=diff(l(i),x);
            end
        end
        
        function Dn = dnlagrangePolynomial(obj,n)
            % Build derivate of Lagrange base functions
            syms x;
            for i=1:obj.nPoints
                l(i)=x/x;
                for j=1:obj.nPoints
                    if(i ~= j)
                        l(i)=l(i)*(x-obj.x0(j))/(obj.x0(i)-obj.x0(j));
                    end
                end
                Dn(i)=diff(l(i),x,n);
            end
        end
        
    end % Method
end % Class
