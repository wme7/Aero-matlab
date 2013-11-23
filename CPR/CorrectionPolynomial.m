classdef CorrectionPolynomial
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CORRECTIONPOLYNOMIAL class
    %   Build Correction Polynomials of degree K for CPR scheme.
    %
    %   by Manuel Diaz, NTU, 2013.10.12
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        pDeg
        pType
    end
    
    properties (Dependent = true, SetAccess = private)
        xi % local coordiante
        P  % Correction Polynomail
        dP % derivate of correction Polynomial
    end
    
    methods (Static)
        function legP = LegendreP(l)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Legendre Polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : kDeg: Polynomial Degree requested
            % Output : legP: symbolic Legendre polynomial
            %
            x = sym('x'); 
            legP = simplify(1/(2^l*factorial(l))*diff((x^2-1)^l,x,l));
        end
                   
    end % Methods
           
    methods
        function obj = CorrectionPolynomial(type,Kdeg)
            obj.pDeg = Kdeg;
            obj.pType = type;
        end
        
        function cpoly = get.P(obj)
            switch obj.pType
                case 'Legendre'
                    cpoly = obj.LegendreP(obj.pDeg);
                case 'LGL'
                    cpoly = obj.LobattoP(obj.pDeg);
                case 'RadauRight'
                    cpoly = obj.RadauRightP(obj.pDeg);
                case 'RadauLeft'
                    cpoly = obj.RadauLeftP(obj.pDeg);
                otherwise
                    error('correction polynomial not available')
            end
        end
        
        function RRP = RadauRightP(obj,kDeg)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load from table the coefs for Radau Polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : kDeg: Polynomial Degree requested
            % Output : RRP: Right Radau polynomial
            %
            RRP = (-1)^(kDeg)/2*(obj.LegendreP(kDeg) - obj.LegendreP(kDeg-1));
        end
            
        function RLP = RadauLeftP(obj,kDeg)
        	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load from table the coefs for Radau Polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : kDeg: Polynomial Degree requested
            % Output : RLP: Left Radau polynomial
            %
            RLP = (1/2)*(obj.LegendreP(kDeg) + obj.LegendreP(kDeg-1) );
        end
            
        function lobP = LobattoP(obj,kDeg)
        	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Lobatto Polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : kDeg: Polynomial Degree requested
            % Output : lobP: Symbolic Lobatto polynomial
            %   
            lobP = obj.LegendreP(kDeg) - obj.LegendreP(kDeg-2);
        end
        
        function dcorrection = get.dP(obj)
            x = sym('x'); dcorrection = diff(obj.P,x);
        end
        
        function gpoints = eval_dP(obj,solutionPoints)
            gpoints = double(subs(obj.dP,solutionPoints));
        end
        
    end % Methods
    
end % Class

