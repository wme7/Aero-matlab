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
        function legP = LegendreP(kDeg)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Legendre Polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : kDeg: Polynomial Degree requested
            % Output : legP: symbolic Legendre polynomial
            %
            x = sym('x');
            switch kDeg
                case(0)
                    legP = x^0; % = 1
                case(1)
                    legP = x;
                case(2)
                    legP = 1/2*(-1+3*x^2);
                case(3)
                    legP = 1/2*(-3*x+5*x^3);
                case(4)
                    legP = 1/8*(3-30*x^2+35*x^4);
                case(5)
                    legP = 1/8*(63*x^5-70*x^3+15*x);
                case(6)
                    legP = 1/16*(231*x^6-315*x^4+105*x^2-5);
                case(7)
                    legP = 1/16*(-35*x + 315*x^3 - 693*x^5 + 429*x^7);
                case(8)
                    legP = 1/128*(-35 + 1260*x^2 - 6930*x^4 + 12012*x^6 - 6435*x^8);
                case(9)
                    legP = 1/128*(315*x - 4620*x^3 + 18018*x^5 - 25740*x^7 + 12155*x^9);
                case(10)
                    legP = 1/256*(63 - 3465*x^2 + 30030*x^4 - 90090*x^6 + 109395*x^8 - 46189*x^10);
                otherwise
                    error('Legendre Polynomial not available');
            end
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

