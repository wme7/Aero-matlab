classdef DGtools
    %DGTOOLS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        kDeg
        nNodes
        xi
    end
    
    properties (Dependent = true, SetAccess = private)
        Vadermonde
        MassMatrix
        invMassMatrix
        CoefDiffMatrix
        legLeftEnd
        legRightEnd
    end
    
    methods (Static)
        function legP = legendreP(x,k)
            % Construct Legendre Polynomials
            %**************************************************************
            % Compute the value of the legendre polinomial of degree 'kDeg'
            % for every column value 'xi'
            %**************************************************************
            
            %number of points
            p = k+1;
            switch k
                case {0} % when k = 0
                    legP = ones(size(x));
                case {1} % when k = 1
                    legP = x;
                otherwise % k >= 2
                    legP0 = ones(size(x));
                    legP1 = x;
                    for m = 1:k-1
                        legPX = ((2*m+1)*x.*legP1 - m*legP0)./(m+1);
                        legP0 = legP1;
                        legP1 = legPX;
                    end
                    legP = legPX;
            end
        end
    end
    
    methods
        function obj = DGtools(localNodes) % Constructor
            obj.nNodes = length(localNodes);
            obj.kDeg = length(localNodes)-1;
            obj.xi = localNodes;
        end
        
        function V = get.Vadermonde(obj)
            % Construct Legendre Vandermonde Matrix, V,
            %**************************************************************
            % V is a matrix array in which each element is a Legendre
            % Polynomial of Degree n, n=0:kDeg evaluated at a column array
            % of points, xi.
            %
            % Coded by Manuel Diaz 2012.12.05
            %**************************************************************
            V = zeros(obj.kDeg);% Allocate
            for l = 0:obj.kDeg 	% All Polynomial Degrees up to kDeg
                j = l + 1;      % Dummy index
                for i = 1:obj.nNodes  % All Element Points in xi
                    V(i,j) = obj.legendreP(obj.xi(i),l);
                end
            end
        end
        
        function legL = get.legLeftEnd(obj)
            % Scaled Legendre polynomials of degree 'l', l=0:kDeg,
            % evaluated at x = -1 
            %**************************************************************
            % same as evaluating legendreP(-1,l)
            %**************************************************************
            l = 0:obj.kDeg;
            legL = (-1).^(l); % LegP @ x_{i-1/2}^(+)
        end
    
        function legR = get.legRightEnd(obj)
            % Scaled Legendre polynomials of degree 'l', l=0:kDeg,
            % evaluated at x = +1 
            %**************************************************************
            % same as evaluating legendreP(+1,l)
            %**************************************************************
            l = 0:obj.kDeg;
            legR = (1).^(l);  % LegP @ x_{i+1/2}^(-)
        end
        
        function M = get.MassMatrix(obj)
            %**************************************************************
            % Mass Matrix
            %**************************************************************
            l = 0:obj.kDeg;
            M = diag(2./(2*l+1));
        end
        
        function invM = get.invMassMatrix(obj)
            %**************************************************************
            % Inverse of Mass Matrix
            %**************************************************************
            invM = inv(obj.MassMatrix);
        end
        
        function Dhat = get.CoefDiffMatrix(obj)
            %**************************************************************
            % Coefficients Differenciation Matrix
            %**************************************************************
            k = obj.kDeg; Dhat = zeros(k+1);
            for m = 0:k
                for j = m+1:2:k
                    Dhat(m+1,j+1) = (2*m+1);
                end
            end
        end
                
        % transform u(x,t) to degress of freedom u(t)_{l,i} for each i-Cell/Element
        % ut = zeros(np,nx);
        % for l = 0:k             % for all degress of freedom
        %     i = l+1;            % Dummy index
        %     for j = 1:nx
        %         ut(i,j) = (2*l+1)/2*sum(w(:,j).*u(:,j).*V(:,i));
        %     end
        % end
        %
        % % transform f(x,t) to degress of freedom f(t)_{l,i} for each i-Cell/Element
        % ft = zeros(np,nx);
        % for l = 0:k             % for all degress of freedom
        %     i = l+1;            % Dummy index
        %     for j = 1:nx
        %         ft(i,j) = (2*l+1)/2*sum(w(:,j).*f(:,j).*V(:,i));
        %     end
        % end
        %
        % % transform s(x,t) to degress of freedom s(t)_{l,i} for each i-Cell/Element
        % st = zeros(np,nx);
        % for l = 0:k             % for all degress of freedom
        %     i = l+1;            % Dummy index
        %     for j = 1:nx
        %         %st(i,j) = (2*l+1)/2.*sum(w(:,j).*s(:,j).*V(:,i));
        %         st(i,j) = dx/2*sum(w(:,j).*s(:,j).*V(:,i));
        %     end
        % end
        %
        % % Volume term (test on methods)
        % % 1. Cell Integral volume by Gauss-type quadrature.
        % v_term1 = zeros(np,nx);
        % for l = 0:k
        %     i = l+1;    % Dummy index
        %     for j = 1:nx
        %         v_term1(i,j) = (2*l+1)/2*sum(w(:,j).*f(:,j).*dsLegendreP(l,xi));
        %     end
        % end
        % % 2. Cell Integral volume by Prof T.W. method
        % v_term2 = D'*ft;
        
    end % Methods
end % Class

