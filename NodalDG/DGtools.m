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
        GradVadermonde
        MassMatrix
        invMassMatrix
        CoefDiffMatrix
        legLeftEnd
        legRightEnd
        nodalMassMatrix
        nodalInvMassMatrix
        nodalCoefDiffMatrix
    end
    
    methods (Static)
        function legP = legendreP(x,k)
            % Construct Legendre Polynomials
            %**************************************************************
            % Compute the value of the legendre polinomial of degree 'kDeg'
            % for every column value 'xi'
            %**************************************************************
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
        
        function h = DGnflux(f,df,u,strategy)
            %**************************************************************
            % General flux subroutine
            % We use 1 of 4 flux strategies/algorithms optimized for matlab.
            % coded by Manuel Diaz, 2012.12.21 (the end of the world)
            %**************************************************************
            
            % Parameters
            flx = f(u);     % flux value at every point of the domain
            dflx = df(u);   % flux slope at every point of the domain
            
            switch strategy
                case{1} % Roe Flux
                    
                    u_ave = (u(1,:) + u(2,:))/2;    % u @ cells boundaries
                    bool_p = df(u_ave) > 0; % if positive
                    bool_n = df(u_ave) <= 0;  % if negative
                    h = bool_p.*flx(1,:) + bool_n.*flx(2,:);
                    
                case{2} % LF
                    
                    alpha = max(max(abs(dflx)));
                    h = 0.5*(flx(1,:) + flx(2,:) - alpha*(u(2,:)- u(1,:)));
                    
                case{3} % LLF
                    
                    %fluxmat = [dflx(1:nx-1);dflx(2:nx)];
                    beta = max(abs(dflx));
                    h = 0.5*(flx(1,:) + flx(2,:) - beta.*(u(2,:)- u(1,:)));
                    
                case{4} % Upwind Flux
                    
                    % For dflux is constant along the domain!
                    a_p = max(max(dflx - abs(dflx))/2) == [0]; %#ok<NBRAK> % boolen operator for a>0
                    a_n = max(max(dflx + abs(dflx))/2) == [0]; %#ok<NBRAK> % boolen operator for a<0
                    h = a_p*flx(1,:) + a_n*flx(2,:);
                    
                otherwise
                    error('strategy not suported')
            end
        end
        
    end %Static Methods
     
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

