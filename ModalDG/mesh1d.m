classdef mesh1d
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %MESH1D class
    %   Build a one-dimensional mesh for polynomial reconstructions
    %
    %   by Manuel Diaz, NTU, 2013.10.12
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        range
        quadratureType
        solutionDegree
        nElements
        nFaces
        nSolutionPoints
    end
    
    properties (Dependent = true, SetAccess = private)
        nNodes
        solutionPoints
        weights
        elementNodes
        elementFaces
        elementCenter
        elementSize
        nodeCoordinates
        Jacobian
    end
    
    methods (Static)
        function [x,w] = GaussLaguerre(n, alpha)
            % This function determines the abscissas (x) and weights (w) for the
            % Gauss-Laguerre quadrature of order n>1, on the interval [0, +infinity].
            % Unlike the function 'GaussLaguerre', this function is valid for
            % n>=34. This is due to the fact that the companion matrix (of the n'th
            % degree Laguerre polynomial) is now constructed as a symmetrical
            % matrix, guaranteeing that all the eigenvalues (roots) will be real.

            % � Geert Van Damme
            % geert@vandamme-iliano.be
            % February 21, 2010

            % Building the companion matrix CM
            % CM is such that det(xI-CM)=L_n(x), with L_n the Laguerree polynomial
            % under consideration. Moreover, CM will be constructed in such a way
            % that it is symmetrical.
            i   = 1:n;
            a   = (2*i-1) + alpha;
            b   = sqrt( i(1:n-1) .* ((1:n-1) + alpha) );
            CM  = diag(a) + diag(b,1) + diag(b,-1);

            % Determining the abscissas (x) and weights (w)
            % - since det(xI-CM)=L_n(x), the abscissas are the roots of the
            %   characteristic polynomial, i.d. the eigenvalues of CM;
            % - the weights can be derived from the corresponding eigenvectors.
            [V,L]   = eig(CM);
            [x,ind] = sort(diag(L));
            V       = V(:,ind)';
            w       = gamma(alpha+1) .* V(:,1).^2;
        end
        
        function [x,w,P]=GaussLobatto(kk)
        	N = kk-1; % Compute for the number of points kk
            % lglnodes.m
            %
            % Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde
            % matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
            % integration and spectral methods.
            %
            % Reference on LGL nodes and weights:
            %   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
            %   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
            %
            % Written by Greg von Winckel - 04/17/2004
            % Contact: gregvw@chtm.unm.edu
            
            % Truncation + 1
            N1=N+1;
            % Use the Chebyshev-Gauss-Lobatto nodes as the first guess
            x=-cos(pi*(0:N)/N)';
            % The Legendre Vandermonde Matrix
            P=zeros(N1,N1);
            % Compute P_(N) using the recursion relation
            % Compute its first and second derivatives and
            % update x using the Newton-Raphson method.
            xold=2;
            while max(abs(x-xold))>eps
                xold=x;
                P(:,1)=1;    P(:,2)=x;
                for k=2:N
                    P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
                end
                x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
            end
            w=2./(N*N1*P(:,N1).^2);
        end
        
        function [x,w] = GaussLegendre(n)
            % This function determines the abscisas (x) and weights (w)  for the        %
            % Gauss-Legendre quadrature, of order n>1, on the interval [-1, +1].        %
            %   Unlike many publicly available functions, 'GaussLegendre_2' is valid    %
            %   for n>=46. This is due to the fact that 'GaussLegendre_2' does not      %
            %   rely on the build-in Matlab routine 'roots' to determine the roots of   %
            %   the Legendre polynomial, but finds the roots by looking for the         %
            %   eigenvalues of an alternative version of the companion matrix of the    %
            %   n'th degree Legendre polynomial. The companion matrix is constructed    %
            %   as a symmetrical matrix, guaranteeing that all the eigenvalues          %
            %   (roots) will be real. On the contrary, the 'roots' function uses a      %
            %   general form for the companion matrix, which becomes unstable at        %
            %   higher orders n, leading to complex roots.                              %
            
            % Geert Van Damme
            % geert@vandamme-iliano.be
            % February 21, 2010
            
            % Building the companion matrix CM
            % CM is such that det(xI-CM)=P_n(x), with P_n the Legendre polynomial
            % under consideration. Moreover, CM will be constructed in such a way
            % that it is symmetrical.
            i  = 1:n-1;
            a  = i./sqrt(4*i.^2-1);
            CM = diag(a,1) + diag(a,-1);
            
            % Determining the abscissas (x) and weights (w)
            % - since det(xI-CM)=P_n(x), the abscissas are the roots of the
            %   characteristic polynomial, i.d. the eigenvalues of CM;
            % - the weights can be derived from the corresponding eigenvectors.
            [V,L] = eig(CM);
            [x,ind] = sort(diag(L));
            % V = V'
            w = 2 * (V(1,:).^2)';
            clear('ind');
        end
        
        function [x,w,P] = GaussRadau(kk)
            % Computes the Legendre-Gauss-Radau nodes, weights and the LGR Vandermonde
            % matrix. The LGR nodes are the zeros of P_N(x)+P_{N+1}(x).
          
            % References on LGR nodes and weights:
            %   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
            %   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
            %
            %   F. B. Hildebrand , "Introduction to Numerical Analysis," Section 8.11
            %   Dover 1987
            %
            % Written by Greg von Winckel - 05/02/2004
            % Contact: gregvw@chtm.unm.edu
            
            % Compute for the number of points kk
            N = kk-1; 
            % Truncation + 1
            N1=N+1;
            % Use Chebyshev-Gauss-Radau nodes as initial guess for LGR nodes
            x=-cos(2*pi*(0:N)/(2*N+1))';
            % The Legendre Vandermonde Matrix
            P=zeros(N1,N1+1);
            % Compute P_(N) using the recursion relation
            % Compute its first and second derivatives and
            % update x using the Newton-Raphson method.
            xold=2;
            % Free abscissae
            free=2:N1;
            while max(abs(x-xold))>eps
                xold=x;
                P(1,:)=(-1).^(0:N1);
                P(free,1)=1;    P(free,2)=x(free);
                for k=2:N1
                    P(free,k+1)=( (2*k-1)*x(free).*P(free,k)-(k-1)*P(free,k-1) )/k;
                end
                x(free)=xold(free)-((1-xold(free))/N1).*(P(free,N1)+P(free,N1+1))...
                    ./(P(free,N1)-P(free,N1+1));
            end
            % The Legendre-Gauss-Radau Vandermonde
            P=P(1:N1,1:N1);
            % Compute the weights
            w=zeros(N1,1);
            w(1)=2/N1^2;
            w(free)=(1-x(free))./(N1*P(free,N1)).^2;
        end
        
        function [x,w] = GaussChebyshev(n)
            j = 1:n;
            x = cos((2*j-1)*pi/(2*n));
            w = pi/n*ones(size(x));
        end
        
        function x = Chebyshev(n)
            x = sin(0.5*pi*linspace(-1,1,n)');
        end
        
        function [x,w,k] = cotes_xw(a,b,N,nodes)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % cotes.m
            %
            % 2 to 11 point Newton-Cotes nodes and weight values.
            %
            % Inputs:
            %
            % a - lower limit of integration
            % b - upper limit of integration
            % N - Requested number of grid points
            % nodes - number of nodes in Newton-Cotes formula
            %
            % Sample Usage:
            %
            % >>[x,w,k]cotes_wx(0,pi/2,20,5)
            %
            % x : Nodes in x
            % w : Weight constants
            % k : Quadrature constant
            %
            % Written by: Manuel Diaz, NTU, 2012
            % Based on http://mathworld.wolfram.com/Newton-CotesFormulas.html
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N = (nodes-1)*ceil(N/(nodes-1));  N1=N+1;
            x = linspace(a,b,N1)';
            h = x(2)-x(1);
            
            switch nodes
                
                case{2} % Trapezoidal Rule
                    w = kron(ones(1,N/1),2)'; w(1)=1; w(N1)=1;
                    k = 1/2*h;
                    
                case{3} % Simpson's Rule
                    w = kron(ones(1,N/2),[2 4])'; w(1)=1; w(N1)=1;
                    k = 1/3*h;
                    
                case{4} % Simpson's 3/8 Rule
                    w = kron(ones(1,N/3),[2 3 3])'; w(1)=1; w(N1)=1;
                    k = 3/8*h;
                    
                case{5} % Boole's 4/90 Rule
                    w = kron(ones(1,N/4),[14 32 12 32])'; w(1)=7; w(N1)=7;
                    k = 2/45*h;
                    
                case{6}
                    w = kron(ones(1,N/5),[38 75 50 50 75])'; w(1)=19; w(N1)=19;
                    k = 5/288*h;
                    
                case{7}
                    w = kron(ones(1,N/6),[82 216 27 272 27 261])'; w(1)=41; w(N1)=41;
                    k = 1/140*h;
                    
                case{8}
                    w = kron(ones(1,N/7),[1502 3577 1323 2989 2989 1323 3577])';
                    w(1)=751; w(N1)=751;
                    k = 7/17280*h;
                    
                case{9}
                    w = kron(ones(1,N/8),[1978 5888 -928 10496 -4540 10496 -928 5888])';
                    w(1)=989; w(N1)=989;
                    k = 4/14175*h;
                    
                case{10}
                    w = kron(ones(1,N/9),[5714 15741 1080 19344 5778 ...
                        5778 19344 1080 15741])';
                    w(1)=2857; w(N1)=2857;
                    k = 9/89600*h;
                    
                case{11}
                    w = kron(ones(1,N/10),[32134 106300 48525 272400 260550 427368 ...
                        260550 272400 48525 106300])';
                    w(1)=16067; w(N1)=16067;
                    k = 5/299367*h;
                    
                otherwise
                    error('Order must be between 2 and 11');
            end
        end
        
    end % Static Methods
    
    methods
        function obj = mesh1d(range,nE,type,KDeg) % The Constuctor
            if nargin > 0 % Support calling with 0 arguments
                obj.range = range;
                obj.quadratureType = type;
                obj.nElements = nE;
                obj.nFaces    = nE+1;
                obj.solutionDegree = KDeg;
                obj.nSolutionPoints = KDeg+1;
            end
        end
        
        function SPs = get.solutionPoints(obj)
            switch obj.quadratureType
                case 'Legendre'
                    SPs = obj.GaussLegendre(obj.nSolutionPoints);
                case 'Laguerre'
                    SPs = obj.GaussLaguerre(obj.nSolutionPoints);
                case 'LGL'
                    SPs = obj.GaussLobatto(obj.nSolutionPoints);
                case 'Radau'
                    SPs = obj.GaussRadau(obj.nSolutionPoints);
                case 'Chebyshev'
                    SPs = obj.GaussChebyshev(obj.nSolutionPoints);
                case 'ChebyshevMod'
                    SPs = obj.Chebyshev(obj.nSolutionPoints);
                otherwise
                    error('quadrature scheme not defined')
            end
        end
        
        function w = get.weights(obj)
            switch obj.quadratureType
                case 'Legendre'
                    [~,w] = obj.GaussLegendre(obj.nSolutionPoints);
                case 'Laguerre'
                    [~,w] = obj.GaussLaguerre(obj.nSolutionPoints);
                case 'LGL'
                    [~,w] = obj.GaussLobatto(obj.nSolutionPoints);
                case 'Radau'
                    [~,w] = obj.GaussRadau(obj.nSolutionPoints);
                case 'Chebyshev'
                    [~,w] = obj.GaussChebyshev(obj.nSolutionPoints);
                otherwise
                    error('quadrature scheme not defined')
            end
        end
                
        function nN = get.nNodes(obj) % Element Nodes
            nN = obj.nElements*obj.nSolutionPoints;
        end
        
        function eFaces = get.elementFaces(obj) % Element Faces
            Faces = linspace(obj.range(1),obj.range(2),obj.nFaces);
            eFaces(1,:) = Faces(1:end-1);
            eFaces(2,:) = Faces(2:end);
        end
        
        function eSize = get.elementSize(obj) % Element Size
            eSize = (obj.range(2)-obj.range(1))/obj.nElements;
        end
        
        function eCenter =  get.elementCenter(obj) % Element Center
               eCenter = (obj.elementFaces(2,:)+obj.elementFaces(1,:))/2;
        end
        
        function J = get.Jacobian(obj) % Jacobian dx/dxi
            J = obj.elementSize/2;
        end
        
        function eNodes = get.elementNodes(obj)
            eNodes = reshape(1:obj.nNodes,obj.nSolutionPoints,obj.nElements);
        end
        
        function nodeCoords = get.nodeCoordinates(obj) % x-Grid
            [xSPs,xc] = meshgrid(obj.elementCenter,obj.Jacobian*obj.solutionPoints);
            nodeCoords = xSPs+xc; % x-Grid
        end
        
    end % Methods
    
end % Class

