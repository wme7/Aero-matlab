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
        
        function x = Chebyshev(n)
            x = sin(0.5*pi*linspace(-1,1,n+1)');
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
                otherwise
                    error('quadrature scheme not defined')
            end
        end
        
        function w = get.weights(obj)
            switch obj.quadratureType
                case 'Legendre'
                    [x,w] = obj.GaussLegendre(obj.nSolutionPoints);
                case 'Laguerre'
                    [x,w] = obj.GaussLaguerre(obj.nSolutionPoints);
                case 'LGL'
                    [x,w] = obj.GaussLobatto(obj.nSolutionPoints);
                case 'Radau'
                    [x,w] = obj.GaussRadau(obj.nSolutionPoints);
                otherwise
                    error('quadrature scheme not defined')
            end
            clear('x');
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

