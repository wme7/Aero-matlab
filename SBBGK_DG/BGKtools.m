classdef BGKtools
    %BGKTOOLS class
    %   This is the collection of fixed subroutines required to compute
    %   Kinetic Boltzmann equation using Discrete Ordinate Method.
    %   coded by Manuel Diaz, NTU, 2012.05.15
    
    properties
    end
    
    methods (Static)
        
        function [a,b,c] = apply_DG_DOM(a,b,c,nv)            % Repeat arrays 'nv' times
            a = repmat(a,[1,1,nv]);     
            c = repmat(c,[1,1,nv]);
            b = repmat(b,[1,1,nv]);
        end
        
        function p = bisection(f,a,b,tol)
            % provide the equation you want to solve with R.H.S = 0 form.
            % Write the L.H.S by using inline function
            % Give initial guesses.
            % Solves it by method of bisection.
            % A very simple code and very handy!
            if f(a)*f(b)>0
                disp('Wrong choice Bro')
            else
                p = (a + b)/2;
                err = abs(b-a);
                while err >= tol
                    if f(a)*f(p)<0
                        b = p;
                    else
                        a = p;
                    end
                    p = (a + b)/2;
                    err = abs(b-a);
                end
            end
        end
        
        function f = f_equilibrium_1d(z,u,c,t,theta)
            % Compute equilibrum in 1d cases for:
            %
            % MB:  Maxwell-Boltzmann, theta =  0
            % FD:  Fermi-Diract,      theta = +1
            % BE:  Bose-Einstein,     theta = -1
            %
            % Inputs:
            % u: macroscopic velocity
            % x: microscopic velocity
            % t: temperature
            % z: fugacity
            f = 1./(exp( (c-u).^2 ./ t)./z + theta);
        end
        
        function f = f_ES_equilibrium_1d(z,p,n,u,c,t,theta)
            % Compute equilibrum in 1d cases for:
            %
            % MB:  Maxwell-Boltzmann, theta =  0
            % FD:  Fermi-Diract,      theta = +1
            % BE:  Bose-Einstein,     theta = -1
            %
            % inputs:
            % u: macroscopic velocity
            % x: microscopic velocity
            % t: temperature
            % r: fugacity
            % p: pressure
            % n: density
            f = 1./(exp( (c-u).^2 ./ t)./z + theta);
        end
        
        function FD = FermiDiract(z,nu)
            % Fermi-Diract function:
            % Implementation of the Fermi-Diract Statistical function
            % in latex words:
            %
            % $$F_\nu(z)=\frac{1}{\Gamma(\nu)} \int_{0}^{\infty}\frac{x^{\nu-1}}{z^{-1}e^x+1}
            %   \approx \sum_{l=1}^{\infty}(-1)^{l-1}\frac{z^l}{l^\nu}$$
            %
            % As we can notice this function is of the order $\nu$
            l   = 1:50; % up to l = 50 to ensure a fair accuaracy
            fd  = (-1).^(l-1) .* (z.^l) ./ (l.^nu);
            FD  = sum(fd);
        end
        
        function BE = BoseEinstein(z,nu)
            % Bose-Einstein function:
            % Implementation of the Bose-Einstein Statistical function
            % in latex words:
            %
            % $$F_\nu(z)=\frac{1}{\Gamma(\nu)}\int_{0}^{\infty}\frac{x^{\nu-1}}{z^{-1}e^x-1}
            %   \approx \sum_{l=1}^{\infty}\frac{z^l}{l^\nu}$$
            %
            % As we can notice this function is of the order $\nu$
            l   = 1:50; % up to l = 50 to ensure a fair accuaracy
            be  = (z.^l) ./ (l.^nu);
            BE  = sum(be);
        end
        
        function [x,w,k] = cotes_xw(a,b,N,nodes)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            % Written by: Manuel Diaz
            % Ref.: http://mathworld.wolfram.com/Newton-CotesFormulas.html
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            N = (nodes-1)*ceil(N/(nodes-1));  N1=N+1;
            x = linspace(a,b,N1)';
            h = x(2)-x(1);
            switch nodes
                case{2} % Trapezoidal Rule
                    w = kron(ones(1,N/1),(2))'; w(1)=1; w(N1)=1;
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
        
        function [x, w] = GaussHermite(n)
            % This function determines the abscisas (x) and weights (w) for the
            % Gauss-Hermite quadrature of order n>1, on the interval [-INF, +INF].
            % This function is valid for any degree n>=2, as the companion matrix
            % (of the n'th degree Hermite polynomial) is constructed as a
            % symmetrical matrix, guaranteeing that all the eigenvalues (roots)
            % will be real.
            %
            % � Geert Van Damme
            % geert@vandamme-iliano.be
            % February 21, 2010
            %
            % Building the companion matrix CM
            % CM is such that det(xI-CM)=L_n(x), with L_n the Hermite polynomial
            % under consideration. Moreover, CM will be constructed in such a way
            % that it is symmetrical.
            i   = 1:n-1;
            a   = sqrt(i/2);
            CM  = diag(a,1) + diag(a,-1);
            % Determining the abscissas (x) and weights (w)
            % - since det(xI-CM)=L_n(x), the abscissas are the roots of the
            %   characteristic polynomial, i.d. the eigenvalues of CM;
            % - the weights can be derived from the corresponding eigenvectors.
            %[V L]   = eig(CM);
            [V L]   = trideigs(diag(CM,0),diag(CM,1),1e-12,100);
            [x ind] = sort(diag(L));
            V       = V(:,ind)';
            w       = sqrt(pi) * V(:,1).^2;
        end
        
    end
    
    methods         
    end
        
end

