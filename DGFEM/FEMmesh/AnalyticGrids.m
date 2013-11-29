classdef AnalyticGrids
    %ANALYTICGRIDS Compute points in other curvilinear systems
    % Coordinates systems available:
    %   PolarCoordinates
    %   ParabolicCylinderCoordinates
    %   EllipticCylinderCoordinates
    %   Horseshoe
    %   ModifiedHorseshoe
    %   BipolarCoordinates
    
    properties
        coordinate
        xi
        eta
    end
    
    properties (Dependent = true, SetAccess = private)
        x
        y
    end
    
    methods (Static)
        function [x y] = PolarCoordinates(Xi,Eta,theta0,theta1,r0,r1)
            % theta0, theta1 are two angles.
            % 0<=theta0<theta1<=2.*pi
            %
            % r0, r1 are two radii
            % 0<=r0<r1
            %
            % Check for the input arguments and set default values
            if nargin == 2
                theta0 = 0.;
                theta1 = -pi/4;
                
                r0 = 1.;
                r1 = 2.;
            end
            theta = (theta0-theta1).*Xi;
            r = r0+(r1-r0).*Eta;
            x = r.*cos(theta);
            y = r.*sin(theta);
        end
        
        function [x y] = BipolarCoordinates(Xi,Eta,a)
            % a is radii that defies height of the domain
            % a > 0
            % Check the input arguments and set default values
            %
            if nargin == 2
                a = 1.;
            end
            r = Xi;
            s = pi.*(Eta-1/2);
            x = (a.*sinh(r))./(cosh(r)+cos(s));
            y = (a.*sin(s))./(cosh(r)+cos(s));
        end
        
        function [x y] = EllipticCylinderCoordinates(Xi,Eta,a)
            % a is radius of ellips
            % a > 0
            %
            if nargin == 2
                a = 1.;
            end
            r = 1+Xi;
            s = pi.*Eta;
            x = a.*cosh(r).*cos(s);
            y = a.*sinh(r).*sin(s);
        end
        
        function [x y] = Horseshoe(Xi,Eta,row,b0,b1)
            % row is aspect ratio of the major to minor axes of ellipse
            % row > 0
            % b0, b1 are two radii
            % 0<b0<b1
            % check for the input arguments and set default values
            if nargin == 2
                row = 3.;
            end
            if nargin == 2
                b0 = 4.;
                b1 = 7.;
            end
            r = b0+(b1-b0).*Eta;
            theta = pi/2.*(1-2.*Xi);
            x = row.*r.*cos(theta);
            y = r.*sin(theta);
        end
        
        function [x y] = ModifiedHorseshoe(Xi,Eta,R)
            % R is aspect ratio of the major to minor axes of ellipse
            % R >= 1
            % Check the input arguments and set default values
            if nargin == 2
                R = 3.;
            end
            x = -(1+Eta).*cos(pi.*Xi);
            y = (1+(2.*R-1).*Eta).*sin(pi.*Xi);
        end
        
        function [x y] = ParabolicCylinderCoordinates(Xi,Eta)
            r = 1+Xi;
            s = 1+Eta;
            x = 1/2.*(r.^2-s.^2);
            y = r.*s;
        end
        
    end % Static Methods
    
    methods
        function obj = AnalyticGrids(xi,eta,coordinate) % constructor
            obj.xi = xi;
            obj.eta = eta;
            obj.coordinate = coordinate;
        end
        
        function x = get.x(obj)
            switch obj.coordinate
                case 'PolarCoordinates'
                    [x,~] = obj.PolarCoordinates(obj.xi,obj.eta);
                case 'ParabolicCylinderCoordinates'
                    [x,~] = obj.ParabolicCylinderCoordinates(obj.xi,obj.eta);
                case 'EllipticCylinderCoordinates'
                    [x,~] = obj.EllipticCylinderCoordinates(obj.xi,obj.eta);
                case 'Horseshoe'
                    [x,~] = obj.Horseshoe(obj.xi,obj.eta);
                case 'ModifiedHorseshoe'
                    [x,~] = obj.ModifiedHorseshoe(obj.xi,obj.eta);
                case 'BipolarCoordinates'
                    [x,~] = obj.BipolarCoordinates(obj.xi,obj.eta);
                otherwise
                    error('coordinate system not defined')
            end
        end
        
        function y = get.y(obj)
            switch obj.coordinate
                case 'PolarCoordinates'
                    [~,y] = obj.PolarCoordinates(obj.xi,obj.eta);
                case 'ParabolicCylinderCoordinates'
                    [~,y] = obj.ParabolicCylinderCoordinates(obj.xi,obj.eta);
                case 'EllipticCylinderCoordinates'
                    [~,y] = obj.EllipticCylinderCoordinates(obj.xi,obj.eta);
                case 'Horseshoe'
                    [~,y] = obj.Horseshoe(obj.xi,obj.eta);
                case 'ModifiedHorseshoe'
                    [~,y] = obj.ModifiedHorseshoe(obj.xi,obj.eta);
                case 'BipolarCoordinates'
                    [~,y] = obj.BipolarCoordinates(obj.xi,obj.eta);
                otherwise
                    error('coordinate system not defined')
            end
        end
    end % Methods

end