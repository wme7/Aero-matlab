classdef TroubleCellDectector
    %TROUBLECELLDECTECTOR 
    %   Like the title says in 1D
    
    properties
        nCellsIn
        detector
        averages
        a
        M
        h
        uL
        uR
    end
    
    properties (Dependent = true, SetAccess = private)
        troubledCells
    end
    
    methods (Static)
        function mfunc = MODminmod(a,M,h)
            % MODminmod function
            % In the following implementation, it is assumed 'a' is a vector with 3
            % elements, 1<=j<=l where l = 3.
            % by Manuel Diaz, NTU, 2013.11.08
            mfunc = a(1,:);
            ids = find(abs(mfunc) > M*h.^2);
            
            if(size(ids,2)>0)
                mfunc(ids) = TroubleCellDectector.minmod(a(:,ids));
            end
%             if abs(a(1,:)) <= M*h^2
%                 m_tilde = a(1,:);
%             else
%                 m_tilde = TroubleCellDectector.minmod(a);
%             end
        end
        
        function mfunc = minmod(a)
            % minmod function
            % In the following implementation, it is assumed 'a' is a vector with 3
            % elements, 1<=j<=l where l = 3.
            % by Manuel Diaz, NTU, 2013.11.08
            m = size(a,1); mfunc = zeros(1,size(a,2));
            s = sum(sign(a),1)/m;
            
            ids = find(abs(s)==1);
            if(~isempty(ids))
                mfunc(ids) = s(ids).*min(abs(a(:,ids)),[],1);
            end
%             signtest = sign(a(1,:)) == sign(a(2,:)) & sign(a(1,:)) == sign(a(3,:));
%             s = sign(a(1,:));
%             m = s.*min(abs(a)).*signtest;
        end
    end % Methods
    
    methods
        function obj = TroubleCellDectector(u,strategy,a,uL,uR,M,h)
            obj.averages = u;
            obj.nCellsIn = length(u);
            obj.detector = strategy;
            obj.a = a;
            obj.M = M;
            obj.h = h;
            obj.uL = uL;
            obj.uR = uR;
        end
        
        function tCells = get.troubledCells(obj)
            switch obj.detector
                case 'TVDtheta'
                    tCells = obj.theta1d();
                case 'MODminmod'
                    tCells = obj.MODminmod1d();
                otherwise
                    error('detection strategy not available')
            end
        end
              
        function tCells = MODminmod1d(obj)
            u_bar = obj.averages;
            
            u_1tilde = obj.uR - u_bar;
            u_2tilde = u_bar - obj.uL;
            
            % Taking into account periodic BCs!
            %Dpu_bar = [u_bar(2:end),u_bar(1)] - u_bar;
            %Dnu_bar = [u_bar(end),u_bar(1:end-1)] - u_bar;
            
            % Taking into account Neumann BCs!
            Dpu_bar = [u_bar(2:end),u_bar(end)] - u_bar;
            Dnu_bar = [u_bar(1),u_bar(1:end-1)] - u_bar;
            
            % modif to take care of unsigned zeros
            A = [u_1tilde;Dpu_bar;Dnu_bar]; A(find([u_1tilde;Dpu_bar;Dnu_bar]==0)) = 1e-16;
            B = [u_2tilde;Dpu_bar;Dnu_bar]; B(find([u_2tilde;Dpu_bar;Dnu_bar]==0)) = 1e-16;
            
            uMOD_1tilde = obj.MODminmod(A,obj.M,obj.h);
            uMOD_2tilde = obj.MODminmod(B,obj.M,obj.h);
            
            % Mark troubled cells:
            tCells = find(uMOD_1tilde ==0 & uMOD_2tilde == 0);
        end
        
        function tCells = theta1d(obj)
            % Theta: smooth measuring factor (r_x)
            % INPUTS
            % u: Vector domain, u(j)
            % a: x-velocity in our domain, a(j)
            % Compute r_x
            % Initialize variables
            u = obj.averages;   a = obj.a;
            nx  = length(u);
            r_x = zeros(1,nx);
            % Test a
            % Check whether a is scalar or a vector of velocities with
            % prescribed velocities in the entire domain
            [m,n,r] = size(a);
            if m == 1 && n == 1 && r == 1
                a = a*ones(1,nx); % map the x-velocity
            else
                %do nothing
            end
            % Main loop
            x = 2:nx-1;
            for j = x
                % smooth measurement factor 'r_x'
                if u(j) == u(j+1)
                    r_x(j) = 1;
                elseif a(j) > 0
                    r_x(j) = (u(j) - u(j-1)) / (u(j+1) - u(j));
                elseif a(j) < 0
                    u(nx+1) = u(nx); % we will need an extra column value
                    r_x(j) = (u(j+2) - u(j+1)) / (u(j+1) - u(j));
                end
            end
            r_x(1) = 1; r_x(nx) = 1;
            % Troubled Cells
            tCells = find(r_x==0);
        end
                     
    end % Methods
    
end % Class

