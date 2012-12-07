function u = u_zero(x,case_no)
% Create vector u with an Initial Condition. 4 cases are available.
%**************************************************************************
%
% 4 cases available: 
% {1} Gaussian Jump
% {2} Square Jump
% {3} Sinusoidal
% {4} Riemann IC
%
% Coded by Manuel Diaz 2012.12.06
%**************************************************************************
%% Allocate vector u
[n m] = size(x); u = zeros(n,m);

%% Compute case_no
switch case_no
    case{1} % Gaussian Jump IC
        for j = 1:n
            for i = 1:m
                if (x(j,i)>0.4) && (x(j,i)<0.6)
                    u(j,i)=(x(j,i)-0.4).^10*(x(j,i)-0.6).^10*10^20;
                else
                    u(j,i)=0;
                end
            end
        end
    case{2} % Square Jump IC
        for j = 1:n
            for i = 1:m
                if (x(j,i)>0.4) && (x(j,i)<0.6)
                    u(j,i)=1;
                else
                    u(j,i)=0;
                end
            end
        end
    case{3} % Sinusoidal IC
        
        u = 1/2+sin(pi*x);
        
    case{4} % Riemann IC
        for j = 1:n
            for i = 1:m
                if (x(j,i)<=1)
                    u(j,i)=1;
                else
                    u(j,i)=0;
                end
            end
        end
    otherwise 
        error('Case not available')
end