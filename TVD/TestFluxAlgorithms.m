%% Testing Fluxes Algorithms
% By Manuel Diaz 2012.12.18

clear all; clc;

%% Initial Parameters
 x = 0:0.1:1;
 u = 1/2+sin(2*pi*x); %IC
 nx =length(u);

%% Define our Flux function
     f = @(w) 1.*w; % w.^2/2
% and the Derivate of the flux function
    df = @(w) w./w; % =1 % w
    
%% base data
flx = f(u);     % flux value at every point of the domain
dflx = df(u);   % flux slope at every point of the domain
 
%% Algortithms to compare h
%**************************************************************************
% % Roe with entropy fix
%         
%         % Choose f(a) if f'(u) >= 0 for u that belongs [min(a,b),max(a,b)]
%         % Choose f(b) if f'(u) <= 0 for u`that belongs [min(a,b),max(a,b)]
%         % Otherwise choose LLF flux
%         for i = 1:nx-1 % all middle poinst
%             u_ave(i) = (u(i) + u(i+1))/2;
%             if df(u_ave(i)) > 0 % if positive
%                 h(i) = f(u(i));
%             elseif df(u_ave(i)) <= 0 % if negative
%                 h(i) = f(u(i+1));
%             else
%                 h(i) = 0;
%             end
%         end
% 
% %  Roe Flux - - matlab optimized - Test Passed!
tic
        u_ave = (u(1:nx-1) + u(2:nx))/2;    % u @ cells boundaries
        bool_p = df(u_ave) > 0; % if positive
        bool_n = df(u_ave) <= 0;  % if negative
        h_roe = bool_p.*flx(1:nx-1) + bool_n.*flx(2:nx)
toc
%**************************************************************************        
% % (Global) Lax Friedrichs
%         
%         % Compute  max(|f'(u)|) 
%         % This max value is computed over the entire region of u, that is
%         % [inf u(x), sup u(x)] where u(x) is the initial function.
%         alpha = max(abs(dflx));
%         % Compute Lax-Friedrichs Flux
%         for i = 1:nx-1 % for all the middle points
%             h(i) = 0.5*( flx(i) + flx(i+1) - alpha *( u(i+1) - u(i) ));
%         end
% 
% % LF - - matlab optimized - Test Passed!
tic
        alpha = max(abs(dflx));
        h_LF = 0.5*(flx(1:nx-1) + flx(2:nx) - alpha*(u(2:nx) - u(1:nx-1)))
toc
%**************************************************************************
% % Local Lax-Friedrichs
%         
%      % Similarly for Lax-Friedrichs:
%      for i = 1:nx-1 % for all the middle points
%         beta = max(abs(dflx(i)),abs(dflx(i+1)));
%         h(i) = 0.5*( flx(i) + flx(i+1) - beta *( u(i+1) - u(i) ));
%      end
% LLF - - matlab optimized - Test Passed!
tic
        fluxmat = [dflx(1:nx-1);dflx(2:nx)];
        beta = max(abs(fluxmat));
        h_LLF = 0.5*(flx(1:nx-1) + flx(2:nx) - beta.*(u(2:nx) - u(1:nx-1)))
toc
%**************************************************************************
% % Upwind
%     % For dflux is constant along the domain!
%     % if dflux > 0, flux to the left, then h(a,b) = h(a)
%     % if dflux > 0, flux to the left, then h(a,b) = h(b)
%     % We evaluate professor yang's strategy:
%     % % a = (a - |a|)/2
%     a = (dflx(1) - abs(dflx(1)))/2; 
%     for i = 1:nx-1 % for all middle points
%         if a == 0 % a > 0 
%             h(i) = flx(i); % Flux to the left
%         else % a < 0 
%             h(i) = flx(i+1); % Flux to the right
%         end
%     end
    
% Upwind - - matlab optimzed
tic
    % For dflux is constant along the domain!
    a_p = max(dflx - abs(dflx))/2 == [0]; %#ok<NBRAK> % boolen operator for a>0
    a_n = max(dflx + abs(dflx))/2 == [0]; %#ok<NBRAK> % boolen operator for a<0
    h_UPWIND  = a_p*flx(1:nx-1) + a_n*flx(2:nx)
toc
%**************************************************************************
% Conclusion: 
% 1.The fastest is LF
% 2.The second fastes is Upwind
% 3.Best results are with LLF
% 4.Roe is very easy to undestand
% 5.Upwind is not conservative therefore only works fine and fast when df=a
% END