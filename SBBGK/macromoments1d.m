function [n, j_x, E] = macromoments1d(k,w,f,v)
%% Compute Macroscopic Moments
% Using Quadrature rules to integrate for:
     n   = k*sum(w .* f);    % Number density
     j_x = k*sum(v .* w .* f);   % Macrospic moment in x
     E   = k*sum(1/2*( v.^2 ).* w .* f); % Energy Density

% Using Quadrature rules to integrate for:
%[nv,nx] = size(f);
% for J = 1:nx
%     SR = 0;
%     SU = 0;
%     SE = 0;
%     %SAV= 0;
%     for K = 1:nv
%         SR = SR + k*w(K,J) * f(K,J);
%         SU = SU + k*w(K,J) * f(K,J) * v(K,J);
%         SE = SE + k*w(K,J) * f(K,J) * (0.5 * v(K,J) * v(K,J));
%         %SAV = SAV + k*w(K,J) * f(K,J) * abs(v(K,J));
%     end
%     n(J)    = SR;
%     j_x(J)  = SU;
%     %U(J)    = SU/SR;
%     E(J)   = SE;
%     %AV(J)   = SAV;
% end