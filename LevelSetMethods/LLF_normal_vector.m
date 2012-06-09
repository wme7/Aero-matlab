function [H, alpha_x, alpha_y] = LLF_normal_vector(dx, dy, Vn_ext, u_ext, v_ext, phi_x_minus, phi_x_plus, phi_y_minus, phi_y_plus)
%
% Estimate H (approximately) using Stencil Local 
% Lax-Friedrichs (SLLF) scheme.
% This scheme is used if phi is not approximately 
% a signed distance function. Adds some amount of
% dissipation, which will cause some smoothing.
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

phi_x = (phi_x_minus + phi_x_plus)/2;
phi_y = (phi_y_minus + phi_y_plus)/2;
abs_grad_phi = sqrt(phi_x.^2 + phi_y.^2);

phi_x_diff = (phi_x_plus - phi_x_minus)/2;
phi_y_diff = (phi_y_plus - phi_y_minus)/2;

%calculate alphas here

max_phi_x = max(phi_x_minus, phi_x_plus);
min_phi_x = min(phi_x_minus, phi_x_plus);
max_phi_y = max(phi_y_minus, phi_y_plus);
min_phi_y = min(phi_y_minus, phi_y_plus);

[r c] = size(Vn_ext);
dummy = max(min_phi_x(:)).*ones(r+6,c+6,7);

dummy(4:end-3, 1:end-6, 1) = min_phi_x;
dummy(4:end-3, 2:end-5, 2) = min_phi_x;
dummy(4:end-3, 3:end-4, 3) = min_phi_x;
dummy(4:end-3, 4:end-3, 4) = min_phi_x;
dummy(4:end-3, 5:end-2, 5) = min_phi_x;
dummy(4:end-3, 6:end-1, 6) = min_phi_x;
dummy(4:end-3, 7:end, 7) = min_phi_x;
min_stencil_phi_x = min(dummy, [], 3);
min_stencil_phi_x = min_stencil_phi_x(4:end-3,4:end-3);

dummy = max(min_phi_y(:)).*ones(r+6,c+6,7);

dummy(1:end-6, 4:end-3, 1) = min_phi_y;
dummy(2:end-5, 4:end-3, 2) = min_phi_y;
dummy(3:end-4, 4:end-3, 3) = min_phi_y;
dummy(4:end-3, 4:end-3, 4) = min_phi_y;
dummy(5:end-2, 4:end-3, 5) = min_phi_y;
dummy(6:end-1, 4:end-3, 6) = min_phi_y;
dummy(7:end, 4:end-3, 7) = min_phi_y;
min_stencil_phi_y = min(dummy, [], 3);
min_stencil_phi_y = min_stencil_phi_y(4:end-3,4:end-3);

clear dummy;


% approximately
app_mag_grad = sqrt(min_phi_x.^2 + min_stencil_phi_y.^2);
alpha_x = abs(u_ext) + abs(Vn_ext.*max_phi_x)./(app_mag_grad + dx*dx*(app_mag_grad<=1e-6));

%alpha_x = abs(u_ext);
app_mag_grad = sqrt(min_stencil_phi_x.^2 + min_phi_y.^2);
alpha_y = abs(v_ext) + abs(Vn_ext.*max_phi_y)./(app_mag_grad + dy*dy*(app_mag_grad<=1e-6));
%alpha_y = abs(v_ext);

% surf(app_mag_grad);
% pause;
%surf(alpha_y);
%pause;

H = u_ext.*phi_x + v_ext.*phi_y + Vn_ext.*abs_grad_phi - alpha_x.*phi_x_diff - alpha_y.*phi_y_diff;

%H = u_ext.*phi_x + v_ext.*phi_y - alpha_x.*phi_x_diff - alpha_y.*phi_y_diff;



