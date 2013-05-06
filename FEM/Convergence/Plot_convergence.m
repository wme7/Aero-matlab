%% Plot Errors

%Run first
quad_disp;
linear_disp;

%Plot
figure(2)
loglog(h,L2_linear,'-sr'); hold on; 
loglog(h,L2_quad,'-sb'); 
loglog(h,en_linear,'--*r'); 
loglog(h,en_quad,'--*b'); hold off;
legend('L_2 Linear','L_2 Quadratic','en Linear','en Quadratic','interpreter','latex',4);
xlabel('log_{10}h','interpreter','tex');
box on;