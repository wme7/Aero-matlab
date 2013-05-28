%%%%%%%%%%%%%%%%%%%%%%%
%    R E S U L T S
%%%%%%%%%%%%%%%%%%%%%%%

%% Using T3 elements,

% Max y-displacement
% Node           UX           UY
%   8   1.916937e-06    -1.065042e-05
T3_disp = -1.065042e-05;

% Stress in element 5 [MPA]
T3_Sigma_xx_e5 = -650.431082; 
T3_Sigma_yy_e5 = -650.431082; 
T3_Sigma_xy_e5 = -650.431082; 
T3_Vonmises_e5 =  1300.862163; 

% Stress in element 6 [MPA]
T3_Sigma_xx_e6 =  650.431082; 
T3_Sigma_yy_e6 =  645.504435; 
T3_Sigma_xy_e6 = -1349.568918; 
T3_Vonmises_e6 =  2425.672941;

% Average stress in Node 8 [MPA]

T3_Sigma_xx = (T3_Sigma_xx_e6+T3_Sigma_xx_e5)/2
T3_Sigma_yy = (T3_Sigma_yy_e6+T3_Sigma_yy_e5)/2 
T3_Sigma_xy = (T3_Sigma_xy_e6+T3_Sigma_xy_e5)/2
T3_Vonmises = (T3_Vonmises_e6+T3_Vonmises_e5)/2

% T3_Sigma_xx = 0
% T3_Sigma_yy = -2.4633
% T3_Sigma_xy = -1000
% T3_Vonmises = 1.8633e+03

%% Using Q4 elements,

% Max y-displacement
% Node           UX           UY
%    8   6.1172e-06    -2.6430e-05
Q4_disp = -2.6430e-05;

% Stress in element 3 

% Stress in node 3 [MPA]
Q4_Sigma_xx = 1947.025060;
Q4_Sigma_yy = -837.861465;
Q4_Sigma_xy = -543.285579;
Q4_Vonmises = 2647.590101;

%% Using Q8 elements,

% Max y-displacement
% Node           UX           UY
%   18   9.1557e-06    -3.8250e-05
Q8_disp = -3.8250e-05;

% Stress in element 3 

% Stress in node 3 [MPA]
Q8_Sigma_xx = 547.919444;
Q8_Sigma_yy = -5070.115307;
Q8_Sigma_xy = -1926.948499;
Q8_Vonmises = 6318.519705;

%% Using Q12 elements,

% Max y-displacement
% Node           UX           UY
%   28   9.3240e-06    -3.8900e-05
Q12_disp = -3.8900e-05;

% Stress in element 3

% Stress in node 3 [MPA]
Q12_Sigma_xx = -209.968713;
Q12_Sigma_yy = -7862.634791;
Q12_Sigma_xy = -2337.689116;
Q12_Vonmises = 8752.632553;

%%%%%%%%%%%%%%%%%%%%%%%
%    C O M P A R E
%%%%%%%%%%%%%%%%%%%%%%%
x = 1:4; 
str = ['T3','Q4','Q8','Q12'];
disp_y   = [T3_disp,Q4_disp,Q8_disp,Q12_disp];
Sigma_xx = [T3_Sigma_xx,Q4_Sigma_xx,Q8_Sigma_xx,Q12_Sigma_xx];
Sigma_yy = [T3_Sigma_yy,Q4_Sigma_yy,Q8_Sigma_yy,Q12_Sigma_yy];
Sigma_xy = [T3_Sigma_xy,Q4_Sigma_xy,Q8_Sigma_xy,Q12_Sigma_xy];
VonMises = [T3_Vonmises,Q4_Vonmises,Q8_Vonmises,Q12_Vonmises];
figure(1) 
subplot(2,3,1); scatter(x,disp_y,'ob'); title('Max y-displacement'); ylabel('[mm]'); xlabel('T3 - Q4 - Q8 - Q12');
subplot(2,3,2); scatter(x,Sigma_xx,'ok'); title('\sigma_{xx}'); ylabel('[MPa]'); xlabel('T3 - Q4 - Q8 - Q12');
subplot(2,3,3); scatter(x,Sigma_yy,'dr'); title('\sigma_{yy}'); ylabel('[MPa]'); xlabel('T3 - Q4 - Q8 - Q12');
subplot(2,3,4); scatter(x,Sigma_xy,'sb'); title('\sigma_{xy}'); ylabel('[MPa]'); xlabel('T3 - Q4 - Q8 - Q12');
subplot(2,3,5);scatter(x,VonMises,'xk'); title('VonMises'); ylabel('[MPa]'); xlabel('T3 - Q4 - Q8 - Q12');
