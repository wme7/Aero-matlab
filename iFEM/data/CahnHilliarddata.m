function pde = CahnHilliarddata(h0)
%\frac{\partial c}{\partial t} - \nabla(\nabla(\frac{df}{dc}-\lambda\nabla^{2}c))= 0 \\Omega,
%(\nabla\left(\frac{d f}{d c} - \lambda \nabla^{2}c\right)\right)= 0 \partial\Omega
%\lambda \nabla c \cdot n &=& 0 \partial\Omega.
%% there are there cases come from 
%A nonconforming finite element method for the Cahn?Hilliard equation; and
%http://fenicsproject.org/documentation/dolfin/dev/python/demo/pde/cahn-hilliard/python/documentation.html
%global efsilonsquare 

pde = struct('initial_u',@initial_u,'initial_w',@initial_w,'dfdu',@dfdu,'df2du2',@df2du2,'Du',@Du);
 

function z = initial_u(p)   % initial value
    

x = p(:,1);
y = p(:,2);
n = size(x,1);
%**************CASE ONE***************************************
% z = zeros(n,1); 
% x0 = 0.5;
% y0 = 0.5;
% w=pi/(4*h0); 
% xys = find((x-x0)<=(8*h0)&(y-y0)<=(8*h0)&(x-x0)>=0&(y-y0)>=0);
% z(xys) = (1.0e-3)*(sin(w*(x(xys)-x0))).^3.*(sin(w*(y(xys)-y0))).^3;
%**************CASE two***************************************
% n = size(x,1);
% z = zeros(n,1); 
% x0 = 0.0;
% y0 = 0.0;
% w=pi/(4*h0); 
% xys = find((x-x0)<=(8*h0)&(y-y0)<=(8*h0)&(x-x0)>=0&(y-y0)>=0);
% z(xys) = (1.0e-3)*(sin(w*(x(xys)-x0))).^3.*(sin(w*(y(xys)-y0))).^3;

%**************CASE Three***************************************

z = 0.63+0.02*(0.5-rand(n,1)); 
end

function z = initial_w(p)   % initial value
    

x = p(:,1);
y = p(:,2);
n = size(x,1);

%**************CASE ONE***************************************

% z = zeros(n,1);
% w = pi/(4*h0);
% x0 = 0.5;
% y0 = 0.5;
% xys = find((x-x0)<=8*h0&(y-y0)<=8*h0&(x-x0)>=0&(y-y0)>=0);
% z(xys) = (sin(2*w*(x(xys)-x0)).*cos(w*(x(xys)-x0))-(sin(w*(x(xys)-x0)).^3)).*(sin(w*(y(xys)-y0))).^3;
% z(xys) = z(xys)+(sin(2*w*(y(xys)-y0)).*cos(w*(y(xys)-y0))-(sin(w*(y(xys)-y0)).^3)).*(sin(w*(x(xys)-x0))).^3;
% z(xys) = 3*w^2*z(xys);
% z = -(1.0e-3)*efsilonsquare*z + dfdu(initial_u(p));


%**************CASE TWO***************************************
% z = zeros(n,1);
% w = pi/(4*h0);
% x0 = 0.0;
% y0 = 0.0;
% xys = find((x-x0)<=8*h0&(y-y0)<=8*h0&(x-x0)>=0&(y-y0)>=0);
% z(xys) = (sin(2*w*(x(xys)-x0)).*cos(w*(x(xys)-x0))-(sin(w*(x(xys)-x0)).^3)).*(sin(w*(y(xys)-y0))).^3;
% z(xys) = z(xys)+(sin(2*w*(y(xys)-y0)).*cos(w*(y(xys)-y0))-(sin(w*(y(xys)-y0)).^3)).*(sin(w*(x(xys)-x0))).^3;
% z(xys) = 3*w^2*z(xys);
% z = -(1.0e-3)*efsilonsquare*z + dfdu(initial_u(p));
%**************CASE Three***************************************

z = zeros(n,1);

end


 function z = dfdu(u) 
 %z = (u.^3-u);
 %**************CASE TWO***************************************
 %z = (u.^3-u);
 %**************CASE Three***************************************
 z = 100*(2*u-6*u.^2+4*u.^3);

 end


 function z = df2du2(u) 
     
%z = (3*u.^2-1);
%**************CASE TWO***************************************
% z = (3*u.^2-1);
%**************CASE Three***************************************

z = 100*(2-12*u+12*u.^2);
 end


end