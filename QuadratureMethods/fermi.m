function xx=fermi(aj,eta)
% Matlab m file to evaluate the Half-order Fermi-Dirac integral
% For semiconductor/solid state applications
% Half-order implies that j (= aj)= -1/2 or 1/2 or 3/2 or 5/2

% Created on 5-Jan-2007 by Natarajan and Mohankumar
%========================================================

% Based on the Published work below 
% Title: The accurate numerical evaluation of half-order
%        Fermi-Dirac integrals
% Authors: N.Mohankumar & A.Natarajan
% Journal: Physica Status Solidi(b) vol.188, 1995, pp. 635-644
%============================================================

% xx is a dummy output
% eta and aj are input values to be given by the user.
% Printed output is the integral value.
% Input variable aj is same as j in the definition of F(j,eta) given below

% Accuracy: For solid state physics applications, eta values
% typically lie in the interval [-5,25]; For eta in [-5,25],
% the code below is guaranteed to yield 14 digit accuracy. 
%============================================================

% Integral definition is given below
% 
% Integral is defined by F(j,eta)
% Integrand is defined by (x^j)/[Gamma(1+j)* (exp(x-eta)+1)]
% Lower limit=0; Upper limit=infinity
% Integrated with respect to x 
% Please check whether your definition of the integral
% needs the Gamma(1+j) factor or not.
%===================================================================

% Trapezoidal Integration in y after the transformation
% from the original integration variable x to y
% where x= y^2
% Residue correction for the poles of the transformed integrand is
% added to the trapezoidal integration sum to expedite convergence
%=========================================================

%Program begins

format long e;
%==============================================================
% Evaluation of Trapezoidal sum begins
range=8.;
if eta > 0.
   range=sqrt(eta+64.);end;
h=0.5;
nmax=range/h;
sum=0.;
if aj== (-0.5)
   sum=1./(1.+exp(-eta));end;
for i=1:nmax
   u=i*h;
   ff=2.*(u^(2.*aj+1))/(1.+exp(u*u-eta));
   sum=sum+ff;end;

%Trapezoidal Summation ends
%==============================================================

% Pole correction for trapezoidal sum begins
pol=0.;
npole=0;
% Fix the starting value of  BK1 to start while loop
bk1=0;
while bk1 <= 14.*pi
   npole=npole+1;
   bk=(2*npole -1)*pi;
   rho=sqrt(eta*eta+bk*bk);
   t1=1;
   t2=0;
   if eta < 0;
      tk=- aj*(atan(-bk/eta)+pi);
   elseif eta ==0;
      tk=0.5*pi*aj;
   else;
      eta > 0;
      tk=aj*atan(bk/eta);
   end;
   rk=- (rho^aj);
   tk=tk+0.5*atan(t2/t1);
   if eta < 0
      rk= -rk;
   end;
   ak=(2.*pi/h)*sqrt(0.5*(rho+eta));
   bk1=(2.*pi/h)*sqrt(0.5*(rho-eta));
   if bk1 <= (14.*pi)
   gama=exp(bk1);
   t1=gama*sin(ak+tk)-sin(tk);
   t2=1.-2.*gama*cos(ak)+gama*gama;
   pol=pol+4.*pi*rk*t1/t2;
   end; %ends if loop above
  end; % Top while loop ends
  npole=npole-1;
  fdp=sum*h+pol;
% Program ends with the following output
%===============================================================================
  disp('Fermi-Dirac Integral Value');
  disp(fdp/gamma(1+aj));
  disp('Number of trapezoidal points & number of poles');
  disp([round(nmax),npole]);