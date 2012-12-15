%-------------------------------------------------------------------
% MATLAB routine for solving Burgers equation 
%   v_t + (v^2/2)_x  (conservation form)
%
% created: April 2010
% author:  Susanne Hoefner, Uppsala University 
%-------------------------------------------------------------------

clear   % clear all variables
clf     % clear plot area

%-------------------------------------------------------------------
% defining spatial grid
%-------------------------------------------------------------------

npt   = 200;     % number of spatial grid points
x_min = 0.0;     % left boundary
x_max = 2.0;     % right boundary

dx    = (x_max-x_min)/(npt-1);  % grid spacing (equidistant)
for i=1:npt
  x(i) = x_min + (i-1)*dx;      % setting spatial grid points 
end

x_diff  = (x_max-x_min);        % length of spatial range
x_mid   = (x_min+x_max)/2.0;    % mid-point of spatial range

%-------------------------------------------------------------------
% setting initial conditions
%-------------------------------------------------------------------

% parameters for initial data

%i_dat = 1;      % Gaussian
%i_dat = 2;      % slope
i_dat = 3;      % triangle

v_min = +0.20;  % maxium velocity
v_max = +0.70;  % minimum velocity

f_tip   = 0.30; % location of gaussian and triangle
f_gauss = 0.40; % width parameter for gaussian
f_width = 0.20; % width parameter for triangle and slope

% creating initial data

for i=1:npt

  % gaussian

  if i_dat == 1
    x_tip   = x_min + f_tip*x_diff;  % center of triangle
    v(i,1) = (v_max-v_min) * exp( -( (x(i)-x_tip) / f_gauss*x_diff )^2 ) + v_min;
  end

  % slope

  if i_dat == 2
    x_left  = x_mid-f_width*x_diff;
    x_right = x_mid+f_width*x_diff;
    v(i,1) = v_max;
    if x(i) > (x_left)
      if x(i) < (x_right)
        v(i,1) = v_max - (v_max-v_min)*((x(i)-x_left)/(x_right-x_left));
      else
        v(i,1) = v_min;
      end
    end
  end

  % triangle

  if i_dat == 3
    x_tip   = x_min + f_tip*x_diff;  % center of triangle
    x_left  = x_tip-f_width*x_diff;  % left border of triangle
    x_right = x_tip+f_width*x_diff;  % right border of triangle
    v(i,1) = v_min;                  % velocity outside triangle region
    if x(i) > x_left
      if x(i) < x_right
        if x(i) < x_tip              % left side of triangle
          v(i,1) = v_min + (v_max-v_min)*((x(i)-x_left)/(x_tip-x_left));
        else                         % right side of triangle
          v(i,1) = v_max - (v_max-v_min)*((x(i)-x_tip)/(x_right-x_tip));
        end
      end
    end
  end

end

%-------------------------------------------------------------------
% computing solution at chosen time levels
%    (using a conservative explicit upwind differencing scheme)
%-------------------------------------------------------------------

t_start = 0.0;  % starting time 
t_stop  = 2.4;  % stopping time 

cfl     = 0.8;  % Courant number (must be < 1 for stability)

t(1)    = t_start; % time level of initial data
age     = t(1);    % control variable for time integration
n       = 0;       % initializing time level counter

while (age <= t_stop) 

  n      = n+1;                      % index of time level
  dt     = dx*cfl/max(abs(v(:,n)));  % time step
  t(n+1) = t(n) + dt;                % time 
  age    = t(n+1);

  v(1,n+1) = v(1,1);       % left boundary condition (fixed value)
  v(npt,n+1) = v(npt,1);   % right boundary condition (fixed value)

  % compute upwind fluxes at cell boundaries

  for i=1:npt-1            
    vp = 0.5 * ( v(i,n) + v(i+1,n) ); % average value
    if (vp > 0.0) % average value positive
      if (v(i,n) >= 0.0) 
        flp(i) = 0.5 * v(i,n)^2;
      end
    elseif (vp <= 0.0) % average value negative
      if (v(i+1,n) <= 0.0) 
        flp(i) = 0.5 * v(i+1,n)^2;
      end
    else
      flp(i) = 0.0;
    end
  end  

  % compute solution v for new time level t(n+1) at locations x(i)
  
  for i=2:npt-1
    v(i,n+1) = v(i,n) - dt/dx * ( flp(i) - flp(i-1) );
  end
  
end

ntst = n+1;  % total number of time levels (including initial data)

%-------------------------------------------------------------------
% plotting initial data and solution
%-------------------------------------------------------------------

vdata_max = max(v(:,1));
vdata_min = min(v(:,1));
dv = vdata_max - vdata_min;

fb    = 0.4;              % parameter for extending y-range
y_min = vdata_min-dv*fb;
y_max = vdata_max+dv*fb;

hold on 

title ('Burgers equation');  
ylabel('velocity'); 
xlabel('x'); 

axis([x_min x_max y_min y_max]);

plot(x,v(:,1),'r-')       % plot initial data
plot(x,v(:,ntst),'b-')    % plot solution at last time level

legend('inital data','numerical solution')

npl  = 5; % number of snapshots to be plotted

if npl > 1
  dn = (ntst-mod(ntst,npl))/npl;
  for i=dn:dn:(ntst-npl)
    plot(x,v(:,i),'b--')   % plot intermediate snapshots
  end
end

hold off

%-------------------------------------------------------------------