function [u] = iBurgers1D(u, FinalTime)

% function [rho, rhou, Ener] = Euler1D(rho, rhou, Ener, FinalTime)
% Purpose  : Integrate 1D Euler equations until FinalTime starting with
%            initial conditions [rho, rhou, Ener]

Globals1D;

% Parameters
CFL = 1; time = 0;

% Prepare for adaptive time stepping
mindx = min(x(2,:)-x(1,:));

% Limit initial solution
u =SlopeLimitN(u);

% outer time step loop 
while(time<FinalTime)
  
  dt = CFL*min(min(mindx./(abs(u))));
  
  if(time+dt>FinalTime)
    dt = FinalTime-time;
  end

  % 3rd order SSP Runge-Kutta
  
  % SSP RK Stage 1.
  [rhsu]  = iBurgersRHS1D(u);
  u1  = u  + dt*rhsu;
  
  % Limit fields
  u1  = SlopeLimitN(u1);

  % SSP RK Stage 2.
  [rhsu]  = iBurgersRHS1D(u1);
  u2   = (3*u  + u1  + dt*rhsu )/4;
  
  % Limit fields
  u2  = SlopeLimitN(u2);

  % SSP RK Stage 3.
  [rhsu]  = iBurgersRHS1D(u2);
  u  = (u  + 2*u2  + 2*dt*rhsu )/3;
  
  % Limit solution
  u =SlopeLimitN(u);
  
  % Increment time and adapt timestep
  time = time+dt;
end
return
