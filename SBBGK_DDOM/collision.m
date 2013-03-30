function [u_next] = collision(u,fu,tau,dt)
% Compute BGK approximate collision operator for Boltzmann Equation.

u_next = u - dt/tau*(fu - u);
        
return