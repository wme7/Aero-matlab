function [u_next] = collision(u,u_eq,tau,dt)
% Compute BGK approximate collision operator for Boltzmann Equation.
u_next = u - dt/tau*(u - u_eq);