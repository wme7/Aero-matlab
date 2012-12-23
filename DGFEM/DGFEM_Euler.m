%**************************************************************************
%* --- Main program for the 1D Euler shock-tube solver.
%*
%* This code solves the Sod's shock tube problem which is described
%* in Section 7.13.3 of "I do like CFD, VOL.1": Sod's problem 1, Figure 7.12.2.
%*
%* - t=0                               - t=tf
%* Density                             Density
%*   ****************|                 *********\
%*                   |                           \
%*                   |                            \
%*                   |                             ****|
%*                   |                                 |
%*                   |                                 ****|
%*                   ***************                       ***********
%*
%* Methods employed:
%*   - Roe's flux
%*   - Minmod limiter
%*   - Two-stage Runge-Kutta time stepping scheme
%*
%* Input ---------------------------------------------------
%*    ncells = # of cells on a grid.
%*        tf = Final time
%*       cfl = CFL number (<1)
%*      xmin = Left end of the domain
%*      xmax = Right end of the domain
%* 
%* Output --------------------------------------------------
%*  "solution.dat" = Data file containing for each cell,
%*                   cell-center coordinate, density, velocity, pressure, 
%*                   entropy, in the following format:
%*
%*       write(*,*) ncells
%*      do i=1,ncells
%*       write(*,*) xc(i), rho(i), u(i), p(i), entropy(i)
%*      end do
%* 
%*     Use the matlab program, oned_euler_v1.m, to plot the solutions.
%*
%*
%*  Note: Explore other shock tube problems by changing the initial condition
%*        and the final time (Other problems are described in Section 7.13.3
%*        of "I do like CFD, VOL.1").
%*
%*  Note: Other limiters may be explored (see CFD textbooks).
%*
%*  Note: Other flux functions may be explored.
%*        Various subroutines are available at cfdbooks.com: Osher, Van Leer, etc.
%*
%*  Note: Boundary condition need to be modified for problems with hard
%*        boundaries.
%*
%*
%*  Based on Katate Masatsuka, December 2010. http://www.cfdbooks.com
%**************************************************************************
%% Clear Work Space
clc; clear all; close all;