% TANDEMO   This file generates the tan(x) - x example from Chapter 1.
%           See Figures 1.2 and 1.3.
%
% C. T. Kelley, October 22, 2002.
%
x0=4.5; tol=1.d-20;
%
% Solve the problem three times.
%
[x,hist]=newtsol(x0,'ftan',tol,tol,1);
[x,histc]=chordsol(x0,'ftan',tol,tol);
[x,hists]=secant(x0,'ftan',tol,tol);
%
% Iteration history for Newton. 
% I use handle graphics to get the axis the way I
% want them. Try commenting out "set(p1,'XTick',[0 1 2 3 4 5]);"
% and see what the difference is.
%
maxit=6;
figure(1);
p1=subplot(1,1,1);
semilogy(hist(1:maxit,1),abs(hist(1:maxit,2)))
set(p1,'XTick',[0 1 2 3 4 5],'FontSize',14);
axis([0 5 1.d-16 1]);
xlabel('Nonlinear iterations'); ylabel('Absolute Nonlinear Residual');

%
% Plot 15 iterations for all three methods.
%
figure(2);
maxit=15;
p2=subplot(1,1,1);
semilogy(hist(1:maxit,1),abs(hist(1:maxit,2)),'-',...
histc(1:maxit,1),abs(histc(1:maxit,2)),'--',...
hists(1:maxit,1),abs(hists(1:maxit,2)),'-.');
set(p2,'XTick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14],'FontSize',14);
axis([0 14 1.d-16 1]);
legend('Newton','Chord','Secant');
xlabel('Nonlinear iterations'); ylabel('Absolute Nonlinear Residual');
