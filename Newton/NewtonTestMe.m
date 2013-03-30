% This MATLAB script runs multiple tests for the 
% function [XIter,XHistory,FXHistory,failureFlag] = Newton(X,F,X0,varargin)
% Newton's method for solving a system of nonlinear equations
% Newton(X,F,X0) solves nonlinear system F(x)=0 by Newton's method, using
% the given initial approximation X0

% License: BSD 
% Copyright (c) 2011-2012 A.V. Knyazev, Andrew.Knyazev@ucdenver.edu
% http://math.ucdenver.edu/~aknyazev/
% $Revision: 1.2 $  $Date: 6-Jan-2012
% Tested in MATLAB 7.13 (R2011b) and its Symbolic Math Toolbox 5.7 (R2011b) 
% Does NOT work in Octave 3.4.2 and below because of poor symbolic toolbox

clear all; close all;
link = 'http://en.wikipedia.org/wiki/Newtons_method';
linktitle = 'Newton''s method';
fprintf('This is a demo of <a href = "%s">%s</a>\n', link, linktitle) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('A scalar example, F=sin(X); X0=1; resulting in convergence');
disp('to the root X=0.');
clear all; close all; X=sym('X'); F=sin(X); X0=1;
Sol=Newton(X,F,X0); disp(Sol); % is equivalent to
Sol=Newton(X,F,X0,struct('tol',1e-15,'maxIter',5)); disp(Sol); % also equivalent to
opts.tol=1e-15; opts.maxIter=5; Sol=Newton(X,F,X0,opts); disp(Sol);
disp('The same example, only using different notaion for the input');
clear all; t=sym('t'); f=sin(t); t0=1;
[Sol,XHistory,FXHistory,failureFlag]=Newton(t,f,t0); 
fprintf('The failureFlag is %d, i.e. the method has converged.\n', failureFlag);
fprintf('The approximate solution is %g.\n',Sol); 
disp('Plotting the symbolic function,');
plotInterval=[min(XHistory) max(XHistory)];
ezplot(f,plotInterval); hold on; grid on;
disp('the numeric approximations to the root,');
plot(XHistory,FXHistory,'rp','MarkerSize',15);
title('Newton"s method approximations'); 
disp('and the corresponding lines.');
plot(plotInterval,[0 0],'-m');  % Plotting the x-axis
for j=1:length(XHistory)-1,
    plot([XHistory(j) XHistory(j+1)],[FXHistory(j) 0],'-g');
    plot([XHistory(j+1) XHistory(j+1)],[FXHistory(j+1) 0],'--g');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' '); disp('Press any key to continue to the next example'); pause;
disp('A scalar example, F=sin(X); X0=1.16556; resulting in convergence');
disp('to the root X=0.');
clear all; close all; t=sym('t'); f=sin(t); t0=1.16556;
[Sol,~,~,failureFlag]=Newton(t,f,t0); 
fprintf('The failureFlag is %d, i.e. the method has NOT converged.\n', ...
    failureFlag);
fprintf('The approximate solution is %g, far away from the  root X=0.\n',Sol);
disp('Let us try to increase the maxIter parameter from its defaut value 5 to 20');
[Sol,XHistory,FXHistory,failureFlag]=Newton(t,f,t0,struct('maxIter',20)); 
fprintf('The failureFlag is %d, i.e. the method has converged.\n', failureFlag);
fprintf('The approximate solution is %g.\n',Sol); 
disp('Plotting the symbolic function,');
plotInterval=[min(XHistory) max(XHistory)];
ezplot(f,plotInterval); hold on; grid on;
disp('the numeric approximations to the root,');
plot(XHistory,FXHistory,'rp','MarkerSize',15);
title('Newton"s method approximations'); 
disp('and the corresponding lines.');
plot(plotInterval,[0 0],'-m');  % Plotting the x-axis
for j=1:length(XHistory)-1,
    plot([XHistory(j) XHistory(j+1)],[FXHistory(j) 0],'-g');
    plot([XHistory(j+1) XHistory(j+1)],[FXHistory(j+1) 0],'--g');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' '); disp('Press any key to continue to the next example'); pause;
disp('A scalar example, F=sin(X); X0=1.165562; resulting in convergence');
disp('to the root X=pi.');
clear all; close all; t=sym('t'); f=sin(t); t0=1.165562;
[Sol,XHistory,FXHistory,failureFlag]=Newton(t,f,t0,struct('maxIter',20)); 
fprintf('The failureFlag is %d, i.e. the method has converged.\n', failureFlag);
fprintf('The approximate solution is %g.\n',Sol); 
disp('Plotting the symbolic function,');
plotInterval=[min(XHistory) max(XHistory)];
ezplot(f,plotInterval); hold on; grid on;
disp('the numeric approximations to the root,');
plot(XHistory,FXHistory,'rp','MarkerSize',15);
title('Newton"s method approximations'); 
disp('and the corresponding lines.');
plot(plotInterval,[0 0],'-m');  % Plotting the x-axis
for j=1:length(XHistory)-1,
    plot([XHistory(j) XHistory(j+1)],[FXHistory(j) 0],'-g');
    plot([XHistory(j+1) XHistory(j+1)],[FXHistory(j+1) 0],'--g');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' '); disp('Press any key to continue to the next example'); pause;
link = 'http://en.wikipedia.org/wiki/Newtons_method#Starting_point_enters_a_cycle';
linktitle = 'Newton''s method cycle';
fprintf('Finding a simple <a href = "%s">%s</a> for F=2+X*(-2+X*X).\n', ...
    link, linktitle) 
disp('A simple cycle of Newton''s method is when after two steps'); 
disp('the iterative vector is back to the starting point');
clear all; close all; X=sym('X'); F=2+X*(-2+X*X);  H = jacobian(F,X); 
N=X-H\F;  NN=compose(N,N); loopPoints=double(solve(NN-X)); 
fprintf('The solution %g is stationary, so also a cycle point. \n',...
    loopPoints(7));
fprintf('One pair of cycle points is %g abd %g \n',...
    loopPoints(5),loopPoints(6));
fprintf('The other pair of cycle points is %g abd %g \n',...
    loopPoints(1),loopPoints(2));
disp('There are also two other pairs, complex conjugate, not diplayed here.');
disp('Let us verify the first pair of cycle points by actual calculations.');
disp('A scalar example, F=2+X*(-2+X*X); X0=1; resulting in a cycle');
clear all; close all; X=sym('X'); F=2+X*(-2+X*X); X0=1; 
[~,XHistory,FXHistory]=Newton(X,F,X0);  
disp('Plotting the symbolic function,');
plotInterval=[-2 2];
ezplot(F,plotInterval); hold on; grid on;
disp('the numeric approximations to the root,');
plot(XHistory,FXHistory,'rp','MarkerSize',15);
title('Newton"s method approximations'); 
disp('and the corresponding lines.');
plot(plotInterval,[0 0],'-m');  % Plotting the x-axis
for j=1:length(XHistory)-1,
    plot([XHistory(j) XHistory(j+1)],[FXHistory(j) 0],'-g');
    plot([XHistory(j+1) XHistory(j+1)],[FXHistory(j+1) 0],'--g');
end
%
disp(' ');
disp('Press any key to continue to a slideshow for the same example'); pause;
plotInterval=[-2 2]; 
hf=figure; axis tight; set(gca,'NextPlot','replacechildren');
for j=1:length(XHistory)-1,
    plot(plotInterval,[0 0],'-m'); hold on; % Plotting the x-axis
    ezplot(F,plotInterval);   % Plotting the symbolic function
    plot(XHistory(j),FXHistory(j),'rp','MarkerSize',15); 
    plot([XHistory(j) XHistory(j+1)],[FXHistory(j) 0],'-g');
    plot([XHistory(j+1) XHistory(j+1)],[FXHistory(j+1) 0],'--g');
    title('Newton"s method approximations'); 
    drawnow; hold off; 
    M(j) = getframe(gcf); %captures entire interior of the current figure 
end
disp('Now you have a slide-show in MATLAB format that can be played in MATLAB');
disp('Press any key to play the MATLAB slide-show'); pause; 
hf = figure; axis off; movie(hf,M,1,1); %play one time with 1 frames/second  
disp('But wait, you want to send it to a friend, who has no MATLAB?');
disp('Let us use  MOVIE2AVI function for the conversion');
movie2avi(M,'MyCoolSlideShowNewton.avi','fps',1,'compression','None');
disp('With a file browser in your current directory'); 
disp(pwd) % shows the location of the current directory
disp('locate the file called MyCoolSlideShowNewton.avi'); 
disp('It will play in most video players.');
disp(' ');
disp('It is also possible to produce the AVI-format slide-show directly.');
disp('Press any key and be patient until your see "DONE"'); 
pause; close all; 
aviobj=avifile('MyCoolSlideShowNewtonDirect.avi','fps',1,'compression','None'); 
hf = figure('visible','off'); %turns visibility of figure off
for j=1:length(XHistory)-1,
    plot(plotInterval,[0 0],'-m'); hold on; % Plotting the x-axis
    ezplot(F,plotInterval);   % Plotting the symbolic function
    plot(XHistory(j),FXHistory(j),'rp','MarkerSize',15); 
    plot([XHistory(j) XHistory(j+1)],[FXHistory(j) 0],'-g');
    plot([XHistory(j+1) XHistory(j+1)],[FXHistory(j+1) 0],'--g');
    title('Newton"s method approximations'); 
    hold off;
    aviobj=addframe(aviobj,hf); %adds frames to the AVI file
end
aviobj=close(aviobj); %closes the AVI file
close(hf); %closes the handle to invisible figure
disp('DONE! With a file browser in your current directory'); 
disp(pwd) % shows the location of the current directory
disp('locate the file called MyCoolSlideShowNewtonDirect.avi'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' '); disp('Press any key to continue'); pause;
disp('Do you like movies? Let us make a cool movie!'); 
disp('A scalar example, F=2+X*(-2+X*X); X0=1; results in a cycle');
disp('as we have already observed. What if we slightly change X0?');
clear all; close all; X=sym('X');  F=2+X*(-2+X*X); 
opts.maxIter=10;
nFrames=50; startingX0=.99; finishingX0=.97; 
hf=figure; axis tight; set(gca,'NextPlot','replacechildren');
stepSize=(finishingX0-startingX0)/(nFrames-1); plotInterval=[-2 2];
for k=1:nFrames,
    X0=startingX0+stepSize*(k-1);
    [~,XHistory,FXHistory]=Newton(X,F,X0,opts);
    plot(plotInterval,[0 0],'-m');  % Plotting the x-axis
    hold on; grid on; ezplot(F,plotInterval); 
    plot(XHistory,FXHistory,'rp','MarkerSize',15);
    title('Newton"s method approximations');
    for j=1:length(XHistory)-1,
        plot([XHistory(j) XHistory(j+1)],[FXHistory(j) 0],'-g');
        plot([XHistory(j+1) XHistory(j+1)],[FXHistory(j+1) 0],'--g');
    end%for 
    drawnow; hold off; 
    M(k) = getframe(gcf); %captures entire interior of the current figure
end%for
disp('Now you have a movie in MATLAB format that can be played in MATLAB');
disp('Press any key to play the MATLAB movie'); pause; 
hf = figure; axis off; movie(hf,M,1,5); %play one time with 5 frames/second  
disp('But wait, you want to send it to a friend, who has no MATLAB?');
disp('Let us use  MOVIE2AVI function for the conversion');
movie2avi(M,'MyCoolMovieNewton.avi','fps',5,'compression','None');
disp('With a file browser in your current directory'); 
disp(pwd) % shows the location of the current directory
disp('locate the file called MyCoolMovieNewton.avi'); 
disp('It will play in most video players.');
disp(' ');
disp('It is also possible to produce the AVI-format movie directly.');
disp('Press any key and be patient until your see "DONE"'); 
pause; close all; 
aviobj=avifile('MyCoolMovieNewtonDirect.avi','fps',5,'compression','None'); 
hf = figure('visible','off'); %turns visibility of figure off
for k=1:nFrames,
    X0=startingX0+stepSize*(k-1);
    [~,XHistory,FXHistory]=Newton(X,F,X0,opts);
    plot(plotInterval,[0 0],'-m');  % Plotting the x-axis
    hold on; grid on; ezplot(F,plotInterval); 
    plot(XHistory,FXHistory,'rp','MarkerSize',15);
    title('Newton"s method approximations');
    for j=1:length(XHistory)-1,
        plot([XHistory(j) XHistory(j+1)],[FXHistory(j) 0],'-g');
        plot([XHistory(j+1) XHistory(j+1)],[FXHistory(j+1) 0],'--g');
    end%for 
    hold off;
    aviobj=addframe(aviobj,hf); %adds frames to the AVI file
end%for
aviobj=close(aviobj); %closes the AVI file
close(hf); %closes the handle to invisible figure
disp('DONE! With a file browser in your current directory'); 
disp(pwd) % shows the location of the current directory
disp('locate the file called MyCoolMovieNewtonDirect.avi'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Press any key to continue to the next example'); pause;
disp('A scalar example, approximating the +infinity');
clear all; close all; X=sym('X'); F=inv(X); X0=10; 
[~,XHistory,FXHistory,~]=Newton(X,F,X0); 
disp('Plotting the symbolic function');
plotInterval=[min(XHistory) max(XHistory)];
ezplot(F,plotInterval); hold on; grid on;
disp('and the numeric approximations to the root.');
plot(XHistory,FXHistory,'rp','MarkerSize',15);
title('Newton"s method approximations'); 
disp('and the corresponding lines.');
plot(plotInterval,[0 0],'-m');  % Plotting the x-axis
for j=1:length(XHistory)-1,
    plot([XHistory(j) XHistory(j+1)],[FXHistory(j) 0],'-g');
    plot([XHistory(j+1) XHistory(j+1)],[FXHistory(j+1) 0],'--g');
end%for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Press any key to continue to the next example'); pause;
disp('A scalar example, resulting in stagnation near the solution');
clear all; close all; X=sym('X'); F=(exp(X)-exp(100))*(X-100)^3; X0=101;
Sol=Newton(X,F,X0,struct('maxIter',200));
fprintf('The absolute error %g is not getting any smaller. \n',...
    abs(Sol-100));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples below need more work for better presentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Press any key to continue to the next example'); pause;
disp('More equations (2) than unknowns (1), system consistent');
clear all; close all; X=sym('X'); F=[sin(X); tan(X)]; X0=1;
[Sol,XHistory,FXHistory]=Newton(X,F,X0)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Press any key to continue to the next example'); pause;
disp('More equations (2) than unknowns (1), system inconsistent');
clear all; close all; X=sym('X'); F=[sin(X); tan(X)-1]; X0=1;
[Sol,XHistory,FXHistory]=Newton(X,F,X0)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Press any key to continue to the next example'); pause;
disp('2 equations and 2 unknowns, system consistent');
clear all; close all; X=sym('X',[2 1]); 
F=[sin(X(1)); tan(X(2))-1]; X0=ones(2,1); 
[Sol,XHistory,FXHistory]=Newton(X,F,X0) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Press any key to continue to the next example'); pause;
disp('2 equations and 2 unknowns, system consistent');
clear all; close all; X=sym('X',[2 1]); X0=2*ones(2,1);
F=[6*X(1)^3+X(1)*X(2)-3*X(2)^3-4; ...
   X(1)^2-18*X(1)*X(2)^2+16*X(2)^3+1]; 
[Sol,XHistory,FXHistory]=Newton(X,F,X0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('The end. Thank you for running me!');