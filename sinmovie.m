%% Sinmovie.m 
% to plot an animation of a sin function with increasing frequency.
% by Manuel Diaz, 2012.05.09

%% Clear: command pront, variables in memory and close any open figure.
clear; clc; close all;
 
t = linspace(0, 2*pi,1000); % creates 1000 points between 0 and 2*pi()

count = 1;
figure
for freqrps=0.1:0.1:2*pi
          y=sin(freqrps*t);
          plot(t,y);
          xlabel('x');
          ylabel('y');
          axis([0,2*pi,-1,1])
          S1=sprintf('y(t)=sin(%.2f t)',freqrps);
          text(2,.6,S1)
          freqcps=freqrps/(2*pi);
          S2=sprintf('frequency=%.2f rads/sec (%.2f cyc/sec)',freqrps,freqcps);
          text(2,.4,S2)
          title('Sinusoidal Function');
          M(count)=getframe;
          count=count+1;
end
 
movie(M,2,10);