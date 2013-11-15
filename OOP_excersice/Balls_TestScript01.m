%--------------------------------------------------------------------------
%TestScript01 for Collision Ball class
%--------------------------------------------------------------------------

%CodeStart-----------------------------------------------------------------
%Resetting MATLAB environment
    close all
    clear
    clc
%Creating Balls object
    ball=Balls();
%Drawing balls
    ball.draw();
%Randomly adding balls
    for i=1:50
        ball.addBall();
        pause(0.001)
    end
%Animating balls
    ball.play(0.05)
%CodeEnd-------------------------------------------------------------------