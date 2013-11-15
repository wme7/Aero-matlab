%--------------------------------------------------------------------------
%TestScript04 for Collision Ball class
%--------------------------------------------------------------------------

%CodeStart-----------------------------------------------------------------
%Resetting MATLAB environment
    close all
    clear
    clc
%Creating Balls object
    ball=Balls();
%Adding lining ball
    for i=1:9
        ball.addBall(0,(21*i)-110,00,0,10,10,[0,0,255]);
    end
%Adding impactor
    ball.addBall(0,99,0,-10,10,1,[0,0,0]);
%Animating ball
    ball.play(0.05)
%CodeEnd-------------------------------------------------------------------