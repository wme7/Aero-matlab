%--------------------------------------------------------------------------
%TestScript02 for Collision Ball class
%--------------------------------------------------------------------------

%CodeStart-----------------------------------------------------------------
%Resetting MATLAB environment
    close all
    clear
    clc
%Creating Balls object
    ball=Balls();
%Adding 8 big balls as wall
    for i=1:9
    for j=1:9
        ball.addBall((20*i)-100,(20*j)-100,10*(i-5),10*(j-5),10,5,[0,0,255]);
    end
    end
%Animating ball
    ball.play(0.05)
%CodeEnd-------------------------------------------------------------------