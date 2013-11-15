%--------------------------------------------------------------------------
%TestScript03 for Collision Balls class
%--------------------------------------------------------------------------

%CodeStart-----------------------------------------------------------------
%Resetting MATLAB environment
    close all
    clear
    clc
%Creating Balls object
    ball=Balls();
%Adding 10 big balls as wall
    for i=1:10
        ball.addBall(0,(20*i)-110,0,0,...
                     100000000000000000000000000000,10,[0,0,0]);
    end
%Deleting balls to create hole in the wall
    ball.removeBall(8)
    ball.removeBall(3)
%Adding small lining balls
    for i=1:19
        ball.addBall(-90,(10*i)-100,100,0,...
                     5,5,[0,0,255]);
    end
%Adding some random balls
    for i=1:10
        ball.addBall()
    end
%Animating ball
    ball.play(0.05)
%CodeEnd-------------------------------------------------------------------