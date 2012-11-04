function [leg_tb1,Dleg_tb1] = MLPN (N_max,LGLCoord,DOF)

for i = 0:N_max
    x = LGLCoord(i+1,N_max+1);
    [PN,PD] = LPN (DOF,x) ;
    leg_tb1(i+1,:) = PN ;
    Dleg_tb1(i+1,:) = PD ;
end


