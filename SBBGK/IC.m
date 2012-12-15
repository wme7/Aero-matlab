RL=1.0;
UL=0.75;
PL=1.0;

ET=PL+0.5*RL*UL^2
TL=4*ET/RL-2*UL^2
ZL=RL/sqrt(pi*TL)

RR=0.125;
UR=0;
PR=0.1;

ET=PR+0.5*RR*UR^2
TR=4*ET/RR-2*UR^2
ZR=RR/sqrt(pi*TR)