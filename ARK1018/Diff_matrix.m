function [Luq] = Diff_matrix(pp,NN,r,dx)
width_Q = pp*NN;
Lu = zeros(width_Q); % do differential on q
Lq = Lu;             % do differential on u
Luq = Lu;

[Au,Bu,Cu,Aq,Bq,Cq] = Flux(pp,dx);
r_sqrt = r^(1/2);

for ik=2:NN-1
    for i=1:pp
        for jk=1:NN
            for j=1:pp
                if (ik-2)*pp+0<j+(jk-1)*pp && j+(jk-1)*pp<(ik-2)*pp+pp+1
                    Lq(i+(ik-1)*pp,j+(jk-1)*pp) = r_sqrt*Aq(i,j);
                    Lu(i+(ik-1)*pp,j+(jk-1)*pp) = r_sqrt*Au(i,j);
                elseif (ik-2)*pp+pp<j+(jk-1)*pp && j+(jk-1)*pp<(ik-2)*pp+2*pp+1
                    Lq(i+(ik-1)*pp,j+(jk-1)*pp) = r_sqrt*Bq(i,j);
                    Lu(i+(ik-1)*pp,j+(jk-1)*pp) = r_sqrt*Bu(i,j);
                elseif (ik-2)*pp+2*pp<j+(jk-1)*pp && j+(jk-1)*pp<(ik-2)*pp+3*pp+1
                    Lq(i+(ik-1)*pp,j+(jk-1)*pp) = r_sqrt*Cq(i,j);
                    Lu(i+(ik-1)*pp,j+(jk-1)*pp) = r_sqrt*Cu(i,j);
                end
            end
        end
    end
end

Lq(1:pp,1:pp)       = r_sqrt*Bq;
Lq(1:pp,1+pp:pp+pp) = r_sqrt*Cq;      
Lu(1:pp,1:pp)       = r_sqrt*Bu;
Lu(1:pp,1+pp:pp+pp) = r_sqrt*Cu;
% {   
Lq(1:pp,(width_Q-pp+1):width_Q) = r_sqrt*Aq; % period B.C.
Lq((width_Q-pp+1):width_Q,1:pp) = r_sqrt*Cq; % period B.C.
Lu(1:pp,(width_Q-pp+1):width_Q) = r_sqrt*Au; % period B.C.
Lu((width_Q-pp+1):width_Q,1:pp) = r_sqrt*Cu; % period B.C.
%}     
Lq((width_Q-pp+1):width_Q,(width_Q-2*pp+1):(width_Q-pp)) = r_sqrt*Aq;     
Lq((width_Q-pp+1):width_Q,(width_Q-pp+1):width_Q)        = r_sqrt*Bq;      
Lu((width_Q-pp+1):width_Q,(width_Q-2*pp+1):(width_Q-pp)) = r_sqrt*Au;   
Lu((width_Q-pp+1):width_Q,(width_Q-pp+1):width_Q)        = r_sqrt*Bu;
%}
Luq = Lu*Lq;