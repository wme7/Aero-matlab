for J = 1:NXP1
    SR = 0;
    SU = 0;
    SE = 0;
    SAV= 0;
    for K = 1:NV
        SR = SR + C(K) * F(K,J);
        SU = SU + C(K) * F(K,J) * V(K);
        SE = SE + C(K) * F(K,J) * (0.5 * V(K) * V(K));
        SAV = SAV + C(K) * F(K,J) * abs(V(K));
    end
    R(J)    = SR;
    U(J)    = (SU/SR);
    ET(J)   = SE;
    AV(J)   = SAV;
end