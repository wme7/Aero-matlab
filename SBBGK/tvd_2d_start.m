% TVD 2D SBBGK

    NX          = 399;
    NY          = 399;
    CFL         = 0.4;
    OUTTIME     = 0.35;
    %sv_point1 = (1./3.)*outtime
    %sv_point2 = (2./3.)*outtime 
    theta       = 1.;
    %! MB = 0, FD = 1, BE = -1
    NV = 20;
    %PI = ATAN2(1.,1.)*4.;

%!
%! GAUSS-HERMITE QUADRATURE
%!
GH= [-5.38748089001 -4.60368244955 -3.94476404012 -3.34785456738 ...
    -2.78880605843 -2.25497400209 -1.73853771212 -1.2340762154 ...
    -0.737473728545 -0.245340708301 0.245340708301 0.737473728545 ...
    1.2340762154 1.73853771212 2.25497400209 2.78880605843 ...
    3.34785456738 3.94476404012 4.60368244955 5.38748089001];
W = [0.898591961453 0.704332961176 0.62227869619 0.575262442852 ...
    0.544851742366 0.524080350949 0.509679027117 0.499920871336 ...
    0.493843385272 0.490921500667 0.490921500667 0.493843385272 ...
    0.499920871336 0.509679027117 0.524080350949 0.544851742366 ...
    0.575262442852 0.62227869619 0.704332961176 0.898591961453];

    DO I = 1, NV
        C(I)    = W(I)
        V(I)    = GH(I)
        WRITE (10,*) C(I),V(I)    
    END DO

    Z1   = 0.02
    UX1  = -0.5
    UY1  = -0.33
    T1   = 5.
    
    Z2   = 0.084
    UX2  = -0.5
    UY2  = 0.33
    T2   = 2.453
    
    Z3   = 0.02
    UX3  = 0.5
    UY3  = 0.33
    T3   = 5.
    
    Z4   = 0.2048
    UX4  = 0.5
    UY4  = -0.33
    T4   = 1.5934
    
    DX = 1./DFLOAT(NX-1)
    DY = 1./DFLOAT(NY-1)
    X(1)  = 0
    Y(1)  = 0
    DO I = 2,NX
        X(I) = X(I-1) + DX
    END DO
    DO J = 2, NY
        Y(J) = Y(J-1) + DY   
    END DO

DO I = 1,NX
DO J = 1,NY    
           if ((Y(J) .lt. 0.5) .and. (x(i) .gt. 0.5))  then
                Z(I,J) = Z4
                UX(I,J)= UX4
                UY(I,J)= UY4
                T(I,J) = T4   
            else if ((Y(J) .gt. 0.5) .and. (x(i) .lt. 0.5)) then
                Z(I,J) = Z2
                UX(I,J)= UX2
                UY(I,J)= UY2
                T(I,J) = T2
            else if ((Y(J) .ge. 0.5) .and. (x(i) .ge. 0.5)) then
                Z(I,J) = Z1
                UX(I,J)= UX1
                UY(I,J)= UY1
                T(I,J) = T1
            else
                Z(I,J) = Z3
                UX(I,J)= UX3
                UY(I,J)= UY3
                T(I,J) = T3    
            end if                
               
    DO K = 1, NV
    DO L = 1, NV
        PP = ( (V(K)-UX(I,J))^2 + (V(L)-UY(I,J))^2 ) / T(I,J)
        F(k,l,i,j) = 1/((EXP(PP)/Z(I,J)) + theta)
    END DO
    END DO    
END DO
END DO

do k = 1,nv
do l = 1,nv
    write (90,*) v(k),v(l),f(k,l,nx*3/4,ny*3/4),f(k,l,nx/4,ny*3/4),f(k,l,nx/4,ny/4),f(k,l,nx*3/4,ny/4)
end do
end do

ITER  = 1
TIME  = 0
ISTOP = 0
1000 CONTINUE
    
    DT = min(DX,DY) * CFL/V(NV)
    TIME = TIME + DT
    DTDX = DT/DX
    DTDY = DT/DY
    
IF (TIME .GT. OUTTIME) THEN
    DTCFL = OUTTIME - (TIME - DTCFL)  
    TIME = OUTTIME
    DT = DTCFL/CFL
    DTDX = DTCFL / DX
    DTDY = DTCFL / DY
    ISTOP = 1
END IF

DO I = 1, NX
DO J = 1, NY
    DO K = 1, NV
    DO L = 1, NV
        PP = ( (V(K)-UX(I,J))^2 + (V(L)-UY(I,J))^2 ) / T(I,J)
        F(k,l,i,j) = 1/((EXP(PP)/Z(I,J))+theta)
    END DO
    END DO  
END DO
END DO

DO K = 1, NV
DO L = 1, NV
    vxp = max(v(k),0.)
    vxm = min(v(k),0.) 
    vyp = max(v(l),0.)
    vym = min(v(l),0.) 
    %!theta
    do i = 2, nx-1
    do j = 2, ny-1
        if (f(k,l,i+1,j) .eq. f(k,l,i,j)) then
            tha(i,j) =  0
        else
            tha(i,j) = (f(k,l,i-sign(1.,v(k))+1,j) - f(k,l,i-sign(1.,v(k)),j))/(f(k,l,i+1,j)-f(k,l,i,j))
        end if 
        if (f(k,l,i,j+1) .eq. f(k,l,i,j)) then
            thb(i,j) =  0
        else
            thb(i,j) = (f(k,l,i,j-sign(1.,v(l))+1) - f(k,l,i,j-sign(1.,v(l))))/(f(k,l,i,j+1)-f(k,l,i,j))
        end if         
    end do
    end do 
       
    do i = 1, nx
        thb(i,1)=1
        thb(i,ny)=1
    end do
    
    do j = 1, ny
        tha(1,j)=1
        tha(nx,j)=1    
    end do
    %!Van Leer
    do i = 1,nx
    do j = 1,ny
        if (tha(i,j) .le. 0) then
            phia(i,j)=0
        else
            phia(i,j)= (abs(tha(i,j))+tha(i,j))/(1+abs(tha(i,j)))
        end if
        if (thb(i,j) .le. 0) then
            phib(i,j)=0
        else
            phib(i,j)= (abs(thb(i,j))+thb(i,j))/(1+abs(thb(i,j)))
        end if    
    end do
    end do
     
    DO I = 2, NX-1
    DO J = 2, NY-1 
        fll(i,j) = vxp*f(k,l,i-1,j)+vxm*f(k,l,i,j)
        flh(i,j) = 0.5d0*v(k)*(f(k,l,i-1,j)+f(k,l,i,j))-0.5d0*dtdx*v(k)^2*(f(k,l,i,j)-f(k,l,i-1,j)) 
        fl(i,j) = fll(i,j)+ (phia(i-1,j)*(flh(i,j)-fll(i,j)))
        
        fhl(i,j) = vxp*f(k,l,i,j)+vxm*f(k,l,i+1,j)
        fhh(i,j) = 0.5d0*v(k)*(f(k,l,i,j)+f(k,l,i+1,j))-0.5d0*dtdx*v(k)^2*(f(k,l,i+1,j)-f(k,l,i,j)) 
        fh(i,j)  = fhl(i,j)+ (phia(i,j)*(fhh(i,j)-fhl(i,j)))
    
        gll(i,j) = vyp*f(k,l,i,j-1)+vym*f(k,l,i,j)
        glh(i,j) = 0.5d0*v(l)*(f(k,l,i,j-1)+f(k,l,i,j))-0.5d0*dtdy*v(l)^2*(f(k,l,i,j)-f(k,l,i,j-1))
        gl(i,j)  = gll(i,j)+ (phib(i,j-1)*(glh(i,j)-gll(i,j)))
    
        ghl(i,j) = vyp*f(k,l,i,j)+vym*f(k,l,i,j+1)
        ghh(i,j) = 0.5d0*v(l)*(f(k,l,i,j)+f(k,l,i,j+1))-0.5d0*dtdy*v(l)^2*(f(k,l,i,j+1)-f(k,l,i,j))
        gh2(i,j)  = ghl(i,j)+ (phib(i,j)*(ghh(i,j)-ghl(i,j)))
    
        fn(i,j) = f(k,l,i,j) - dtdx*(fh(i,j)-fl(i,j)) - dtdy*(gh2(i,j)-gl(i,j))
    
        f(k,l,i,j) = fn(i,j)
        f(k,l,i,1) = fn(i,3)
        f(k,l,i,ny)= fn(i,ny-2)
        f(k,l,1,j) =fn(3,j)
        f(k,l,nx,j)=fn(nx-2,j)   
    END DO
    END DO  
END DO
END DO

DO I = 1, NX
DO J = 1, NY
    SR  = 0
    SUX = 0
    SUY = 0
    SE  = 0
      DO K = 1, NV
      DO L = 1, NV
        SR  = SR + C(K)*C(L) * F(k,l,i,j)
        SUX = SUX + C(K)*C(L) * F(k,l,i,j) * V(K)
        SUY = SUY + C(K)*C(L) * F(k,l,i,j) * V(L)
        SE  = SE + C(K)*C(L) * F(k,l,i,j) * (0.5 * (V(K)*V(K) + V(L)*V(L)))
      END DO
      END DO
    R(I,J)    = SR
    UX(I,J)   = SUX/SR 
    UY(I,J)   = SUY/SR
    ET(I,J)   = SE       
END DO
END DO

if (theta .eq. 0.) go to 1100
      
DO I = 1, NX
DO J = 1, NY
  ZA = 0.001
  ZB = 0.9
  DO WHILE (ABS(ZA-ZB) .GT. 0.0001)
    GA1 = 0
    GB1 = 0
    GA2 = 0
    GB2 = 0
        DO L = 1, 50
            if (theta .eq. 1.) then      
            GA1 = GA1 + (ZA^L) * (-1)^(L-1)/L
            GB1 = GB1 + (ZB^L) * (-1)^(L-1)/L
            GA2 = GA2 + (ZA^L) * (-1)^(L-1)/(L^2)
            GB2 = GB2 + (ZB^L) * (-1)^(L-1)/(L^2)
            else
            GA1 = GA1 + (ZA^L) /L
            GB1 = GB1 + (ZB^L) /L
            GA2 = GA2 + (ZA^L) /(L^2)
            GB2 = GB2 + (ZB^L) /(L^2)
            end if    
        END DO
    PSIA = 2*ET(I,J) - (GA2*(R(I,J)/GA1)^2)/PI - R(I,J)*(UX(I,J)*UX(I,J)+UY(I,J)*UY(I,J))
    PSIB = 2*ET(I,J) - (GB2*(R(I,J)/GB1)^2)/PI - R(I,J)*(UX(I,J)*UX(I,J)+UY(I,J)*UY(I,J))
    ZC = (ZA+ZB)/2
    GC1 = 0
    GC2 = 0
        DO L = 1, 50
            if (theta .eq. 1.) then
            GC1 = GC1 + (ZC^L) * (-1)^(L-1)/L
            GC2 = GC2 + (ZC^L) * (-1)^(L-1)/(L^2)
            else 
            GC1 = GC1 + (ZC^L)/L
            GC2 = GC2 + (ZC^L)/(L^2)
            end if
        END DO
    PSIC = 2*ET(I,J) - (GC2*(R(I,J)/GC1)^2)/PI - R(I,J)*(UX(I,J)*UX(I,J)+UY(I,J)*UY(I,J))
    IF ((PSIA*PSIC) .LT. 0) THEN
        ZB = ZC
    ELSE
        ZA = ZC
    END IF
  END DO
  Z(I,J) = ZC
  T(I,J) = R(I,J)/(PI*GC1)
  P(I,j) = ET(i,j) - 0.5*r(i,j)*(ux(i,j)^2+uy(i,j)^2)
END DO
END DO 
go to 1300

1100 continue
DO I = 1, NX
DO J = 1, NY
  T(I,J) = (2* ET(I,J)/R(I,J))-(UX(I,J)^2+Uy(I,J)^2)
  Z(I,J) = R(I,J)/(pi*T(i,j))
  P(I,j) = r(i,j)*T(i,j)/2
END DO
END DO    

1300 continue
DO I = 1, NX
DO J = 1, NY
   if (iter.eq.floor(sv_point1/DT)) then
     WRITE (20,*) X(I), Y(J), r(I,J), p(I,J), z(I,J)
   else if (iter.eq.floor(sv_point2/dt)) then
     WRITE (30,*) X(I), Y(J), r(I,J), p(I,J), z(I,J)
   end if
END DO 
END DO
     
WRITE (*,*) ITER,time, R(NX/2,NY/2), F(NV/2,NV/2,NX/2,NY/2)
IF (ISTOP .EQ. 1) GOTO 2000
ITER = ITER + 1
GOTO 1000

2000 CONTINUE
DO I = 1, NX
DO J = 1, NY
   WRITE (60,*) X(I), Y(J), r(I,J), p(I,J), z(I,J)
END DO 
END DO
