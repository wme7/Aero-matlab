program Square
implicit real(8) (a-h,o-z)
real, dimension (20,20,400,400)     :: f,feq
real, dimension (400,400)           :: ux,uy,z,t,r,et,p,fn,tha,thb,phia,phib
real, dimension (400,400)           :: fll,flh,fl,fhl,fhh,fh,gll,glh,gl,ghl,ghh,gh2
real, dimension (400)               :: x,y
real, dimension (20)                :: c,v,gh,w
    nx          = 40
    ny          = 40
    cfl         = 0.45
    outtime     = 5.3
    sv_point1 = (1./5.)*outtime   
    sv_point2 = (2./5.)*outtime
	sv_point2b = (5./10.)*outtime
	sv_point3 = (3./5.)*outtime
	sv_point3b = (7./10.)*outtime  
	sv_point4 = (4./5.)*outtime 
	sv_point4b = (9./10.)*outtime
    theta       = 0
    ! mb = 0, fd = 1, be = -1
    nv = 20
	r_time = 0.1 
    pi = atan2(1.,1.)*4.
open (unit = 10, file = 'abscissas.tec', status = 'unknown')
open (unit = 15, file = 'initial.tec', status = 'unknown')
open (unit = 20, file = 'sv_point(1).tec', status = 'unknown')
open (unit = 25, file = 'sv_point(2).tec', status = 'unknown')
open (unit = 27, file = 'sv_point(2b).tec', status = 'unknown')	 !
open (unit = 30, file = 'sv_point(3).tec', status = 'unknown')
open (unit = 33, file = 'sv_point(3b).tec', status = 'unknown')	 !
open (unit = 35, file = 'sv_point(4).tec', status = 'unknown')
open (unit = 37, file = 'sv_point(4b).tec', status = 'unknown')
open (unit = 40, file = 'results.tec', status = 'unknown')

888 format ('variables = "x","y","n","p","z"')  
!
! gauss-hermite quadrature
!
gh(:) =(/-5.38748089001,-4.60368244955,-3.94476404012,-3.34785456738,-2.78880605843,-2.25497400209,-1.73853771212,-1.2340762154,-0.737473728545,-0.245340708301,0.245340708301,0.737473728545,1.2340762154,1.73853771212,2.25497400209,2.78880605843,3.34785456738,3.94476404012,4.60368244955,5.38748089001/)
w(:)  =(/0.898591961453,0.704332961176,0.62227869619,0.575262442852,0.544851742366,0.524080350949,0.509679027117,0.499920871336,0.493843385272,0.490921500667,0.490921500667,0.493843385272,0.499920871336,0.509679027117,0.524080350949,0.544851742366,0.575262442852,0.62227869619,0.704332961176,0.898591961453/)
    do i = 1, nv
        c(i)    = w(i)
        v(i)    = gh(i)
        write (10,*) c(i),v(i)    
    end do
!
    z1   = 0.078
    ux1  = 0.7156
    uy1  = 0.
    t1   = 1.2586
!    
    z2   = 0.1
    ux2  = 0.
    uy2  = 0.
    t2   = 0.5 
!
	zw   = z2
    uxw  = 0.
    uyw  = 0.
    tw   = t2
! 	    
    dx = 10./dfloat(nx-1)
    dy = 10./dfloat(ny-1)
    x(1)  = 0
    y(1)  = 0
    do i = 2,nx
        x(i) = x(i-1) + dx
    end do
    do j = 2, ny
        y(j) = y(j-1) + dy   
    end do 
!initial condition
do k = 1, nv
do l = 1, nv 
	do i = 1,nx
	do j = 1,ny    
           if  (x(i).le.1)then
                z(i,j) = z1
                ux(i,j)= ux1
                uy(i,j)= uy1
                t(i,j) = t1
				pp = ( (v(k)-ux(i,j))**2 + (v(l)-uy(i,j))**2 ) / t(i,j)
				f(k,l,i,j) = 1/((exp(pp)/z(i,j)) + theta)  
            else if ((x(i).le.4.5) .or. (x(i).ge.5.5) .or. (y(j).le.4.5) .or. (y(j).ge.5.5)) then
                z(i,j) = z2
                ux(i,j)= ux2
                uy(i,j)= uy2
                t(i,j) = t2
				pp = ( (v(k)-ux(i,j))**2 + (v(l)-uy(i,j))**2 ) / t(i,j)
				f(k,l,i,j) = 1/((exp(pp)/z(i,j)) + theta)
			else if ((x(i).le.(4.5+dx)) .or. (x(i).ge.(5.5-dx)) .or. (y(j).le.(4.5+dy)) .or. (y(j).ge.(5.5-dy)))	then
				z(i,j) = zw
                ux(i,j)= uxw
                uy(i,j)= uyw
                t(i,j) = tw
				pp = ( (v(k)-ux(i,j))**2 + (v(l)-uy(i,j))**2 ) / t(i,j)
				f(k,l,i,j) = 1/((exp(pp)/z(i,j)) + theta)
            end if
    end do
    end do    
end do
end do
!
iter  = 1
time  = 0
istop = 0
1000 continue
    dt = min(dx,dy) * cfl/v(nv)
    time = time + dt
    dtdx = dt/dx
    dtdy = dt/dy     
if (time .gt. outtime) then
    dtcfl = outtime - (time - dtcfl)  
    time = outtime
    dt = dtcfl/cfl
    dtdx = dtcfl / dx
    dtdy = dtcfl / dy
    istop = 1
end if		
! projection to equilibrium
do k = 1, nv
do l = 1, nv	
	do i = 1, nx
	do j = 1, ny
		if ((x(i).le.4.5) .or. (x(i).ge.5.5) .or. (y(j).le.4.5) .or. (y(j).ge.5.5)) then 
        pp = ( (v(k)-ux(i,j))**2 + (v(l)-uy(i,j))**2 ) / t(i,j)
        feq(k,l,i,j) = 1/((exp(pp)/z(i,j))+theta)
		f(k,l,i,j)=feq(k,l,i,j)	
		end if
    end do
    end do  
end do
end do
!
do k = 1, nv
do l = 1, nv
    vxp = max(v(k),0.)
    vxm = min(v(k),0.) 
    vyp = max(v(l),0.)
    vym = min(v(l),0.) 
    !theta
    do i = 2, nx-1
    do j = 2, ny-1
	if ((x(i).le.4.5) .or. (x(i).ge.5.5) .or. (y(j).le.4.5) .or. (y(j).ge.5.5)) then
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
    !van leer
    do i = 1,nx
    do j = 1,ny
	if ((x(i).le.4.5) .or. (x(i).ge.5.5) .or. (y(j).le.4.5) .or. (y(j).ge.5.5)) then
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
	end if   
    end do
    end do
! l part
	do i = 2,nx-1
	do j = 2,ny-1	
		if ((x(i).le.4.5+dx).and.(x(i).ge.4.5-dx).and.(y(j).le.5.5+dy) .and. (y(j).ge.4.5-dy)) then
		ppi = ( (v(k)+ux(i-1,j))**2 + (v(l)-uy(i-1,j))**2 ) / t(i-1,j)
        f(k,l,i,j) = 1/((exp(ppi)/z(i-1,j))+theta) 
		end if 
	end do
	end do
	do i = 2,nx-1
	do j = 2,ny-1	 	
		if (x(i).le.4.5) then	 
        fll(i,j) = vxp*f(k,l,i-1,j)+vxm*f(k,l,i,j)
        flh(i,j) = 0.5d0*v(k)*(f(k,l,i-1,j)+f(k,l,i,j))-0.5d0*dtdx*v(k)**2*(f(k,l,i,j)-f(k,l,i-1,j)) 
        fl(i,j) = fll(i,j)+ (phia(i-1,j)*(flh(i,j)-fll(i,j)))
        fhl(i,j) = vxp*f(k,l,i,j)+vxm*f(k,l,i+1,j)
        fhh(i,j) = 0.5d0*v(k)*(f(k,l,i,j)+f(k,l,i+1,j))-0.5d0*dtdx*v(k)**2*(f(k,l,i+1,j)-f(k,l,i,j)) 
        fh(i,j)  = fhl(i,j)+ (phia(i,j)*(fhh(i,j)-fhl(i,j)))
        gll(i,j) = vyp*f(k,l,i,j-1)+vym*f(k,l,i,j)
        glh(i,j) = 0.5d0*v(l)*(f(k,l,i,j-1)+f(k,l,i,j))-0.5d0*dtdy*v(l)**2*(f(k,l,i,j)-f(k,l,i,j-1))
        gl(i,j)  = gll(i,j)+ (phib(i,j-1)*(glh(i,j)-gll(i,j)))
        ghl(i,j) = vyp*f(k,l,i,j)+vym*f(k,l,i,j+1)
        ghh(i,j) = 0.5d0*v(l)*(f(k,l,i,j)+f(k,l,i,j+1))-0.5d0*dtdy*v(l)**2*(f(k,l,i,j+1)-f(k,l,i,j))
        gh2(i,j)  = ghl(i,j)+ (phib(i,j)*(ghh(i,j)-ghl(i,j)))
        fn(i,j) = f(k,l,i,j) - dtdx*(fh(i,j)-fl(i,j)) - dtdy*(gh2(i,j)-gl(i,j))+ dt*(feq(k,l,i,j)-f(k,l,i,j))/r_time
		end if
	end do
	end do
! b part
	do i = 2,nx-1
	do j = 2,ny-1
		if ((y(j).le.4.5+dy).and.(y(j).ge.4.5-dy).and.(x(i).ge.4.5-dx).and.(x(i).le.5.5+dx)) then
		ppii = ( (v(k)-ux(i,j-1))**2 + (v(l)+uy(i,j-1))**2 ) / t(i,j-1)
        f(k,l,i,j) = 1/((exp(ppii)/z(i,j-1))+theta)
		end if
	end do
	end do
	do i = 2,nx-1
	do j = 2,ny-1	 	
		if ((y(j).le.4.5).and.(x(i).ge.4.5-dx).and.(x(i).le.5.5+dx)) then	 
        fll(i,j) = vxp*f(k,l,i-1,j)+vxm*f(k,l,i,j)
        flh(i,j) = 0.5d0*v(k)*(f(k,l,i-1,j)+f(k,l,i,j))-0.5d0*dtdx*v(k)**2*(f(k,l,i,j)-f(k,l,i-1,j)) 
        fl(i,j) = fll(i,j)+ (phia(i-1,j)*(flh(i,j)-fll(i,j)))
        fhl(i,j) = vxp*f(k,l,i,j)+vxm*f(k,l,i+1,j)
        fhh(i,j) = 0.5d0*v(k)*(f(k,l,i,j)+f(k,l,i+1,j))-0.5d0*dtdx*v(k)**2*(f(k,l,i+1,j)-f(k,l,i,j)) 
        fh(i,j)  = fhl(i,j)+ (phia(i,j)*(fhh(i,j)-fhl(i,j)))
        gll(i,j) = vyp*f(k,l,i,j-1)+vym*f(k,l,i,j)
        glh(i,j) = 0.5d0*v(l)*(f(k,l,i,j-1)+f(k,l,i,j))-0.5d0*dtdy*v(l)**2*(f(k,l,i,j)-f(k,l,i,j-1))
        gl(i,j)  = gll(i,j)+ (phib(i,j-1)*(glh(i,j)-gll(i,j)))
        ghl(i,j) = vyp*f(k,l,i,j)+vym*f(k,l,i,j+1)
        ghh(i,j) = 0.5d0*v(l)*(f(k,l,i,j)+f(k,l,i,j+1))-0.5d0*dtdy*v(l)**2*(f(k,l,i,j+1)-f(k,l,i,j))
        gh2(i,j)  = ghl(i,j)+ (phib(i,j)*(ghh(i,j)-ghl(i,j)))
        fn(i,j) = f(k,l,i,j) - dtdx*(fh(i,j)-fl(i,j)) - dtdy*(gh2(i,j)-gl(i,j))+ dt*(feq(k,l,i,j)-f(k,l,i,j))/r_time
		end if
	end do
	end do
! t part
	do i = 2,nx-1
	do j = 2,ny-1
		if ((y(j).ge.5.5-dy).and.(y(j).le.5.5+dy).and.(x(i).ge.4.5-dx).and.(x(i).le.5.5+dx)) then
		ppiv = ( (v(k)-ux(i,j+1))**2 + (v(l)+uy(i,j+1))**2 ) / t(i,j+1)
        f(k,l,i,j) = 1/((exp(ppiv)/z(i,j+1))+theta)
		end if
	end do
	end do
	do i = 2,nx-1
	do j = 2,ny-1	 	
		if ((y(j).ge.5.5).and.(x(i).ge.4.5-dx).and.(x(i).le.5.5+dx)) then	 
        fll(i,j) = vxp*f(k,l,i-1,j)+vxm*f(k,l,i,j)
        flh(i,j) = 0.5d0*v(k)*(f(k,l,i-1,j)+f(k,l,i,j))-0.5d0*dtdx*v(k)**2*(f(k,l,i,j)-f(k,l,i-1,j)) 
        fl(i,j) = fll(i,j)+ (phia(i-1,j)*(flh(i,j)-fll(i,j)))
        fhl(i,j) = vxp*f(k,l,i,j)+vxm*f(k,l,i+1,j)
        fhh(i,j) = 0.5d0*v(k)*(f(k,l,i,j)+f(k,l,i+1,j))-0.5d0*dtdx*v(k)**2*(f(k,l,i+1,j)-f(k,l,i,j)) 
        fh(i,j)  = fhl(i,j)+ (phia(i,j)*(fhh(i,j)-fhl(i,j)))
        gll(i,j) = vyp*f(k,l,i,j-1)+vym*f(k,l,i,j)
        glh(i,j) = 0.5d0*v(l)*(f(k,l,i,j-1)+f(k,l,i,j))-0.5d0*dtdy*v(l)**2*(f(k,l,i,j)-f(k,l,i,j-1))
        gl(i,j)  = gll(i,j)+ (phib(i,j-1)*(glh(i,j)-gll(i,j)))
        ghl(i,j) = vyp*f(k,l,i,j)+vym*f(k,l,i,j+1)
        ghh(i,j) = 0.5d0*v(l)*(f(k,l,i,j)+f(k,l,i,j+1))-0.5d0*dtdy*v(l)**2*(f(k,l,i,j+1)-f(k,l,i,j))
        gh2(i,j)  = ghl(i,j)+ (phib(i,j)*(ghh(i,j)-ghl(i,j)))
        fn(i,j) = f(k,l,i,j) - dtdx*(fh(i,j)-fl(i,j)) - dtdy*(gh2(i,j)-gl(i,j))	+ dt*(feq(k,l,i,j)-f(k,l,i,j))/r_time
		end if
	end do
	end do
! r part
	do i = 2,nx-1
	do j = 2,ny-1
		if ((x(i).ge.5.5-dx).and.(x(i).le.5.5+dx).and.(y(j).le.5.5+dy) .and. (y(j).ge.4.5-dy)) then
		ppiii = ( (v(k)+ux(i+1,j))**2 + (v(l)-uy(i+1,j))**2 ) / t(i+1,j)
        f(k,l,i,j) = 1/((exp(ppiii)/z(i+1,j))+theta)
		end if
	end do
	end do
	do i = 2,nx-1
	do j = 2,ny-1	 	
		if (x(i).ge.5.5) then	 
        fll(i,j) = vxp*f(k,l,i-1,j)+vxm*f(k,l,i,j)
        flh(i,j) = 0.5d0*v(k)*(f(k,l,i-1,j)+f(k,l,i,j))-0.5d0*dtdx*v(k)**2*(f(k,l,i,j)-f(k,l,i-1,j)) 
        fl(i,j) = fll(i,j)+ (phia(i-1,j)*(flh(i,j)-fll(i,j)))
        fhl(i,j) = vxp*f(k,l,i,j)+vxm*f(k,l,i+1,j)
        fhh(i,j) = 0.5d0*v(k)*(f(k,l,i,j)+f(k,l,i+1,j))-0.5d0*dtdx*v(k)**2*(f(k,l,i+1,j)-f(k,l,i,j)) 
        fh(i,j)  = fhl(i,j)+ (phia(i,j)*(fhh(i,j)-fhl(i,j)))
        gll(i,j) = vyp*f(k,l,i,j-1)+vym*f(k,l,i,j)
        glh(i,j) = 0.5d0*v(l)*(f(k,l,i,j-1)+f(k,l,i,j))-0.5d0*dtdy*v(l)**2*(f(k,l,i,j)-f(k,l,i,j-1))
        gl(i,j)  = gll(i,j)+ (phib(i,j-1)*(glh(i,j)-gll(i,j)))
        ghl(i,j) = vyp*f(k,l,i,j)+vym*f(k,l,i,j+1)
        ghh(i,j) = 0.5d0*v(l)*(f(k,l,i,j)+f(k,l,i,j+1))-0.5d0*dtdy*v(l)**2*(f(k,l,i,j+1)-f(k,l,i,j))
        gh2(i,j)  = ghl(i,j)+ (phib(i,j)*(ghh(i,j)-ghl(i,j)))
        fn(i,j) = f(k,l,i,j) - dtdx*(fh(i,j)-fl(i,j)) - dtdy*(gh2(i,j)-gl(i,j)) + dt*(feq(k,l,i,j)-f(k,l,i,j))/r_time
		end if
	end do
	end do
	do i = 2, nx-1
    do j = 2, ny-1
	    f(k,l,i,j) = fn(i,j)
        f(k,l,i,1) = fn(i,2)
        f(k,l,i,ny)= fn(i,ny-1)
        f(k,l,1,j) =fn(2,j)
        f(k,l,nx,j)=fn(nx-1,j)
	end do
	end do
end do
end do
!
do i = 1, nx
do j = 1, ny
if ((x(i).le.4.5) .or. (x(i).ge.5.5) .or. (y(j).le.4.5) .or. (y(j).ge.5.5)) then
    sr  = 0
    sux = 0
    suy = 0
    se  = 0
      do k = 1, nv
      do l = 1, nv
        sr  = sr + c(k)*c(l) * f(k,l,i,j)
        sux = sux + c(k)*c(l) * f(k,l,i,j) * v(k)
        suy = suy + c(k)*c(l) * f(k,l,i,j) * v(l)
        se  = se + c(k)*c(l) * f(k,l,i,j) * (0.5 * (v(k)*v(k) + v(l)*v(l)))
      end do
      end do
    r(i,j)    = sr
    ux(i,j)   = sux/sr 
    uy(i,j)   = suy/sr
    et(i,j)   = se 
end if      
end do
end do
!
if (theta .eq. 0.) go to 1100
!      
do i = 1, nx
do j = 1, ny
if ((x(i).le.4.5) .or. (x(i).ge.5.5) .or. (y(j).le.4.5) .or. (y(j).ge.5.5)) then
  za = 0.01
  zb = 0.9
  do while (abs(za-zb) .gt. 0.0001)
    ga1 = 0
    gb1 = 0
    ga2 = 0
    gb2 = 0
        do l = 1, 50
            if (theta .eq. 1.) then      
            ga1 = ga1 + (za**l) * (-1)**(l-1)/l
            gb1 = gb1 + (zb**l) * (-1)**(l-1)/l
            ga2 = ga2 + (za**l) * (-1)**(l-1)/(l**2)
            gb2 = gb2 + (zb**l) * (-1)**(l-1)/(l**2)
            else
            ga1 = ga1 + (za**l) /l
            gb1 = gb1 + (zb**l) /l
            ga2 = ga2 + (za**l) /(l**2)
            gb2 = gb2 + (zb**l) /(l**2)
            end if    
        end do
    psia = 2*et(i,j) - (ga2*(r(i,j)/ga1)**2)/pi - r(i,j)*(ux(i,j)*ux(i,j)+uy(i,j)*uy(i,j))
    psib = 2*et(i,j) - (gb2*(r(i,j)/gb1)**2)/pi - r(i,j)*(ux(i,j)*ux(i,j)+uy(i,j)*uy(i,j))
    zc = (za+zb)/2
    gc1 = 0
    gc2 = 0
        do l = 1, 50
            if (theta .eq. 1.) then
            gc1 = gc1 + (zc**l) * (-1)**(l-1)/l
            gc2 = gc2 + (zc**l) * (-1)**(l-1)/(l**2)
            else 
            gc1 = gc1 + (zc**l)/l
            gc2 = gc2 + (zc**l)/(l**2)
            end if
        end do
    psic = 2*et(i,j) - (gc2*(r(i,j)/gc1)**2)/pi - r(i,j)*(ux(i,j)*ux(i,j)+uy(i,j)*uy(i,j))
    if ((psia*psic) .lt. 0) then
        zb = zc
    else
        za = zc
    end if
  end do
  z(i,j) = zc
  t(i,j) = r(i,j)/(pi*gc1)
  p(i,j) = et(i,j) - 0.5*r(i,j)*(ux(i,j)**2+uy(i,j)**2)
end if
end do
end do 
go to 1300
!
1100 continue
do i = 1, nx
do j = 1, ny
if ((x(i).le.4.5) .or. (x(i).ge.5.5) .or. (y(j).le.4.5) .or. (y(j).ge.5.5)) then
  t(i,j) = (2* et(i,j)/r(i,j))-(ux(i,j)**2+uy(i,j)**2)
  z(i,j) = r(i,j)/(pi*t(i,j))
  p(i,j) = r(i,j)*t(i,j)/2
end if
end do
end do    
!
1300 continue
!@@@@@@@@@@@@@!
! write block !
!@@@@@@@@@@@@@!
if 	(iter.eq.floor(sv_point1/dt)) then
write(20,888)
write(20,*) 'zone t="left",i=',nx,',j=',floor((ny/(10/4.5))+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).le.4.5 ) then
   write (20,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(20,*) 'zone t="right",i=',nx,',j=',floor(ny/(10/4.5)+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).ge.5.5 ) then
   write (20,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(20,*) 'zone t="top",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).ge.5.5)) then
   write (20,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(20,*) 'zone t="bottom",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).le.4.5)) then
   write (20,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
end if
!
if 	(iter.eq.floor(sv_point2/dt)) then
write(25,888)
write(25,*) 'zone t="left",i=',nx,',j=',floor((ny/(10/4.5))+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).le.4.5 ) then
   write (25,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(25,*) 'zone t="right",i=',nx,',j=',floor(ny/(10/4.5)+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).ge.5.5 ) then
   write (25,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(25,*) 'zone t="top",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).ge.5.5)) then
   write (25,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(25,*) 'zone t="bottom",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).le.4.5)) then
   write (25,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
end if
!
if 	(iter.eq.floor(sv_point2b/dt)) then
write(27,888)
write(27,*) 'zone t="left",i=',nx,',j=',floor((ny/(10/4.5))+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).le.4.5 ) then
   write (27,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(27,*) 'zone t="right",i=',nx,',j=',floor(ny/(10/4.5)+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).ge.5.5 ) then
   write (27,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(27,*) 'zone t="top",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).ge.5.5)) then
   write (27,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(27,*) 'zone t="bottom",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).le.4.5)) then
   write (27,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
end if	  
!
if 	(iter.eq.floor(sv_point3/dt)) then
write(30,888)
write(30,*) 'zone t="left",i=',nx,',j=',floor((ny/(10/4.5))+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).le.4.5 ) then
   write (30,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(30,*) 'zone t="right",i=',nx,',j=',floor(ny/(10/4.5)+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).ge.5.5 ) then
   write (30,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(30,*) 'zone t="top",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).ge.5.5)) then
   write (30,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(30,*) 'zone t="bottom",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).le.4.5)) then
   write (30,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
end if
!
if 	(iter.eq.floor(sv_point3b/dt)) then
write(33,888)
write(33,*) 'zone t="left",i=',nx,',j=',floor((ny/(10/4.5))+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).le.4.5 ) then
   write (33,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(33,*) 'zone t="right",i=',nx,',j=',floor(ny/(10/4.5)+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).ge.5.5 ) then
   write (33,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(33,*) 'zone t="top",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).ge.5.5)) then
   write (33,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(33,*) 'zone t="bottom",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).le.4.5)) then
   write (33,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
end if

!
if 	(iter.eq.floor(sv_point4/dt)) then
write(35,888)
write(35,*) 'zone t="left",i=',nx,',j=',floor((ny/(10/4.5))+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).le.4.5 ) then
   write (35,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(35,*) 'zone t="right",i=',nx,',j=',floor(ny/(10/4.5)+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).ge.5.5 ) then
   write (35,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(35,*) 'zone t="top",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).ge.5.5)) then
   write (35,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(35,*) 'zone t="bottom",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).le.4.5)) then
   write (35,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
end if
!
if 	(iter.eq.floor(sv_point4b/dt)) then
write(37,888)
write(37,*) 'zone t="left",i=',nx,',j=',floor((ny/(10/4.5))+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).le.4.5 ) then
   write (37,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(37,*) 'zone t="right",i=',nx,',j=',floor(ny/(10/4.5)+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).ge.5.5 ) then
   write (37,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(37,*) 'zone t="top",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).ge.5.5)) then
   write (37,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(37,*) 'zone t="bottom",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).le.4.5)) then
   write (37,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
end if
!
write(*,777) time, r(nx/3,ny/2)
777	format (1X,'Elapsed time:',F7.4,4X, 'Density at x=4.0,y=5.:',F7.4)
if (istop .eq. 1) goto 2000
iter = iter + 1
goto 1000
pause
!
2000 continue
!
write(40,888)
write(40,*) 'zone t="left",i=',nx,',j=',floor((ny/(10/4.5))+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).le.4.5 ) then
   write (40,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(40,*) 'zone t="right",i=',nx,',j=',floor(ny/(10/4.5)+1)
do i = 1, nx
do j = 1, ny
   if ( x(i).ge.5.5 ) then
   write (40,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(40,*) 'zone t="top",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).ge.5.5)) then
   write (40,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
write(40,*) 'zone t="bottom",i=',floor(nx/(10/4.5)+1),',j=',floor(ny/10.+2)
do i = 1, nx
do j = 1, ny
   if (( x(i).ge.(4.5-dx) ) .and. ( x(i).le.(5.5+dx) ) .and. (y(j).le.4.5)) then
   write (40,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
end do 
end do
stop

end program    
