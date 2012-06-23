% TVD

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