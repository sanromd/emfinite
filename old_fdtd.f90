	integer, parameter :: np, npml, norder, ns
	integer :: nmax, istep, jstep, ie, je, ib, jb, ip, jp, is, je
	double precision, parameter :: xx, yy, ddx, ddy, tp, top, alambda, xnbg, deltan, xoffe, yoffe, xoffu, yoffu, pvelx, pvely
	double precision, parameter :: rpex, rpey, rpux, rpuy
	double precision :: eguid, pi, epso, c o, dt, Ro, k0, k0e, rs0, rg0, rx0, muo, dt2, muo, vwave
	double precision :: vxe, vye, vxy, vyu, w, ms
	double precision, dimension(ib,jb,np) :: Px, Pxn, Pxn1, Py, Pyn, Pyn1, Mz, Mzn, Mzn1
	double precision, dimension(ib,jb) :: ex, ey, exn, eyn, hz, hzn, Dx, Dy, Bz, Bzx, Bzy, Sx, Sy
	double precision, dimension(ib,jb) :: caDx,cbDx, caDy, cbDy, daBzx, daBzy, dbBzx, dbBzy
	double precision, dimension(ib,jb) :: epsx, epsy, tfac e, tfac u, uoz, see
	double precision, dimension(np) :: C1, C2, C3, c1, c2, c3
	double precision, dimension(ib) :: x, vxe, vxu, fmode
	double precision, dimension(jb) :: y, vye, vyu, Pfi
	double precision, dimension(npml) :: sigex, sigmx, sigey, sigmy, cax, cbx, cay, cby, dax, dbx, day, dby
	double precision, dimension(np) :: eoo
	double precision, dimension(6) :: fi, gi, wi	

!	basi! parameters for the grid
	np = 3, xx = 1e-6, yy = 500e-6, ddx = 5e-9, ddy = 5e-9
	istep = 1, jstep = 2, ns = 60
	npml = 8, norder = 3, Ro = 1.0e-6
	nmax = 250000

!	wave propagation parameters
	tp = 15.0e-15, top = 3.2*tp, alambda = 1.0e-6

!	refra!tive index perturbation parameters
	xnbg = 1.5, deltan = 0.1*xnbg 
	xoffe = xx/2, yoffe = 458e-6, xoffu = xoffe, yoffu = yoffe
	rpex = .5*alambda, rpey = 5.0*alambda
	rpux = .5*alambda, rpuy = 5.0*alambda
	pvelx = 0.0, pvely = -0.61

!	pre-!al!ulations
	
	pi = 4.0*atan(1.0), epso = 8.854e-12
	co = 1.0/sqrt(4.0*pi*1.0e-7*epso)
	dt = 0.90/(co*sqrt(1.0/(ddx**2)  +  1.0/(ddy**2)))
	ie = floor(xx/ddx), je = floor(yy/ddy)
	ib = ie  + 1, jb = he  +  1
	ip = ie - npml, jp = je - npml
	is = floor(ie/2), js = floor(je/2)
	ms = nmax/ns
	eguid = xnbg**2
	dt2 = dt*dt
	muo = 4.0*pi*1.0e-7
	vwave = co/xnbg
	vxe = pvelx*co
	vye = pvely*co
	vxu = vxe
	vyu = vye	
	w = 2.0*pi*co/alambda          

! Free spa!e/diele!tri! (1)
	eoo(1) = 1.0
	c1(1) = 0.0
	c2(1) = 0.0
	c3(1) = 0.0

! Lorentz (single pole, Diele!tri!) Pole 2-e
	eoo(2) = xnbg**2
	wo2 = 0.7d16**2
	es = xnbg**2
	tampe = deltan*(2.0*xnbg  +  deltan)
	dlta = 0.0*8e11
	ag = epso*(es - eoo(2))*wo2
	bg = wo2
	cg = 2.0*dlta
	dg = 1.0
	c1(2) = 0.0*(4.0*dg - 2.0*bg*dt2)/(2.0*dg  +  cg*dt)
	c2(2) = 0.0*(-2.0*dg  +  cg*dt)/(2.0*dg  +  cg*dt)
	c3(2) = 0.0*(2.0*ag*dt2)/(2.0*dg  +  cg*dt)

! Lorentz (single pole, Diele!tri!) Pole 3-u
	eoo(3) = xnbg**2
	wo2 = 0.7d16**2
	es = xnbg**2
	tampu = 0.0
	dlta = 0.0*8e11
	ag = muo*(es - eoo(3))*wo2
	bg = wo2
	cg = 2.0*dlta
	dg = 1.0
	c1(3) = 0.0*(4.0*dg - 2.0*bg*dt2)/(2.0*dg  +  cg*dt)
	c2(3) = 0.0*(-2.0*dg  +  cg*dt)/(2.0*dg  +  cg*dt)
	c3(3) = 0.0*(2.0*ag*dt2)/(2.0*dg  +  cg*dt)

!	Data Files
	open(10,file = 'profile')
	open(11,file = 'parameters')
	open(12,file = 'tdata')
	open(13,file = 'infodata')
	open(14,file = 'tfield1')
	open(24,file = 'tfield2')
	open(34,file = 'tfield3')
	open(44,file = 'tfield4')
	open(15,file = 'parameter')
	open(16,file = 'mode')
	write(11,*)ib,jb
	write(13,*)ddx/1e-6,ddy/1e-6,dt/1e-15,nmax,ms

! Visualization

	write(*,*)ib,jb,nmax,dt/1e-15,tampe,tampu
	write(*,*)xnbg,deltan

!	Material Definition

	do i = 1,ib
		x(i) = (i - 1)*ddx
	enddo
    do j = 1,jb
		y(j) =  (j - 1)*ddy
	enddo

!	Guided Mode		
	do i = 1,ie
		fmode0(i) = cos(pi*x(i)/xx)
		write(16,*)x(i),fmode0(i)
	enddo
!
	write(*,*)'stop = 1, !ontinue = 0'
	read(*,*)istop

	if(istop.eq.0) then

		epsx  = eguid*epso
		epsy  = eguid*epso
		uoz  = muo

!	PML cal!ulations

		sigmex = -(norder  +  1.0)*co*log(Ro)/(2.0*ddx*npml)
		sigmmx = sigmex
		sigmey = -(norder  +  1.0)*co*log(Ro)/(2.0*ddy*npml)
		sigmmy = sigmey

		do m = 1,npml
			sigex(m) = sigmex*((m - 0.5)/(npml  +  0.5))**norder
			sigmx(m) = sigmmx*(m/(npml  +  0.5))**norder
			sigey(m) = sigmey*((m - 0.5)/(npml + 0.5))**norder
			sigmy(m) = sigmmy*(m/(npml  +  0.5))**norder
		enddo

		do m = 1,npml
			rex = sigex(m)*dt
			rmx = sigmx(m)*dt
			cax(m) = exp(-rex)
			cbx(m) = -(exp(-rex)-1.0)/sigex(m)/ddx
			dax(m) = exp(-rmx)
			dbx(m) = -(exp(-rmx)-1.0)/sigmx(m)/ddx
			rey = sigey(m)*dt
			rmy = sigmy(m)*dt
			cay(m) = exp(-rey)
			cby(m) = -(exp(-rey)-1.0)/sigey(m)/ddy
			day(m) = exp(-rmy)
			dby(m) = -(exp(-rmy)-1.0)/sigmy(m)/ddy
		enddo

		Pxn1 = 0.0
		Pxn = 0.0
		Px = 0.0
		ex = 0.0
		exn = 0.0
		Dx = 0.0
		caDx = 1.0
		cbDx = dt/ddy
		Sx = 0.0
		Pyn1 = 0.0
		Pyn = 0.0
		Py = 0.0
		ey = 0.0
		eyn = 0.0
		Dy = 0.0
 		caDy = 1.0
		cbDy = dt/ddx
		Sy = 0.0
		hz = 0.0
		Bz = 0.0
		Bzx = 0.0
		Bzy = 0.0
		daBzx = 1.0
		dbBzx = dt/ddx
		daBzy = 1.0
		dbBzy = dt/ddy

		Pfi = 0.0

! 	PML
! 	Left (Dx)
		do 100 i = 1,ie
!			do j = 2,npml + 1
!				m = npml + 2-j
!				caDx(i,j) = cay(m)
!				cbDx(i,j) = cby(m)
!			enddo	
			do j = jp + 1,je
				m = j-jp
				caDx(i,j) = cay(m)
				cbDx(i,j) = cby(m)
			enddo
		enddo

! 	Bottom(Dy)
!		do j = 1,je
!			do i = 2,npml + 1
!				m = npml + 2-i
!				caDy(i,j) = cax(m)
!				cbDy(i,j) = cbx(m)
!			enddo
! 	Top(Dy)
!			do i = ip + 1,ie
!				m = i-ip
!				caDy(i,j) = cax(m)
!				cbDy(i,j) = cbx(m)
!			enddo
!		enddo

! 	Left(Bz)
		do i = 1,ie
!			do 610 j = 1,npml
!				m = npml + 1-j
!				daBzy(i,j) = day(m)
!				dbBzy(i,j) = dby(m)
!			enddo
! 	Right(Bz)
			do j = jp + 1,je
				m = j-jp
				daBzy(i,j) = day(m)
				dbBzy(i,j) = dby(m)
			enddo
		enddo

! 	Bottom(Bz)
!		do j = 1,je
!			do i = 1,npml
!				m = npml + 1-i
!				daBzx(i,j) = dax(m)
!				dbBzx(i,j) = dbx(m)
!			enddo
! 	Top(Bz)
!			do i = ip + 1,ie
!				m = i-ip
!				daBzx(i,j) = dax(m)
!				dbBzx(i,j) = dbx(m)
!			enddo
!		enddo

		write(6,*)'Initialization !omplete'
!
! 	-------------- Time-stepping loop
!
		do n = 1,nmax

			if(mod(n,100).eq.0) then
				write(*,*)n
			endif
			
			time = dt*float(n)

	        do i = 1,ie
	        	do j = 1,je
					tface(i,j) = 1
					tfacu(i,j) = 1
	    			tface(i,j) = epso*(tampe*exp(-(0.0*((x(i) - xoffe) - vxe(i)*time)**2/rpex**2 + ((y(j) - yoffe-rpey/2.0) -&
								 vye(j)*time)**2/rpey**2))) + epso*(-1.0*tampe*exp(-(0.0*((x(i) - xoffe) -&
								 vxe(i)*time)**2/rpex**2 + ((y(j)-yoffe + rpey/2.0) - vye(j)*time)**2/rpey**2)))
	    			tfacu(i,j) = muo*(tampu*exp(-(0.0*((x(i) - xoffu) - vxu(i)*time)**2/rpux**2 + ((y(j) - yoffu - rpuy/2.0) -&
								 vyu(j)*time)**2/rpuy**2))) + muo*(-1.0*tampu*exp(-(0.0*((x(i) - xoffu) -&
								 vxu(i)*time)**2/rpux**2 + ((y(j) - yoffu + rpuy/2.0) - (j)*time)**2/rpuy**2)))
				enddo
			enddo

!			tpulse = sin(w*time)
   			tpulse = (1.0 - exp(-(time/tp)**2))*cos(w*time)
!  			tpulse = exp(-((time - to)/tp)**2)*cos(w*time)

!   -------------- (Bz, Hz)

			Bzx(1:ie,1:je) = daBzx(1:ie,1:je)*Bzx(1:ie,1:je) + dbBzx(1:ie,1:je)*(ey(1:ie,1:je) - ey(2:ib,1:je))
			Bzy(1:ie,1:je) = daBzy(1:ie,1:je)*Bzy(1:ie,1:je) + dbBzy(1:ie,1:je)*(ex(1:ie,2:jb) - ex(1:ie,1:je))
			Bz(1:ie,1:je) = Bzx(1:ie,1:je) + Bzy(1:ie,1:je)
			hz(1:ie,1:je) = Bz(1:ie,1:je)/uoz(1:ie,1:je)
!			hz(1:ie,1) = fmode0(1:ie)*tpulse

!   ------------- (Dx, Px, Ex)

			Dx(1:ie,2:je) = caDx(1:ie,2:je)*Dx(1:ie,2:je) + cbDx(1:ie,2:je)*(hz(1:ie,2:je) - hz(1:ie,1:je-1))
			ex(1:ie,2:je) = Dx(1:ie,2:je)/(epsx(1:ie,2:je) + tface(1:ie,2:je))
			ex(1:ie,1) = tpulse*fmode0(1:ie)

!	-------------- (Dy, Py, Ey)

			Dy(2:ie,1:je) = caDy(2:ie,1:je)*Dy(2:ie,1:je) + cbDy(2:ie,1:je)*(hz(1:ie-1,1:je) - hz(2:ie,1:je))
			ey(2:ie,1:je) = Dy(2:ie,1:je)/(epsy(2:ie,1:je) + tface(2:ie,1:je))

!	-------------- Sx, Sy
		
			Sx(2:ie,2:je) = Hz(2:ie,2:je)*(Ey(3:ib,2:je) + Ey(2:ie,2:je))/2.0
			Sy(2:ie,2:je) = -1.*Hz(2:ie,2:je)*(Ex(2:ie,2:je) + Ex(2:ie,3:jb))/2.

!	Movie 
	 		if(mod(n,ms).eq.0) then
				do i = 1,ie,istep
					do j = 1,je,jstep
	   					write(14,*)ey(i,j),ex(i,j)
	   					write(24,*)tface(i,j),hz(i,j)
						write(34,*)Sx(i,j),Sy(i,j)
					enddo
				enddo
			endif
			write(14,94)n,ey(100, 10000)
			write(24,94)n,ey(100, 35000)
			write(34,94)n,ey(100, 60000)
			write(44,94)n,ey(100, 85000)
94			format(i7,2x,1d24.13)
		enddo

	write(15,1)i3,j3,ns,iw1,iw2,jw1,jw2,is,js,istep,jstep
	format(11(i5,x))

	endif
