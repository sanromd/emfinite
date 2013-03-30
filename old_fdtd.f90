
c	Metallic Waveguides    TM
c   1D Gaussian - collision
c
c       2D FDTD implementing the General Algorithm
c       and MIPML
c
c						    x(i)|	 
c	Thickness of PML layer=npml 			|  
c	Total number of time steps=nmax			|_____ y(j)
c	Domain size= ib x jb	   		   	
c	with movie							 
c									     
c	Simulation solves the RIP problem assuming mu=cte
c
c	M. A. Alsunaidi
c	KFUPM
c	August 9, 2011
c
c	------------------------------------------
c	cosine wave - no trap - metallic  
c	------------------------------------------
c
c	Set basic intialization of definitions
	parameter (np=3)! specify poles
	parameter (xx=1e-6, yy=500e-6)
	parameter (ddx=5e-9, ddy=5e-9)
	parameter (tp=15.0e-15, to=3.2*tp, alambda=1.0e-6)
c 	Time independent material definition, e.g. background waveguide   	
	parameter (xnbg=1.5)
	parameter (deltan=0.1*xnbg)
	parameter (eguid=(xnbg)**2)!dielectric 
c	Time-dependent refractive index definition
  	parameter (rpex=.5*alambda, rpey=5.0*alambda)!epsilon dispersion constant, i.e. sigma epsilon
   	parameter (rpux=.5*alambda, rpuy=5.0*alambda)!mu dispersion constant, i.e. sigma mu
	parameter (xoffe=xx/2, yoffe=458e-6, xoffu=xoffe, yoffu=yoffe)!x,y - offset
	parameter (pi=4.0*atan(1.0), epso=8.854e-12)
	parameter (co=1.0/sqrt(4.0*pi*1.0e-7*epso))
	parameter (dt=0.90/(co*sqrt(1.0/(ddx**2)+1.0/(ddy**2))))
	parameter (pvelx=0.0, pvely=-0.61, nmax=250000)
c Computational Domain
	parameter (npml=8, norder=3, Ro=1.0e-6)
   	parameter (ie=xx/ddx, je=yy/ddy)
	parameter (ib=ie+1, jb=je+1)
	parameter (ip=ie-npml, jp=je-npml)
	parameter (is=ie/2, js=je/2)
	parameter (istep=1, jstep=2, ns=60, ms=nmax/ns)
c
	dimension Px(ib,jb,np),Pxn(ib,jb,np),Pxn1(ib,jb,np)
	dimension Py(ib,jb,np),Pyn(ib,jb,np),Pyn1(ib,jb,np)
	real Mz(ib,jb,np),Mzn(ib,jb,np),Mzn1(ib,jb,np)
	dimension C1(np),C2(np),C3(np)
	dimension ex(ib,jb),ey(ib,jb)
	dimension exn(ib,jb),eyn(ib,jb),hzn(ib,jb)
	dimension Dx(ib,jb),Dy(ib,jb)
	dimension Sx(ib,jb),Sy(ib,jb)
	dimension caDx(ib,jb),cbDx(ib,jb)
	dimension caDy(ib,jb),cbDy(ib,jb)
	dimension hz(ib,jb),tface(ib,jb)
	dimension Bz(ib,jb),tfacu(ib,jb)
	dimension Bzx(ib,jb),Bzy(ib,jb)
	dimension daBzx(ib,jb),dbBzx(ib,jb)
	dimension daBzy(ib,jb),dbBzy(ib,jb)
	dimension epsx(ib,jb),epsy(ib,jb),uoz(ib,jb)
	dimension Pfi(jb),x(ib),y(jb),see(ib,jb)
	dimension vxe(ib),vxu(ib),vye(jb),vyu(jb)
	dimension sigex(npml),sigmx(npml)
	dimension sigey(npml),sigmy(npml)
	dimension cax(npml),cbx(npml)
	dimension cay(npml),cby(npml)
	dimension dax(npml),dbx(npml)
	dimension day(npml),dby(npml)
	dimension eoo(np),fi(6),gi(6),wi(6)
	real k0, k0e, rs0, rg0, rc0
	real fmode0(ib)
	real muo 	

c	Set basic parameters of simulation, dimensions, and step size
	dt2=dt*dt
	muo=4.0*pi*1.0e-7
	vwave=co/xnbg
	vxe=pvelx*co
	vye=pvely*co
	vxu=vxe
	vyu=vye

c 	Source defintion, rasing time (tp), and wavelenth  	
	w=2.0*pi*co/alambda          

c Free space/dielectric (1)
	eoo(1)=1.0
	c1(1)=0.0
	c2(1)=0.0
	c3(1)=0.0

c Lorentz (single pole, Dielectric) Pole 2-e
	eoo(2)=(xnbg)**2
	wo2=(0.7d16)**2
	es=(xnbg)**2
        tampe=deltan*(2.0*xnbg+deltan)
	dlta=0.0*8e11
	ag=epso*(es-eoo(2))*wo2
	bg=wo2
	cg=2.0*dlta
	dg=1.0
	c1(2)=0.0*(4.0*dg-2.0*bg*dt2)/(2.0*dg+cg*dt)
	c2(2)=0.0*(-2.0*dg+cg*dt)/(2.0*dg+cg*dt)
	c3(2)=0.0*(2.0*ag*dt2)/(2.0*dg+cg*dt)

c Lorentz (single pole, Dielectric) Pole 3-u
	eoo(3)=(xnbg)**2
	wo2=(0.7d16)**2
	es=(xnbg)**2
        tampu=0.0
	dlta=0.0*8e11
	ag=muo*(es-eoo(3))*wo2
	bg=wo2
	cg=2.0*dlta
	dg=1.0
	c1(3)=0.0*(4.0*dg-2.0*bg*dt2)/(2.0*dg+cg*dt)
	c2(3)=0.0*(-2.0*dg+cg*dt)/(2.0*dg+cg*dt)
	c3(3)=0.0*(2.0*ag*dt2)/(2.0*dg+cg*dt)
c
c	Data Files
	open(10,file='profile')
	open(11,file='parameters')
	open(12,file='tdata')
	open(13,file='infodata')
	open(14,file='tfield1')
	open(24,file='tfield2')
	open(34,file='tfield3')
	open(44,file='tfield4')
	open(15,file='parameter')
	open(16,file='mode')

        write(11,*)ib,jb
        write(13,*)ddx/1e-6,ddy/1e-6,dt/1e-15,nmax,ms
c
c Visualization

	write(*,*)ib,jb,nmax,dt/1e-15,tampe,tampu
	write(*,*)xnbg,deltan
c
c	Material Definition

        do 22 i=1,ib
	x(i)=(i-1)*ddx
22	continue
        do 21 j=1,jb
	y(j)= (j-1)*ddy
21	continue


c	Guided Mode
		
	do 101 i=1,ie
	fmode0(i)=cos(pi*x(i)/xx)
	write(16,*)x(i),fmode0(i)
101	continue
c c
	write(*,*)'stop=1, continue=0'
	read(*,*)istop
	if(istop.eq.1) goto 999

	epsx=eguid*epso
	epsy=eguid*epso
	uoz=muo

c	PML Calculations

	sigmex=-(norder+1.0)*co*log(Ro)/(2.0*ddx*npml)
	sigmmx=sigmex
	sigmey=-(norder+1.0)*co*log(Ro)/(2.0*ddy*npml)
	sigmmy=sigmey
c
	do 10 m=1,npml
	sigex(m)=sigmex*((m-0.5)/(npml+0.5))**norder
	sigmx(m)=sigmmx*(m/(npml+0.5))**norder
	sigey(m)=sigmey*((m-0.5)/(npml+0.5))**norder
	sigmy(m)=sigmmy*(m/(npml+0.5))**norder
10	enddo
c
	do 20 m=1,npml
	rex=sigex(m)*dt
	rmx=sigmx(m)*dt
	cax(m)=exp(-rex)
	cbx(m)=-(exp(-rex)-1.0)/sigex(m)/ddx
	dax(m)=exp(-rmx)
	dbx(m)=-(exp(-rmx)-1.0)/sigmx(m)/ddx
	rey=sigey(m)*dt
	rmy=sigmy(m)*dt
	cay(m)=exp(-rey)
	cby(m)=-(exp(-rey)-1.0)/sigey(m)/ddy
	day(m)=exp(-rmy)
	dby(m)=-(exp(-rmy)-1.0)/sigmy(m)/ddy
20	enddo
c
	Pxn1=0.0
	Pxn=0.0
	Px=0.0
	ex=0.0
	exn=0.0
	Dx=0.0
	caDx=1.0
	cbDx=dt/ddy
	Sx=0.0
c
	Pyn1=0.0
	Pyn=0.0
	Py=0.0
	ey=0.0
	eyn=0.0
	Dy=0.0
	caDy=1.0
	cbDy=dt/ddx
	Sy=0.0
c
	hz=0.0
	Bz=0.0
	Bzx=0.0
	Bzy=0.0
	daBzx=1.0
	dbBzx=dt/ddx
	daBzy=1.0
	dbBzy=dt/ddy

	Pfi=0.0
c
c 	PML
c 	Left (Dx)
	do 100 i=1,ie
c	do 110 j=2,npml+1
c	m=npml+2-j
c	caDx(i,j)=cay(m)
c	cbDx(i,j)=cby(m)
c 	110	continue
c. Right (Dx)
	do 120 j=jp+1,je
	m=j-jp
	caDx(i,j)=cay(m)
	cbDx(i,j)=cby(m)
120	continue
100	continue
c
c 	Bottom(Dy)
c	do 200 j=1,je
c	do 210 i=2,npml+1
c	m=npml+2-i
c	caDy(i,j)=cax(m)
c	cbDy(i,j)=cbx(m)
c 	210	continue
c 	Top(Dy)
c	do 220 i=ip+1,ie
c	m=i-ip
c	caDy(i,j)=cax(m)
c	cbDy(i,j)=cbx(m)
c 	220	continue
c 	200	continue
c
c
c 	Left(Bz)
	do 600 i=1,ie
c	do 610 j=1,npml
c	m=npml+1-j
c	daBzy(i,j)=day(m)
c	dbBzy(i,j)=dby(m)
c 	610	continue
c 	Right(Bz)
	do 620 j=jp+1,je
	m=j-jp
	daBzy(i,j)=day(m)
	dbBzy(i,j)=dby(m)
620	continue
600	continue
c 	Bottom(Bz)
c	do 630 j=1,je
c	do 640 i=1,npml
c	m=npml+1-i
c	daBzx(i,j)=dax(m)
c	dbBzx(i,j)=dbx(m)
c 	640	continue
c 	Top(Bz)
c	do 650 i=ip+1,ie
c	m=i-ip
c	daBzx(i,j)=dax(m)
c	dbBzx(i,j)=dbx(m)
c 	650	continue
c 	630	continue
c
	write(6,*)'Initialization Complete'
c
c 	-------------- Time-stepping loop
c
	do 99 n=1,nmax
	if(mod(n,100).eq.0)write(6,*)n
	time=dt*float(n)
c
        do 291 i=1,ie
        do 291 j=1,je
		tface(i,j)=1
		tfacu(i,j)=1
    	tface(i,j)=epso*(tampe*exp(-(0.0*((x(i)-xoffe)-
     .       vxe(i)*time)**2/rpex**2+((y(j)-yoffe-rpey/2.0)-
     .       vye(j)*time)**2/rpey**2)))+
     .	     epso*(-1.0*tampe*exp(-(0.0*((x(i)-xoffe)-
     .       vxe(i)*time)**2/rpex**2+((y(j)-yoffe+rpey/2.0)-
     .       vye(j)*time)**2/rpey**2)))
    	tfacu(i,j)=muo*(tampu*exp(-(0.0*((x(i)-xoffu)-
     .       vxu(i)*time)**2/rpux**2+((y(j)-yoffu-rpuy/2.0)-
     .       vyu(j)*time)**2/rpuy**2)))+
     .		 muo*(-1.0*tampu*exp(-(0.0*((x(i)-xoffu)-
     .       vxu(i)*time)**2/rpux**2+((y(j)-yoffu+rpuy/2.0)-
     .       vyu(j)*time)**2/rpuy**2)))
291	continue

c	tpulse=sin(w*time)
   	tpulse=(1.0-exp(-(time/tp)**2))*cos(w*time)
c  	tpulse=exp(-((time-to)/tp)**2)*cos(w*time)

c   -------------- (Bz, Hz)
c
	Bzx(1:ie,1:je)=daBzx(1:ie,1:je)*Bzx(1:ie,1:je)+
     .   dbBzx(1:ie,1:je)*(ey(1:ie,1:je)-ey(2:ib,1:je))
	Bzy(1:ie,1:je)=daBzy(1:ie,1:je)*Bzy(1:ie,1:je)+
     .   dbBzy(1:ie,1:je)*(ex(1:ie,2:jb)-ex(1:ie,1:je))
	Bz(1:ie,1:je)=Bzx(1:ie,1:je)+Bzy(1:ie,1:je)

	hz(1:ie,1:je)=Bz(1:ie,1:je)/uoz(1:ie,1:je)
c	hz(1:ie,1)=fmode0(1:ie)*tpulse
c
c   ------------- (Dx, Px, Ex)
c
	Dx(1:ie,2:je)=caDx(1:ie,2:je)*Dx(1:ie,2:je)+
     .   cbDx(1:ie,2:je)*(hz(1:ie,2:je)-hz(1:ie,1:je-1))

	ex(1:ie,2:je)=Dx(1:ie,2:je)/
     .         (epsx(1:ie,2:je)+tface(1:ie,2:je))

	ex(1:ie,1)=tpulse*fmode0(1:ie)
                              
cc
c      -------------- (Dy, Py, Ey)
c
	Dy(2:ie,1:je)=caDy(2:ie,1:je)*Dy(2:ie,1:je)+
     .   cbDy(2:ie,1:je)*(hz(1:ie-1,1:je)-hz(2:ie,1:je))

	ey(2:ie,1:je)=Dy(2:ie,1:je)/
     .         (epsy(2:ie,1:je)+tface(2:ie,1:je))

cc
c	Update Samples
	
c
	Sx(2:ie,2:je)=Hz(2:ie,2:je)*(Ey(3:ib,2:je)+Ey(2:ie,2:je))/2.0
	Sy(2:ie,2:je)=-1.*Hz(2:ie,2:je)*(Ex(2:ie,2:je)+Ex(2:ie,3:jb))/2.

c	Movie 
c 	
c 	if(mod(n,ms).eq.0)then
c	i3=0
c   do 777 i=1,ie,istep
c		j3=0	
c		i3=i3+1
c		do 778 j=1,je,jstep
c			j3=j3+1
c   		write(14,*)ey(i,j),ex(i,j)
c   		write(24,*)tface(i,j),hz(i,j)
c			write(34,*)Sx(i,j),Sy(i,j)
c 778	continue
c 777continue
c	endif
c
	write(14,94)n,ey(100, 10000)
	write(24,94)n,ey(100, 35000)
	write(34,94)n,ey(100, 60000)
	write(44,94)n,ey(100, 85000)
94	format(i7,2x,1d24.13)
c
99 	continue
	write(15,1)i3,j3,ns,iw1,iw2,jw1,jw2,is,js,istep,jstep
c
1	format(11(i5,x))
c	
999	end
