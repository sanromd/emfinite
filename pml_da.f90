! ========================================================
subroutine pml_da(q,s,aux,pml,dx,dy,dt,xnp,ynp,nx,ny,nq,num_aux,&
	num_pml)
! ========================================================
!
!	Calculates the differential Maxwell equations for the inside
! 	grid points. PML and source are calculated outside this loop.
!
!	dx, dy and dt are the grid and time steps, respectively.
!	xi, xf and yi,yf are the intial and end points for the inner
!	grid points; whereas nx and ny are the total size of the grid.
!	
!	if TM mode
!	q = (E1, E2, H3), s = (D1, D2, B3x, B3y), aux = (eps1, eps2, mu3)
!
!	if TE mode
!	q = (H1, H2, E3), s = (B1, B2, E3x, E3y), aux = (mu1, mu2, eps3)

	integer, intent(in) :: xi, xf, yi, yf, nq, num_aux, nx, ny
	double precision, intent(in) :: dx, dy, dt
	double precision, dimension(nx,ny,num_aux), intent(in) :: aux
	double precision, dimension(ny,nx,nq), intent(inout) :: q
	double precision, dimension(ny,nx,nq+1), intent(inout) :: s
	double precision :: qz

	do j = yi+1,yf-1
		do i = xi+1,xf-1
			s(j,i,1) = s(j,i,1) + (dt/dy)*(q(j,i,3)-q(j-1,i,3))
			s(j,i,2) = s(j,i,2) - (dt/dx)*(q(j,i,3)-q(j,i-1,3))
			q(j,i,1) = s(j,i,1)/aux(j,i,1)
			q(j,i,2) = s(j,i,2)/aux(j,i,2)
			s(j,i,3) = s(j,i,3) - (dt/dx)*(q(j,i,2) - q(j,i-1,2)) 
			s(j,i,4) = s(j,i,4) + (dt/dy)*(q(j,i,1) - q(j-1,i,1))
			qz = s(j,i,3) + s(j,i,4)
			q(j,i,3) = qz/aux(j,i,3)
		enddo
	enddo
	return
end subroutine em_da