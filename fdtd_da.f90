! ========================================================
subroutine em_da(q,s,aux,pml,nq,num_aux,num_pml,dx,dy,dt,nx,ny,&
	             xi,xf,yi,yf,nxi_pml,nxf_pml,nyi_pml,nyf_pml)
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

	implicit none 
	integer, intent(in) :: xi, xf, yi, yf, nq, num_aux, nx, ny
	integer, intent(in) :: num_pml, nxi_pml, nxf_pml, nyi_pml, nyf_pml
	double precision, intent(in) :: dx, dy, dt
	double precision, dimension(ny,nx,num_aux), intent(in) :: aux
	double precision, dimension(ny,nx,nq), intent(inout) :: q
	double precision, dimension(ny,nx,nq+1), intent(inout) :: s
	double precision, dimension(ny,nx,num_pml), intent(in) :: pml

	if((xi.gt.nxi_npml).and.(xf.lt.nx-nxf_npml).and.(yi.gt.nyi_npml).and.(yf.lt.ny-nyf_npml) then
		do j = xi+1,xf-1
			do i = yi+1,yf-1
				s(i,j,3) = s(i,j,3) - (dt/dx)*(q(i,j+1,2) - q(i,j,2)) 
				s(i,j,4) = s(i,j,4) + (dt/dy)*(q(i+1,j,1) - q(i,j,1))
				q(i,j,3) = (s(i,j,3) + s(i,j,4))/aux(i,j,3)
			enddo
		enddo
		do j = xi+1,xf-1
			do i = yi+1,yf-1
				s(i,j,1) = s(i,j,1) + (dt/dy)*(q(i,j,3)-q(i-1,j,3))
				q(i,j,1) = s(i,j,1)/aux(i,j,1)
			enddo
		enddo
		do j = xi+1,xf-1
			do i = yi+1,yf-1
				s(i,j,2) = s(i,j,2) - (dt/dx)*(q(i,j,3)-q(i,j-1,3))
				q(i,j,2) = s(i,j,2)/aux(i,j,2)
			enddo
		enddo
	else
		do j = xi+1,xf-1
			do i = yi+1,yf-1
				s(i,j,3) = pml(i,j,1)*s(i,j,3) - pml(i,j,2)*(q(i,j+1,2) - q(i,j,2)) 
				s(i,j,4) = pml(i,j,3)*s(i,j,4) + pml(i,j,4)*(q(i+1,j,1) - q(i,j,1))
				q(i,j,3) = (s(i,j,3) + s(i,j,4))/aux(i,j,3)
			enddo
		enddo
		do j = xi+1,xf-1
			do i = yi+1,yf-1
				s(i,j,1) = pml(i,j,5)*s(i,j,1) + pml(i,j,6)*(q(i,j,3)-q(i-1,j,3))
				q(i,j,1) = s(i,j,1)/aux(i,j,1)
			enddo
		enddo
		do j = xi+1,xf-1
			do i = yi+1,yf-1
				s(i,j,2) = pml(i,j,7)*s(i,j,2) - pml(i,j,8)*(q(i,j,3)-q(i,j-1,3))
				q(i,j,2) = s(i,j,2)/aux(i,j,2)
			enddo
		enddo
	endif
	return
end subroutine em_da