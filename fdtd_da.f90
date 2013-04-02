! ========================================================
subroutine em_da(q,s,aux,pml,nq,num_aux,num_pml,dx,dy,dt,nx,ny,xi,xf,yi,yf,nxi_pml,nxf_pml,nyi_pml,nyf_pml)
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
!	q = (Ex, Ey, Hz), s = (Dx, Dy, B3x, B3y), aux = (eps1, eps2, mu3)
!
!	if TE mode
!	q = (Hx, Hy, Ez), s = (B1x, B2y, E3x, E3y), aux = (mu1, mu2, eps3)

	implicit none 
	integer, intent(in) :: xi, xf, yi, yf, nq, num_aux, nx, ny
	integer, intent(in) :: num_pml, nxi_pml, nxf_pml, nyi_pml, nyf_pml
	double precision, intent(in) :: dx, dy, dt
	double precision, dimension(num_aux,ny,nx), intent(in) :: aux
	double precision, dimension(nq,ny,nx), intent(inout) :: q
	double precision, dimension(nq+1,ny,nx), intent(inout) :: s
	double precision, dimension(num_pml,ny,nx), intent(in) :: pml
	integer :: i,j,xfj,yfi
	
	if ((xi.gt.nxi_pml).and.(xf.lt.nx-nxf_pml+1).and.(yi.gt.nyi_pml).and.(yf.lt.ny-nyf_pml+1)) then
		do j = xi,xf-1
			do i = yi,yf-1
				s(3,i,j) = s(3,i,j) + (dt/dx)*(q(1,i,j) - q(1,i+1,j)) 
				s(4,i,j) = s(4,i,j) + (dt/dy)*(q(2,i,j+1) - q(2,i,j))
				q(3,i,j) = (s(3,i,j) + s(4,i,j))/aux(3,i,j)
			enddo
		enddo
		

		if (xf.eq.nx) then
			yfi = ny-1
		else
			yfi = yf
		endif

		do j = xi,xf-1
			do i = yi+1,yfi
				s(1,i,j) = s(1,i,j) + (dt/dy)*(q(3,i,j)-q(3,i-1,j))
				q(1,i,j) = s(1,i,j)/aux(1,i,j)
			enddo
		enddo

		if (xf.eq.nx) then
			xfj = nx-1
		else
			xfj = xf
		endif
		
		do j = xi+1,xfj
			do i = yi,yf-1
				s(2,i,j) = s(2,i,j) + (dt/dx)*(q(3,i,j-1)-q(3,i,j))
				q(2,i,j) = s(2,i,j)/aux(2,i,j)
			enddo
		enddo
	
	else
		
		do j = xi,xf-1
			do i = yi,yf-1
				s(3,i,j) = pml(1,i,j)*s(3,i,j) + pml(2,i,j)*(q(1,i,j) - q(1,i+1,j))
				s(4,i,j) = pml(3,i,j)*s(4,i,j) + pml(4,i,j)*(q(2,i,j+1) - q(2,i,j))
				q(3,i,j) = (s(3,i,j)+s(4,i,j))/aux(3,i,j)
			enddo
		enddo


		if (xf.eq.nx) then
			xfj = nx-1
		else
			xfj = xf
		endif

		do j = xi+1,xfj
			do i = yi,yf-1
				s(2,i,j) = pml(5,i,j)*s(2,i,j) + pml(6,i,j)*(q(3,i,j-1) - q(3,i,j))
				q(2,i,j) = s(2,i,j)/aux(2,i,j)
			enddo
		enddo


		if (yf.eq.ny) then
			yfi = ny-1
		else
			yfi = yf
		endif
		do j = xi,xf-1
			do i = yi+1,yfi
				s(1,i,j) = pml(7,i,j)*s(1,i,j) + pml(8,i,j)*(q(3,i,j) - q(3,i-1,j))
				q(1,i,j) = s(1,i,j)/aux(1,i,j)
			enddo
		enddo

	end if
	return
end subroutine em_da