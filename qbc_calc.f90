! ========================================================
subroutine qbc_calc(q,s,aux,pml,xsrc,ysrc,nq,num_aux,num_pml,nx,ny,ixls,ixrs,iyus,iyls)
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
	integer, intent(in) :: nq, num_aux, nx, ny
	integer, intent(in) :: num_pml,ixls,ixrs,iyus,iyls
	double precision, dimension(1:ny-1,2), intent(in) :: xsrc
	double precision, dimension(1:nx-1,2), intent(in) :: ysrc
	double precision, dimension(num_aux,ny,nx), intent(in) :: aux
	double precision, dimension(nq,ny,nx), intent(inout) :: q
	double precision, dimension(nq+1,ny,nx), intent(inout) :: s
	double precision, dimension(num_pml,ny,nx), intent(in) :: pml
	integer :: i,j


! ----- left boundary (j=1)		
	do i = 1,ny-1
		s(3,i,1) = pml(1,i,1)*s(3,i,1) - pml(2,i,1)*(q(2,i,2) - q(2,i,1)) 
		s(4,i,1) = pml(3,i,1)*s(4,i,1) + pml(4,i,1)*(q(1,i+1,1) - q(1,i,1))
		q(3,i,1) = (s(3,i,1) + s(4,i,1))/aux(2,i,1)
	enddo

	do i = 2,ny-1
		s(1,i,1) = pml(5,i,1)*s(1,i,1) + pml(6,i,1)*(q(3,i,1)-q(3,i-1,1))
		q(1,i,1) = s(1,i,1)/aux(1,i,1)
	enddo

	if (ixls.eq.1) then
		q(2,1:ny-1,1) = xsrc(1:ny-1,1)
	else
		do i = 1,ny-1
			s(2,i,1) = pml(7,i,1)*s(2,i,1) - pml(8,i,1)*(q(3,i,1))
			q(2,i,1) = s(2,i,1)/aux(1,i,1)
		enddo
	end if

! ----- right boundary

	do i = 1,ny-1
		s(3,i,nx) = pml(1,i,nx)*s(3,i,nx) - pml(2,i,nx)*(-q(2,i,nx)) 
		s(4,i,nx) = pml(3,i,nx)*s(4,i,nx) + pml(4,i,nx)*(q(1,i+1,nx) - q(1,i,nx))
		q(3,i,nx) = (s(3,i,nx) + s(4,i,nx))/aux(2,i,nx)
	enddo
!
	do i = 2,ny-1
		s(1,i,nx) = pml(5,i,nx)*s(1,i,nx) + pml(6,i,nx)*(q(3,i,nx)-q(3,i-1,nx))
		q(1,i,nx) = s(1,i,nx)/aux(1,i,nx)
	enddo
!
!	if (ixrs.eq.1) then
!		q(2,1:ny-1,nx) = xsrc(1:ny-1,2)
!	else
		do i = 1,ny-1
			s(2,i,nx) = pml(7,i,nx)*s(2,i,nx) - pml(8,i,nx)*(q(3,i,nx)-q(3,i,nx-1))
			q(2,i,nx) = s(2,i,nx)/aux(1,i,nx)
		enddo
!	end if

! ----- upper boundary (i=1)

	do j = 2,nx-1
		s(3,1,j) = pml(1,1,j)*s(3,1,j) - pml(2,1,j)*(q(2,1,j+1) - q(2,1,j)) 
		s(4,1,j) = pml(3,1,j)*s(4,1,j) + pml(4,1,j)*(q(1,2,j) - q(1,1,j))
		q(3,1,j) = (s(3,1,j) + s(4,1,j))/aux(2,1,j)
	enddo

	if (iyus.eq.1) then
		q(1,1,1:nx-1) = ysrc(1:nx-1,1)
	else
		do j = 2 ,nx-1
			s(1,1,j) = pml(5,1,j)*s(1,1,j) + pml(6,1,j)*(q(3,1,j))
			q(1,1,j) = s(1,1,j)/aux(1,1,j)
		enddo
	end if

	do j = 2,nx-1
		s(2,1,j) = pml(7,1,j)*s(2,1,j) - pml(8,1,j)*(q(3,1,j)-q(3,1,j-1))
		q(2,1,j) = s(2,1,j)/aux(1,1,j)
	enddo

! ----- lower boundary

	do j = 2,nx-1
		s(3,ny,j) = pml(1,ny,j)*s(3,ny,j) - pml(2,ny,j)*(q(2,ny,j+1) - q(2,ny,j)) 
		s(4,ny,j) = pml(3,ny,j)*s(4,ny,j) + pml(4,ny,j)*(-q(1,ny,j))
		q(3,ny,j) = (s(3,ny,j) + s(4,ny,j))/aux(2,ny,j)
	enddo

!	if (iyls.eq.1) then
!		q(1,ny,1:nx-1) = ysrc(1:nx-1,2)
!	else
		do j = 2,nx-1
			s(1,ny,j) = pml(5,ny,j)*s(1,ny,j) + pml(6,ny,j)*(q(3,ny,j)-q(3,ny-1,j))
			q(1,ny,j) = s(1,ny,j)/aux(1,ny,j)
		enddo
!	end if

	do j = 2,nx-1
		s(2,ny,j) = pml(7,ny,j)*s(2,ny,j) - pml(8,ny,j)*(q(3,ny,j)-q(3,ny,j-1))
		q(2,ny,j) = s(2,ny,j)/aux(1,ny,j)
	enddo
end subroutine qbc_calc