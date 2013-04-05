subroutine em_da_q3(aux,pml,s3,s4,q1,q2,q3,na,np,xi,xf,yi,yf,gxi,gxf,gyi,gyf)
  ! f2py intent(inplace) s3, s4
  implicit none
  integer, intent(in) :: na, np
  integer, intent(in) :: xi, xf, yi, yf, gxi, gxf, gyi, gyf
  double precision, dimension(na,gxi:gxf,gyi:gyf), intent(in)    :: aux
  double precision, dimension(np,gxi:gxf,gyi:gyf), intent(in)    :: pml
  double precision, dimension(   gxi:gxf,gyi:gyf), intent(inout) :: s3,s4
  double precision, dimension(   gxi:gxf,gyi:gyf), intent(in)    :: q1,q2
  double precision, dimension(   gxi:gxf,gyi:gyf), intent(inout) :: q3
  integer :: i,j, is, js, ie, je

  is = xi; ie = xf
  js = yi; je = yf

  if (gxf == xf) ie = xf-1
  if (gyf == yf) je = yf-1
  
  do j = js,je
     do i = is,ie
        s3(i,j) = pml(1,i,j)*s3(i,j) + pml(2,i,j)*(q1(i,j+1) - q1(i,j))
        s4(i,j) = pml(3,i,j)*s4(i,j) + pml(4,i,j)*(q2(i,j) - q2(i+1,j))
        q3(i,j) = (s3(i,j)+s4(i,j))/aux(3,i,j)
     enddo
  enddo

end subroutine em_da_q3

subroutine em_da_q12(aux,pml,s1,s2,q1,q2,q3,na,np,xi,xf,yi,yf,gxi,gxf,gyi,gyf)

  implicit none
  integer, intent(in) :: na, np
  integer, intent(in) :: xi, xf, yi, yf, gxi, gxf, gyi, gyf
  double precision, dimension(na,gxi:gxf,gyi:gyf), intent(in)    :: aux
  double precision, dimension(np,gxi:gxf,gyi:gyf), intent(in)    :: pml
  double precision, dimension(   gxi:gxf,gyi:gyf), intent(inout) :: s1,s2
  double precision, dimension(   gxi:gxf,gyi:gyf), intent(inout) :: q1,q2
  double precision, dimension(   gxi:gxf,gyi:gyf), intent(in)    :: q3
  integer :: i,j, is, js, ie, je

  is = xi; ie = xf
  js = yi; je = yf
  if (gxf == xf)  ie = xf-1
  if (gyi ==  1)  js = 2
  if (gyf == yf)  je = yf-1

  do j = js,je
     do i = is,ie
        s1(i,j) = pml(7,i,j)*s1(i,j) + pml(8,i,j)*(q3(i,j) - q3(i,j-1))
        q1(i,j) = s1(i,j)/aux(1,i,j)
     enddo
  enddo

  is = xi; ie = xf
  js = yi; je = yf
  if (gxi ==  1)  is = 2
  if (gxf == xf)  ie = xf-1
  if (gyf == yf)  je = yf-1 

  do j = js,je
     do i = is,ie
        s2(i,j) = pml(5,i,j)*s2(i,j) + pml(6,i,j)*(q3(i-1,j) - q3(i,j))
        q2(i,j) = s2(i,j)/aux(2,i,j)
     enddo
  enddo

end subroutine em_da_q12
