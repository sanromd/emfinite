subroutine fdtd2D(aux,pml,dxdt,dydt,s1,s2,s3,s4,q1,q2,q3,na,np,xi,xf,yi,yf,gxi,gxf,gyi,gyf,q_type,pml_decomp)

  implicit none
  integer, intent(in) :: na, np, q_type, pml_decomp
  integer, intent(in) :: xi, xf, yi, yf, gxi, gxf, gyi, gyf
  double precision, dimension(na, xi:xf,  yi:yf),  intent(in)    :: aux
  double precision, dimension(np, xi:xf,  yi:yf),  intent(in)    :: pml
  double precision, dimension(    xi:xf,  yi:yf),  intent(inout) :: s1,s2,s3,s4
  double precision, dimension(   gxi:gxf,gyi:gyf), intent(inout) :: q1,q2,q3
  double precision, intent(in) :: dxdt,dydt
  integer :: i,j, is, js, ie, je

  select case(q_type) 
    case(0) !calculate q3
      is = xi; ie = xf
      js = yi; je = yf

      if (gxf == xf) ie = xf-1
      if (gyf == yf) je = yf-1
      select case(pml_decomp)
        case(1)
          do j = js,je
            do i = is,ie
              s3(i,j) = pml(1,i,j)*s3(i,j) + pml(2,i,j)*(q1(i,j+1) - q1(i,j))
              s4(i,j) = pml(3,i,j)*s4(i,j) + pml(4,i,j)*(q2(i,j) - q2(i+1,j))
              q3(i,j) = (s3(i,j)+s4(i,j))/aux(3,i,j)
            enddo
          enddo
        case default
          do j = js,je
            do i = is,ie
              s3(i,j) = s3(i,j) + dydt*(q1(i,j+1) - q1(i,j))
              s4(i,j) = s4(i,j) + dxdt*(q2(i,j) - q2(i+1,j))
              q3(i,j) = (s3(i,j)+s4(i,j))/aux(3,i,j)
            enddo
          enddo
        end select
    case(1) !calculate q1 and q2 
      ! q1
      is = xi; ie = xf
      js = yi; je = yf
      if (gxf == xf)  ie = xf-1
      if (gyi ==  1)  js = 2
      if (gyf == yf)  je = yf-1
      select case(pml_decomp)
        case(1)
          do j = js,je
            do i = is,ie
              s1(i,j) = pml(7,i,j)*s1(i,j) + pml(8,i,j)*(q3(i,j) - q3(i,j-1))
              q1(i,j) = s1(i,j)/aux(1,i,j)
            enddo
          enddo
        case default
          do j = js,je
            do i = is,ie
                s1(i,j) = s1(i,j) + dydt*(q3(i,j) - q3(i,j-1))
                q1(i,j) = s1(i,j)/aux(1,i,j)
            enddo
          enddo
        end select  

      !q2
      is = xi; ie = xf
      js = yi; je = yf
      if (gxi ==  1)  is = 2
      if (gxf == xf)  ie = xf-1
      if (gyf == yf)  je = yf-1 
      select case(pml_decomp)
        case(1)
          do j = js,je
             do i = is,ie
              s2(i,j) = pml(5,i,j)*s2(i,j) + pml(6,i,j)*(q3(i-1,j) - q3(i,j))
              q2(i,j) = s2(i,j)/aux(2,i,j)
             enddo
          enddo
        case default
          do j = js,je
            do i = is,ie
              s2(i,j) = s2(i,j) + dxdt*(q3(i-1,j) - q3(i,j))
              q2(i,j) = s2(i,j)/aux(2,i,j)
            enddo
          enddo
        end select
    end select

end subroutine fdtd2D

subroutine fdtdDispersion2D(aux,pml,dxdt,dydt,s1,s2,s3,s4,q1,q2,q3,psum,na,np,xi,xf,yi,yf,gxi,gxf,gyi,gyf,q_type,pml_decomp)

  implicit none
  integer, intent(in) :: na, np, q_type, pml_decomp
  integer, intent(in) :: xi, xf, yi, yf, gxi, gxf, gyi, gyf
  double precision, dimension(na, xi:xf,  yi:yf ),  intent(in)    :: aux
  double precision, dimension(np, xi:xf,  yi:yf ),  intent(in)    :: pml
  double precision, dimension(    xi:xf,  yi:yf ),  intent(inout) :: s1,s2,s3,s4
  double precision, dimension(   gxi:gxf,gyi:gyf),  intent(inout) :: q1,q2,q3
  double precision, dimension(3,  xi:xf, yi:yf  ),  intent(in)    :: psum
  double precision, intent(in) :: dxdt,dydt
  integer :: i,j, is, js, ie, je

  select case(q_type) 
    case(0) !calculate q3
      is = xi; ie = xf
      js = yi; je = yf

      if (gxf == xf) ie = xf-1
      if (gyf == yf) je = yf-1
      select case(pml_decomp)
        case(1)
          do j = js,je
            do i = is,ie
              s3(i,j) = pml(1,i,j)*s3(i,j) + pml(2,i,j)*(q1(i,j+1) - q1(i,j))
              s4(i,j) = pml(3,i,j)*s4(i,j) + pml(4,i,j)*(q2(i,j) - q2(i+1,j))
              q3(i,j) = (s3(i,j)+s4(i,j)-psum(3,i,j))/aux(3,i,j)
            enddo
          enddo
        case default
          do j = js,je
            do i = is,ie
              s3(i,j) = s3(i,j) + dydt*(q1(i,j+1) - q1(i,j))
              s4(i,j) = s4(i,j) + dxdt*(q2(i,j) - q2(i+1,j))
              q3(i,j) = (s3(i,j)+s4(i,j)-psum(3,i,j))/aux(3,i,j)
            enddo
          enddo
        end select
    case(1) !calculate q1 and q2 
      ! q1
      is = xi; ie = xf
      js = yi; je = yf
      if (gxf == xf)  ie = xf-1
      if (gyi ==  1)  js = 2
      if (gyf == yf)  je = yf-1
      select case(pml_decomp)
        case(1)
          do j = js,je
            do i = is,ie
              s1(i,j) = pml(7,i,j)*s1(i,j) + pml(8,i,j)*(q3(i,j) - q3(i,j-1))
              q1(i,j) = (s1(i,j)-psum(1,i,j))/aux(1,i,j)
            enddo
          enddo
        case default
          do j = js,je
            do i = is,ie
                s1(i,j) = s1(i,j) + dydt*(q3(i,j) - q3(i,j-1))
                q1(i,j) = (s1(i,j)-psum(1,i,j))/aux(1,i,j)
            enddo
          enddo
        end select  

      !q2
      is = xi; ie = xf
      js = yi; je = yf
      if (gxi ==  1)  is = 2
      if (gxf == xf)  ie = xf-1
      if (gyf == yf)  je = yf-1 
      select case(pml_decomp)
        case(1)
          do j = js,je
             do i = is,ie
              s2(i,j) = pml(5,i,j)*s2(i,j) + pml(6,i,j)*(q3(i-1,j) - q3(i,j))
              q2(i,j) = (s2(i,j)-psum(2,i,j))/aux(2,i,j)
             enddo
          enddo
        case default
          do j = js,je
            do i = is,ie
              s2(i,j) = s2(i,j) + dxdt*(q3(i-1,j) - q3(i,j))
              q2(i,j) = (s2(i,j)-psum(2,i,j))/aux(2,i,j)
            enddo
          enddo
        end select
    end select

end subroutine fdtdDispersion2D

subroutine CalcDispersion2D(np,q1,q2,q3,c1,c2,c3,p1,p2,p3,psum,xi,xf,yi,yf,gxi,gxf,gyi,gyf,q_type,pml_decomp)

  implicit none
  integer, intent(in) :: np, q_type, pml_decomp
  integer, intent(in) :: xi, xf, yi, yf, gxi, gxf, gyi, gyf
  double precision, dimension(       gxi:gxf,gyi:gyf), intent(in)    :: q1,q2,q3
  double precision, dimension(np, 3,  xi:xf,yi:yf   ), intent(inout) :: p1,p2,p3
  double precision, dimension(3,      xi:xf,yi:yf   ), intent(inout) :: psum
  double precision, dimension(np,     xi:xf,yi:yf   ), intent(in)    :: c1,c2,c3
  integer :: i,j, is, js, ie, je, k

  select case(q_type) 
    case(0) !calculate q3
      is = xi; ie = xf
      js = yi; je = yf

      if (gxf == xf) ie = xf-1
      if (gyf == yf) je = yf-1
      select case(pml_decomp)
        case(1)
          do j = js,je
            do i = is,ie
              psum(3,i,j) = sum(p3(:,1,i,j))
              do k = 1,np
                p3(k,1,i,j) = c1(k,i,j)*p3(k,2,i,j) + c2(k,i,j)*p3(k,3,i,j) + c3(k,i,j)*q3(i,j)
                p3(k,3,i,j) = p3(k,2,i,j)
                p3(k,2,i,j) = p3(k,1,i,j)
              enddo
            enddo
          enddo
        case default
          do j = js,je
            do i = is,ie
              psum(3,i,j) = sum(p3(:,1,i,j))
              do k = 1,np
                p3(k,1,i,j) = c1(k,i,j)*p3(k,2,i,j) + c2(k,i,j)*p3(k,3,i,j) + c3(k,i,j)*q3(i,j)
                p3(k,3,i,j) = p3(k,2,i,j)
                p3(k,2,i,j) = p3(k,1,i,j)
              enddo
            enddo
          enddo
        end select
    case(1) !calculate q1 and q2 
      ! q1
      is = xi; ie = xf
      js = yi; je = yf
      if (gxf == xf)  ie = xf-1
      if (gyi ==  1)  js = 2
      if (gyf == yf)  je = yf-1
      select case(pml_decomp)
        case(1)
          do j = js,je
            do i = is,ie
              psum(1,i,j) = sum(p1(:,1,i,j))
              do k = 1,np
                p1(k,1,i,j) = c1(k,i,j)*p1(k,2,i,j) + c2(k,i,j)*p1(k,3,i,j) + c3(k,i,j)*q1(i,j)
                p1(k,3,i,j) = p1(k,2,i,j)
                p1(k,2,i,j) = p1(k,1,i,j)
              enddo
            enddo
          enddo
        case default
          do j = js,je
            do i = is,ie
              psum(1,i,j) = sum(p1(:,1,i,j))
              do k = 1,np
                p1(k,1,i,j) = c1(k,i,j)*p1(k,2,i,j) + c2(k,i,j)*p1(k,3,i,j) + c3(k,i,j)*q1(i,j)
                p1(k,3,i,j) = p1(k,2,i,j)
                p1(k,2,i,j) = p1(k,1,i,j)
              enddo
            enddo
          enddo
        end select  

      !q2
      is = xi; ie = xf
      js = yi; je = yf
      if (gxi ==  1)  is = 2
      if (gxf == xf)  ie = xf-1
      if (gyf == yf)  je = yf-1 
      select case(pml_decomp)
        case(1)
          do j = js,je
             do i = is,ie
              psum(2,i,j) = sum(p2(:,1,i,j))
              do k = 1,np
                p2(k,1,i,j) = c1(k,i,j)*p2(k,2,i,j) + c2(k,i,j)*p2(k,3,i,j) + c3(k,i,j)*q2(i,j)
                p2(k,3,i,j) = p2(k,2,i,j)
                p2(k,2,i,j) = p2(k,1,i,j)
              enddo
             enddo
          enddo
        case default
          do j = js,je
            do i = is,ie
              psum(2,i,j) = sum(p2(:,1,i,j))
              do k = 1,np
                p2(k,1,i,j) = c1(k,i,j)*p2(k,2,i,j) + c2(k,i,j)*p2(k,3,i,j) + c3(k,i,j)*q2(i,j)
                p2(k,3,i,j) = p2(k,2,i,j)
                p2(k,2,i,j) = p2(k,1,i,j)
              enddo
            enddo
          enddo
        end select
    end select

end subroutine CalcDispersion2D