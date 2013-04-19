subroutine build_pml(pml,num_pml,pml_type,ddx,ddy,dt,norder,Ro,co,xi,xf,yi,yf,gxi,gxf,gyi,gyf,nx,ny)

    implicit none
    integer, intent(in) :: num_pml,pml_type,norder,nx,ny
    integer, intent(in) :: xi, xf, yi, yf, gxi, gxf, gyi, gyf
    double precision, intent(in) :: ddx, ddy,dt,Ro,co
    double precision, dimension(num_pml,gxi:gxf,gyi:gyf), intent(inout) :: pml
    double precision, dimension(num_pml) :: sigex, sigmx, sigey, sigmy, cax, cbx, cay, cby, dax, dbx, day, dby
    integer :: i,j, is, js, ie, je
    integer :: m
    
!   precalculations
    sigmex = -(norder + 1.0)*co*log(Ro)/(2.0*ddx*num_pml)
    sigmmx = sigmex
    sigmey = -(norder + 1.0)*co*log(Ro)/(2.0*ddy*num_pml)
    sigmmy = sigmey

    do m = 1,num_pml
        sigex(m) = sigmex*((m - 0.5)/(num_pml + 0.5))**norder
        sigmx(m) = sigmmx*(m/(num_pml + 0.5))**norder
        sigey(m) = sigmey*((m - 0.5)/(num_pml + 0.5))**norder
        sigmy(m) = sigmmy*(m/(num_pml + 0.5))**norder
    end do

    do m = 1,num_pml
        rex = sigex(m)*dt
        rmx = sigmx(m)*dt
        cax(m) = exp(-rex)
        cbx(m) = -(exp(-rex) - 1.0)/sigex(m)/ddx
        dax(m) = exp(-rmx)
        dbx(m) = -(exp(-rmx) - 1.0)/sigmx(m)/ddx
        rey = sigey(m)*dt
        rmy = sigmy(m)*dt
        cay(m) = exp(-rey)
        cby(m) = -(exp(-rey) - 1.0)/sigey(m)/ddy
        day(m) = exp(-rmy)
        dby(m) = -(exp(-rmy) - 1.0)/sigmy(m)/ddy
    end do

    pml(:,:,:) = 1.0
    pml(2,:,:) = dt/ddy
    pml(4,:,:) = dt/ddx
    pml(6,:,:) = dt/ddx
    pml(8,:,:) = dt/ddy

!   fill pml array
    select case(pml_type)

        case(0) !left boundary PML
            ! left q2
            is = 1; ie = num_pml
            js = yi; je = yf
            if (gyf == yf) je = yf-1
            do j = js,je
                do i = is,ie
                    pml(3,i,j) = day(i)
                    pml(4,i,j) = dby(i)
                end do 
            end do

            is = xi; ie = num_pml+1
            js = yi; je = yf
            if (gxi ==  1)  is = 2
            if (gyf == yf)  je = yf-1 
            do j = js,je
                do i = is,ie
                    pml(5,i,j) = cay(i-1)
                    pml(6,i,j) = cby(i-1)
                end do
            end do

        case(1) !right boundary 
        is = nx - num_pml, ie = nx - 1
        js = yi; je = yf
        if (gyf == yf) je = yf - 1
        do j = js,je
            do i = is,ie
                m = i - nx + num_pml + 1
                pml(3,i,j) = day(m)
                pml(4,i,j) = dby(m)
                pml(5,i,j) = cay(m)
                pml(6,i,j) = cby(m)
            end do 
        end do

        case(2) !bottom boundary
        is = xi; ie = xf
        js = 1; je = num_pml
        if (gxf == xf) ie = xf-1
        do j = js,je
            do i = is, ie 
                pml(1,i,j) = dax(i)
                pml(2,i,j) = dbx(i)
            end do
        end do

        is = xi; ie = xf
        js = 1; je = num_pml + 1
        if (gxf == xf)  ie = xf-1
        if (gyi ==  1)  js = 2

        do j = js,je
            do i = is, ie 
                pml(7,i,j) = cax(i-1)
                pml(8,i,j) = cbx(i-1)
            end do
        end do

        case(3) !top boundary
        is = xi; ie = xf
        js = ny - num_pml; je = ny - 1
        if (gxf == xf)  ie = xf-1
        do j = js, je
            do i = is, ie
                m = j - ny + num_pml + 1
                pml(1,i,j) = dax(m)
                pml(2,i,j) = dbx(m)
                pml(7,i,j) = cax(m)
                pml(8,i,j) = cbx(m)
            end do
        end do

        case default
            pml(:,:,:) = 1.0
            pml(2,:,:) = dt/ddy
            pml(4,:,:) = dt/ddx
            pml(6,:,:) = dt/ddx
            pml(8,:,:) = dt/ddy

    end select
end subroutine build_pml