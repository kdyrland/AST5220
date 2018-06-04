module cl_mod
    use healpix_types
    use evolution_mod
    use sphbess_mod
    implicit none

    integer(i4b), allocatable,  dimension(:)        :: ls
    real(dp),     pointer,      dimension(:)        :: khires, xhires
    real(dp),     allocatable,  dimension(:,:)      :: Thetal, Thetal2
    real(dp),     allocatable,  dimension(:)        :: lsdp, lhires, clhires
    real(dp),     allocatable,  dimension(:)        :: integrand, integrandx
    real(dp),     allocatable,  dimension(:,:)      :: integrandk


contains

! finally computing the CMB power spectrum
subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, k, l, lnum, nhires, nspline
    integer(i4b), allocatable, dimension(:)         :: lw
    real(dp),     pointer,      dimension(:,:)      :: jl, jl2
    real(dp),     pointer,      dimension(:)        :: x_arg, int_arg, cls, cls2
    real(dp),     pointer,      dimension(:,:,:,:)  :: S_coeff
    real(dp),     pointer,      dimension(:,:)      :: S, S2
    real(dp),     allocatable,  dimension(:)        :: z_spline, jl_spline, jl_spline2
    real(dp)                                        :: integralx, integralk, h1, h2
    real(dp)                                        :: start, finish
    logical(lgt)       :: exist
    character(len=128) :: filename
    CHARACTER(*), PARAMETER :: loc = "/uio/hume/student-u06/kariandy/AST5220/code/data"

    ! set up which l's to compute
    lnum = 44
    allocate(ls(lnum))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)


    ! ------- get source function from evolution_mod -------
    nhires = 5000

    write(*,*) "We call"
    call get_hires_source_function(xhires, khires, S)



    ! ------- calculate spherical bessel functions --------
    nspline = 5400
    allocate(z_spline(nspline))    ! j_l(z)
    allocate(jl(nspline, lnum))
    allocate(jl2(nspline, lnum))
    allocate(integrand(nhires))

    write(*,*) "We bessel"

    ! initialize array
    do i = 1,nspline
        z_spline(i) = (i-1.d0)*3500.d0/(nspline-1.d0)
    end do

    ! check for binary file
    filename = 'jl.bin'
    inquire(file = filename, exist = exist)
    if (exist) then
        ! read
        write(*,*) "We read"
        open(10, form = 'unformatted', file = filename)
        read(10) jl, jl2
        close(10)
    else
        ! compute spherical bessel functions
        write(*,*) "We find"
        do i = 1,nspline
            do l = 1,lnum
                if (z_spline(i) > 0.d0) then
                    call sphbes(ls(l),z_spline(i),jl(i,l))
                end if
            end do
        end do
        ! spline
        write(*,*) "We spline"
        do l = 1,lnum
            call spline(z_spline, jl(:,l), 1.d30, 1.d30, jl2(:,l))
        end do
        ! write to cache
        write(*,*) "We cache"
        open(10, form = 'unformatted', file = filename)
        write(10) jl, jl2
        close(10)
    end if


    ! write integrand to file for sanity checks
    write(*,*) "We sanity check"
    j = locate_dp(khires,340.d0*H_0/c)          ! find k = 340*H/c = k(1700)
    l = 17                                      ! ls(17) = 100

    do i = 1,nhires
        integrand(i) = S(i,j) * splint(z_spline, jl(:,l), jl2(:,l), khires(j)*(get_eta(0.d0)-get_eta(xhires(i))))
    end do

!    write(*,*) "We write test"
!    open(34,file = loc//"integrand.dat", action="write", status="replace")
!    do i = 1,nhires
!        write(34,*) integrand(i)
!    end do
!    close(34)


    ! -------- compute the C_l's for each given l --------

    allocate(Thetal(lnum,nhires))
    allocate(integrandx(nhires))
    allocate(integrandk(lnum,nhires))
    allocate(cls(lnum))
    allocate(cls2(lnum))

    ! for integration
    h1 = (xhires(nhires) - xhires(1))/nhires
    h2 = (khires(nhires) - khires(1))/nhires

    write(*,*) "we find cls"

!    open (unit=123, file="integrand1.dat", action="write", status="replace")
!    open (unit=124, file="integrand2.dat", action="write", status="replace")
!    open (unit=125, file="integrand3.dat", action="write", status="replace")
!    open (unit=126, file="integrand4.dat", action="write", status="replace")
!    open (unit=127, file="integrand5.dat", action="write", status="replace")
!    open (unit=128, file="integrand6.dat", action="write", status="replace")

    do l = 1, lnum
        write(*,*) "l = ", l
        ! compute transfer function
        integralk = 0                   ! reset for each cl
        do k = 1, nhires
            integralx = 0               ! reset for each theta
            ! integrate theta
            do i = 1, nhires
                integrandx(i) = S(i,k) * splint(z_spline,jl(:,l),jl2(:,l),khires(k)*(get_eta(0.d0)-get_eta(xhires(i))))
                integralx = integralx + integrandx(i)
            end do
            ! store theta
            Thetal(l,k) = h1 * (integralx - 0.5d0*(integrandx(1)+integrandx(nhires)))

            ! integrate cl
            integrandk(l,k) = (c*khires(k)/H_0)**(n_s-1.d0) * Thetal(l,k)**2/khires(k)
            integralk = integralk + integrandk(l,k)
        end do

        ! write integrand to file
!        if  (ls(l) == 6) then
!        write(*,*) "writing int 1"
!            do j = 1,nhires
!                write(123,'(*(2X, ES14.6E3))') integrandk(l,j)
!            end do
!        end if
!
!        if  (ls(l) == 100) then
!        write(*,*) "writing int 2"
!            do j = 1,nhires
!                write(124,'(*(2X, ES14.6E3))') integrandk(l,j)
!            end do
!        end if
!
!        if  (ls(l) == 200) then
!        write(*,*) "writing int 3"
!            do j = 1,nhires
!                write(125,'(*(2X, ES14.6E3))') integrandk(l,j)
!            end do
!        end if
!
!        if  (ls(l) == 500) then
!        write(*,*) "writing int 4"
!            do j = 1,nhires
!                write(126,'(*(2X, ES14.6E3))') integrandk(l,j)
!            end do
!        end if
!
!        if  (ls(l) == 800) then
!            write(*,*) "writing int 5"
!            do j = 1,nhires
!                write(127,'(*(2X, ES14.6E3))') integrandk(l,j)
!            end do
!        end if
!
!        if  (ls(l) == 1200) then
!        write(*,*) "writing int 6"
!            do j = 1,nhires
!                write(128,'(*(2X, ES14.6E3))') integrandk(l,j)
!            end do
!        end if

        ! store cl
        integralk = h2 * (integralk - 0.5d0*(integrandk(l,1)+integrandk(l,nhires)))
        cls(l) = integralk * ls(l) * (ls(l)+1.d0) / (2.d0*pi)

        ! timer
        call cpu_time(finish)
        print '("time = ",f7.2," sec")', finish-start
    end do

!    open (unit=121, file="xkhires.dat", action="write", status="replace")
!    do i = 1,nhires
!        write(121,'(*(2X, ES14.6E3))') xhires(i), khires(i)
!    end do
!
!    close(121)
!    close(123)
!    close(124)
!    close(125)
!    close(126)
!    close(127)
!    close(128)

    allocate(lsdp(lnum))
    allocate(clhires(int(maxval(ls))))
    allocate(lhires(int(maxval(ls))))

    ! get double precision
    do l = 1,lnum
        lsdp(l) = ls(l)
    end do

    ! spline cl
    write(*,*) "we spline cl"
    call spline(lsdp, cls, 1.d30, 1.d30, cls2)


    write(*,*) "we highres l"
    do l = 1,ls(lnum)
        lhires(l) = lsdp(1) + (l-1.d0) * (lsdp(lnum)-lsdp(1)) / (lsdp(lnum)-2.d0)
    end do

    ! write cl
    write(*,*) "we highres cl + write"
    open(unit=129,file = "cls.dat")
    do l = 1,ls(lnum)
        clhires(l) = splint(lsdp, cls, cls2, lhires(l))
        write(129,'(7F20.7)') lhires(l), clhires(l)
    end do
    close(129)

    ! ------------- write --------------

    !call write_mk4         ! does not work :(

end subroutine compute_cls


! ------------- write all -------------
subroutine write_mk4    ! does not work
    use healpix_types
    implicit none

    integer(i4b)                                :: i
    integer(i4b), allocatable, dimension(:)     :: lw
    CHARACTER(*), PARAMETER :: loc = "/uio/hume/student-u06/kariandy/AST5220/code/data"

    write(*,*) "we write mk4"

    lw(1:6) = (/6, 100, 200, 500, 800, 1200/)
    !l(1:6) = (/4, 17, 22, 30, 36, 44/)  ! locations in ls

    ! open files
    open(1, file = loc//"lvalues.dat")
    open(2, file = loc//"xkhires.dat")
    open(3, file = loc//"thetals.dat")
    open(4, file = loc//"cl_ints.dat")
    open(5, file = loc//"cls.dat")

    ! write files
    do i = 1,6
        write(*,*) "we write l"
        write(1,*) lw(i)
    end do
    do i = 1,nhires
        write(*,*) "we write x,theta,ints"
        write (2,'(7F20.7)') xhires(i), khires(i)
        write (3,'(7F20.7)') &
            Thetal(lw(1), i),Thetal(lw(2),i),Thetal(lw(3),i),Thetal(lw(4),i),Thetal(lw(5),i),Thetal(lw(6),i)
        write (4,'(7F20.7)') &
            integrandk(lw(1),i) / (c*khires(i)/H_0)**(n_s-1.d0), &
            integrandk(lw(2),i) / (c*khires(i)/H_0)**(n_s-1.d0), &
            integrandk(lw(3),i) / (c*khires(i)/H_0)**(n_s-1.d0), &
            integrandk(lw(4),i) / (c*khires(i)/H_0)**(n_s-1.d0), &
            integrandk(lw(5),i) / (c*khires(i)/H_0)**(n_s-1.d0), &
            integrandk(lw(6),i) / (c*khires(i)/H_0)**(n_s-1.d0)
    end do
    do i = 1,1200
        write(*,*) "we write cls"
        write (5,'(7F20.7)') lhires(i), clhires(i)
    end do

    ! close files
    do i = 1,5
        close(i)
    end do

    deallocate(lw)

end subroutine write_mk4

end module cl_mod
