module evolution_mod
    use healpix_types
    use params
    use time_mod
    use ode_solver
    use rec_mod
    implicit none

    ! Accuracy parameters
    real(dp),     parameter, private :: a_init  = 1.d-8
    real(dp),     parameter, private :: x_init  = log(a_init)
    real(dp),     parameter, private :: kmin    = 0.1d0 * H_0 / c
    real(dp),     parameter, private :: kmax    = 1.d3  * H_0 / c
    integer(i4b), parameter          :: nk      = 100
    !integer(i4b), parameter          :: nf      = 1000
    integer(i4b), parameter, private :: lmax_int = 6

    ! Perturbation quantities
    real(dp), allocatable, dimension(:,:,:) :: Theta
    real(dp), allocatable, dimension(:,:)   :: delta
    real(dp), allocatable, dimension(:,:)   :: delta_b
    real(dp), allocatable, dimension(:,:)   :: Phi
    real(dp), allocatable, dimension(:,:)   :: Psi
    real(dp), allocatable, dimension(:,:)   :: v
    real(dp), allocatable, dimension(:,:)   :: v_b
    real(dp), allocatable, dimension(:,:)   :: dPhi
    real(dp), allocatable, dimension(:,:)   :: dPsi
    real(dp), allocatable, dimension(:,:)   :: dv_b
    real(dp), allocatable, dimension(:,:,:) :: dTheta

    ! Fourier mode list
    real(dp), allocatable, dimension(:) :: ks, xf

    ! Book-keeping variables
    real(dp),     private :: k_current
    integer(i4b), private :: npar = 6 + lmax_int

    real(dp),     private                   :: ck, ckHp, dt, a
    logical(lgt)                            :: firsttime

contains


  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  subroutine get_hires_source_function(k, x, S)
    implicit none

    real(dp), pointer, dimension(:),   intent(out) :: k, x
    real(dp), pointer, dimension(:,:), intent(out) :: S

    integer(i4b) :: i, j
    real(dp)     :: g, dg, ddg, tau, dt, ddt, Hp, dHp, ddHHp, Pi, dPi, ddPi
    real(dp), allocatable, dimension(:,:) :: S_lores

    ! Task: Output a pre-computed 2D array (over k and x) for the 
    !       source function, S(k,x). Remember to set up (and allocate) output 
    !       k and x arrays too. 
    !
    ! Substeps:
    !   1) First compute the source function over the existing k and x
    !      grids
    !   2) Then spline this function with a 2D spline
    !   3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays

  end subroutine get_hires_source_function





! ----------- solving boltzmann-einstein -----------
subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b) :: l, i, k

    ! initializing k-grid
    allocate(ks(nk))
    do i = 1, nk
        ks(i) = kmin + (kmax-kmin) * ((i-1.d0) / (nk-1.d0))**2
    end do

    ! expand x array
    nt = 1000
    allocate(xf(nt))
    do i = 0, nt
        if (i <= 500) then
            xf(i) = x_init + (i-1.d0)*(-log(1.d0+1630.4d0) - x_init)/(500.d0+1.d0)
        else
            xf(i) = xt(i-500)
        end if
    end do

    ! allocating arrays for perturbation quantities
    allocate(Theta(0:nt, 0:lmax_int, nk))
    allocate(delta(0:nt,    nk))
    allocate(delta_b(0:nt,  nk))
    allocate(v(0:nt,        nk))
    allocate(v_b(0:nt,      nk))
    allocate(Phi(0:nt,      nk))
    allocate(Psi(0:nt,      nk))
    allocate(dPhi(0:nt,     nk))
    allocate(dPsi(0:nt,     nk))
    allocate(dv_b(0:nt,     nk))
    allocate(dTheta(0:nt, 0:lmax_int, nk))

    ! initial conditions for the Boltzmann and Einstein equations
    Phi(0,:)     = 1d0
    delta(0,:)   = 1.5d0 * Phi(0,:)
    delta_b(0,:) = delta(0,:)
       
    do i = 1, nk
        !constants
        ckHp         = c*ks(i) / get_Hp(x_init)
        dt           = get_dtau(x_init)

        !initial conditions
        v(0,i)       = ckHp/(2.d0) * Phi(0,i)
        v_b(0,i)     = v(0,i)
        Theta(0,0,i) = 0.5d0 * Phi(0,i)
        Theta(0,1,i) = -ckHp/(6.d0) * Phi(0,i)
        Theta(0,2,i) = -20.d0 * ckHp/(45.d0*dt) * Theta(0,1,i)
        do l = 3, lmax_int
            Theta(0,l,i) = -l/((2.d0*l+1)*dt) * ckHp * Theta(0,l-1,k)
        end do
    end do

end subroutine initialize_perturbation_eqns


! -------------- integrating --------------
subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b)    :: i, j, k, l
    real(dp)        :: x1, x2
    real(dp)        :: eps, hmin, h1, xtc, Hp, dt

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, dydx

    eps    = 1.d-8
    hmin   = 0.d0
    h1     = 1.d-5

    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))

    ! propagate each k-mode independently
    do k = 1, nk
        firsttime   = .true.
        x1          = x_init
        k_current   = ks(k)           ! Store k_current as a global module variable
        ck          = k_current * c

        ! initialize equation set for tight coupling
        y_tight_coupling(1) = delta(0,k)
        y_tight_coupling(2) = delta_b(0,k)
        y_tight_coupling(3) = v(0,k)
        y_tight_coupling(4) = v_b(0,k)
        y_tight_coupling(5) = Phi(0,k)
        y_tight_coupling(6) = Theta(0,0,k)
        y_tight_coupling(7) = Theta(0,1,k)

        ! find the time to which tight coupling is assumed
        xtc = get_tight_coupling_time(k_current)

        do i = 2, nt
            x2 = xf(i)
            if (x1 < xtc) then
            ! integrate from x(init) until the end of tight coupling
                call odeint(y_tight_coupling, x1, x2, eps, h1, hmin, derivstc, bsstep, output)
                x1 = x2

                dt = get_dtau(x1)
                Hp = get_Hp(x1)
                ckHp = ck*Hp

                delta(i,k)      = y_tight_coupling(1)
                delta_b(i,k)    = y_tight_coupling(2)
                v(i,k)          = y_tight_coupling(3)
                v_b(i,k)        = y_tight_coupling(4)
                Phi(i,k)        = y_tight_coupling(5)
                Theta(i,0,k)    = y_tight_coupling(6)
                Theta(i,1,k)    = y_tight_coupling(7)
                Theta(i,2,k)    = -20.d0 * ckHp/(45.d0*dt) * Theta(i,1,k)

                do l = 3, lmax_int
                    Theta(i,l,k) = -l/(2.d0*l+1.d0) * ckHp/dt * Theta(i,l-1,k)
                end do

                Psi(i,k) = -Phi(i,k) - 12.d0*H_0**2 / (ck*exp(xf(i)))**2.d0 * (omega_r*Theta(i,2,k))

            else
            ! integrate from end of tight coupling until today
                if (firsttime) then
                    dt = get_dtau(x1)
                    Hp = get_Hp(x1)
                    ckHp = ck*Hp

                    y(1:7)  = y_tight_coupling(1:7)
                    y(8)    = -20.d0 * ckHp/45.d0/dt * y(7)     ! theta2
                    do l = 3, lmax_int
                        y(6+l)  = -l/(2.d0*l+1.d0) * ckHp/dt * y(6+l-1)
                    end do
                    firsttime = .false.
                end if

                call odeint(y, x1 ,x2, eps, h1, hmin, derivs, bsstep, output)
                x1 = x2

                ! store variables at i
                delta(i,k)   = y(1)
                delta_b(i,k) = y(2)
                v(i,k)       = y(3)
                v_b(i,k)     = y(4)
                Phi(i,k)     = y(5)
                do l = 0, lmax_int
                    Theta(i,l,k) = y(6+l)
                end do
                Psi(i,k) = -y(5) - 12.d0*H_0**2 / (ck*exp(xf(i)))**2.d0 * omega_r * y(8)

                ! store derivatives required for Cl
                call derivs(xf(i), y, dydx)

                dv_b(i,k)     = dydx(4)
                dPhi(i,k)     = dydx(5)
                do l = 0, lmax_int
                    dTheta(i,l,k) = dydx(6+l)
                end do
                dPsi(i,k) = -dPhi(i,k) - 12.d0*H_0**2.d0 / (ck*exp(xf(i)))**2.d0 * omega_r * (-2.d0*Theta(i,2,k) + dTheta(i,2,k))
            end if
        end do
    end do

    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)


end subroutine integrate_perturbation_eqns


! ------------- subroutines for odeint -------------
! from x_init to end of tight coupling
subroutine derivstc(x,ytc,dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: ytc
    real(dp), dimension(:), intent(out) :: dydx

    real(dp) :: ddelta, ddelta_b
    real(dp) :: dv, dv_b
    real(dp) :: q, R
    real(dp) :: delta, delta_b, v, v_b, Phi, Theta0, Theta1, Theta2
    real(dp) :: Psi, dPhi, dTheta0, dTheta1
    real(dp) :: dt, ddt, a, Hp, dHp, ckHp

    ! store variables for derivation
    delta   = ytc(1)
    delta_b = ytc(2)
    v       = ytc(3)
    v_b     = ytc(4)
    Phi     = ytc(5)
    Theta0  = ytc(6)
    Theta1  = ytc(7)

    ! calculate derivatives
    dt      = get_dtau(x)
    ddt     = get_ddtau(x)
    a       = exp(x)
    Hp      = get_Hp(x)
    dHp     = get_dHp(x)
    ckHp    = ck/Hp

    Theta2    = -20.d0*ckHp / (45.d0*dt)*Theta1
    R         = (4.d0*omega_r) / (3.d0*omega_b*a)
    Psi       = -Phi - 12.d0*(H_0/(ck*a))**2.d0 * omega_r*Theta2
    dPhi      = Psi - (ckHp**2.d0)/3.d0*Phi + H_0**2.d0 / (2.d0*Hp**2.d0) * (omega_m/a*delta + omega_b/a*delta_b + 4.d0*omega_r * Theta0/a**2.d0)
    dTheta0   = -ckHp*Theta1 - dPhi
    ddelta    = ckHp*v - 3.d0*dPhi
    ddelta_b  = ckHp*v_b - 3.d0*dPhi
    dv        = -v - ckHp * Psi
    q         = (-((1.d0-2.d0*R)*dt + (1.d0+R)*ddt)*(3.d0*Theta1+v_b) - ckHp*Psi + (1.d0-dHp/Hp)*ckHp*(-Theta0+2.d0*Theta2) - ckHp*dTheta0) / ((1.d0+R)*dt+dHp/Hp -1.d0)
    dv_b      = (1.d0/(1.d0+R)) * (-v_b - ckHp*Psi + R*(q+ckHp*(-Theta0+2.d0*Theta2) - ckHp*Psi))
    dTheta1   = (1.d0/3.d0) * (q - dv_b)

    ! store derivatives
    dydx(1) = ddelta
    dydx(2) = ddelta_b
    dydx(3) = dv
    dydx(4) = dv_b
    dydx(5) = dPhi
    dydx(6) = dTheta0
    dydx(7) = dTheta1

    !write(*,*) ytc,dydx
    !stop
end subroutine derivstc


! from end of tight coupling until today
subroutine derivs(x,y,dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx

    integer(i4b) :: l
    real(dp) :: ddelta, ddelta_b, dv, R
    real(dp) :: delta, delta_b, v, v_b, Phi
    real(dp) :: Theta0, Theta1, Theta2, Theta3, Theta4, Theta5, Theta6
    real(dp) :: Psi, dPhi, dv_b, dTheta0, dTheta1, dTheta2
    real(dp) :: a, Hp, ckHp, dt

    ! store variables for derivation
    delta   = y(1)
    delta_b = y(2)
    v       = y(3)
    v_b     = y(4)
    Phi     = y(5)
    Theta0  = y(6)
    Theta1  = y(7)
    Theta2  = y(8)
    Theta3  = y(9)
    Theta4  = y(10)
    Theta5  = y(11)
    Theta6  = y(12)

    ! calculate derivatives
    a = exp(x)
    Hp = get_Hp(x)
    ckHp = ck/Hp
    dt = get_dtau(x)

    R         = (4.d0*omega_r) / (3.d0*omega_b*a)
    Psi       = -Phi - 12.d0*(H_0/(ck*a))**2.d0 * omega_r * Theta2
    dPhi      = Psi - (ckHp**2.d0)/3.d0 * Phi + H_0**2.d0 / (2.d0*Hp**2.d0) * (omega_m/a*delta + omega_b/a*delta_b + 4.d0*omega_r * Theta0/a**2.d0)
    dTheta0   = -ckHp * Theta1 - dPhi
    ddelta    = ckHp * v - 3.d0 * dPhi
    ddelta_b  = ckHp * v_b - 3.d0 * dPhi
    dv        = -v - ckHp * Psi
    dv_b      = -v_b - ckHp * Psi + dt*R*(3.d0*Theta1+v_b)
    dTheta1   = ckHp/3.d0*Theta0 - 2.d0/3.d0*ckHp*Theta2 + ckHp/3.d0*Psi + dt*(Theta1+v_b/3.d0)
    !dTheta2   = 2.d0/5.d0*ckHp*Theta1 - 3.d0/5.d0*ckHp*Theta3 + dt*1.d0/10.d0*Theta2

    dydx(1) = ddelta
    dydx(2) = ddelta_b
    dydx(3) = dv
    dydx(4) = dv_b
    dydx(5) = dPhi
    dydx(6) = dTheta0
    dydx(7) = dTheta1
    !dydx(8) = dTheta2

    ! dTheta2,3,4,5
    do l = 2, lmax_int-1
        dydx(6+l) = l*ckHp/(2.d0*l+1.d0) * y(6+l-1) - (l+1.d0)*ckHp/(2.d0*l+1.d0) * y(6+l+1) + dt*(y(6+l) - 1.d0/10.d0*y(6+l)*abs(l==2))    !true=-1
    end do
    ! dTheta6
    dydx(6+lmax_int) = ckHp*y(6+lmax_int-1) - c*(lmax_int+1.d0)/(Hp*get_eta(x))*y(6+lmax_int) + dt*y(6+lmax_int)

    !write(*,*) y(8), y
    !stop


  end subroutine derivs




! Task: Complete the following routine, such that it returns the time at which
!       tight coupling ends. In this project, we define this as either when
!       dtau < 10 or c*k/(Hp*dt) > 0.1 or x > x(start of recombination)

! ------------- functions -------------
function get_tight_coupling_time(k)
    implicit none

    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time
    integer(i4b)          :: i
    real(dp)              :: x, Hp, dt, x_start_rec

    x_start_rec = -log(1.d0+1630.4d0)
    x = x_init

    do i = 1, nt
        x = x + (x_start_rec-x_init)/(nt-1)
        Hp = get_Hp(x)
        dt = get_dtau(x)
        if ((x < x_start_rec) &
        .and. (abs(c*k/(Hp*dt)) <= 0.1d0) &
        .and. abs(dt) > 10.d0) &
        then
            get_tight_coupling_time = x
        end if
    end do

end function get_tight_coupling_time



! ------------- write all -------------
subroutine write_mk3
    use healpix_types
    implicit none

    CHARACTER(*), PARAMETER :: loc = "/uio/hume/student-u06/kariandy/AST5220/code/data/"

    integer(i4b)                             :: i
    integer(i4b), allocatable, dimension(:)  :: kw
    allocate(kw(nk))
    ! choose k values that show all 3 regimes
    kw(1:6) = (/ 1, 12, 30, 40, 85, 100 /)
    !kw(1:6) = (/ 1, 5, 10, 40, 60, 100 /)

    open(1, file = loc//"phi.dat")
    open(2, file = loc//"psi.dat")
    open(3, file = loc//"delta.dat")
    open(4, file = loc//"delta_b.dat")
    open(5, file = loc//"v.dat")
    open(6, file = loc//"v_b.dat")
    open(7, file = loc//"theta0.dat")
    open(8, file = loc//"theta1.dat")
    open(9, file = loc//"k_values.dat")

    do i=1,nt
        write(1,'(7F20.7)') xf(i), Phi(i,kw(1)), Phi(i,kw(2)), Phi(i,kw(3)), Phi(i,kw(4)), Phi(i,kw(5)), Phi(i,kw(6))
        write(2,'(7F20.7)') xf(i), Psi(i,kw(1)), Psi(i,kw(2)), Psi(i,kw(3)), Psi(i,kw(4)), Psi(i,kw(5)), Psi(i,kw(6))
        write(3,'(7F20.7)') xf(i), delta(i,kw(1)), delta(i,kw(2)), delta(i,kw(3)), delta(i,kw(4)), delta(i,kw(5)), delta(i,kw(6))
        write(4,'(7F20.7)') xf(i), delta_b(i,kw(1)), delta_b(i,kw(2)), delta_b(i,kw(3)), delta_b(i,kw(4)), delta_b(i,kw(5)), delta_b(i,kw(6))
        write(5,'(7F20.7)') xf(i), v(i,kw(1)), v(i,kw(2)), v(i,kw(3)), v(i,kw(4)), v(i,kw(5)), v(i,kw(6))
        write(6,'(7F20.7)') xf(i), v_b(i,kw(1)), v_b(i,kw(2)), v_b(i,kw(3)), v_b(i,kw(4)), v_b(i,kw(5)), v_b(i,kw(6))
        write(7,'(7F20.7)') xf(i), Theta(i,0,kw(1)), Theta(i,0,kw(2)), Theta(i,0,kw(3)), Theta(i,0,kw(4)), Theta(i,0,kw(5)), Theta(i,0,kw(6))
        write(8,'(7F20.7)') xf(i), Theta(i,1,kw(1)), Theta(i,1,kw(2)), Theta(i,1,kw(3)), Theta(i,1,kw(4)), Theta(i,1,kw(5)), Theta(i,1,kw(6))
    end do

    do i = 1,6
        write(9,*) ks(i), kw(i)
    end do

    ! close files
    do i=1,9
        close(i)
    end do

    deallocate(kw)

end subroutine write_mk3


end module evolution_mod
