module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b)                        :: n                      ! grid steps
  real(dp), allocatable, dimension(:) :: xrec                   ! conformal times
  real(dp), allocatable, dimension(:) :: Xe                     ! fractional electron density ne/nH
  real(dp), allocatable, dimension(:) :: Hrec                   ! hubble
  real(dp), allocatable, dimension(:) :: tau, tau2, tau22       ! splined tau and second derivatives
  real(dp), allocatable, dimension(:) :: dtau,logtau            ! first derivative of tau
  real(dp), allocatable, dimension(:) :: ne, ne2                ! splined (log of) electron density ne
  real(dp), allocatable, dimension(:) :: g, g2, g22             ! splined visibility function and second derivatives
  real(dp), allocatable, dimension(:) :: dg                     ! first derivative of g
  real(dp)                            :: yp1, ypn, eps, hmin

contains

  subroutine initialize_rec_mod
    implicit none

    CHARACTER(*), PARAMETER :: fileplace = "/uio/hume/student-u06/kariandy/AST5220/code/data/"

    integer(i4b) :: i, j, k
    real(dp)     :: saha_lim, y, Tb, nb, dydx, xmin, xmax, dx, f, ne0, Xe0, xstart, xstop, xstep, hstep
    real(dp)     :: C_r
    logical(lgt) :: use_saha

    saha_lim    = 0.99d0         ! switch from saha to peebles when Xe < 0.99
    xstart      = log(1.d-10)    ! start at a = 1e-10
    xstop       = 0.d0           ! end at a = 1
    n           = 1000           ! grid points

    ! spline variables
    yp1 = 1.d30
    ypn = 1.d30

    ! integration variables
    eps = 1.d-10
    hmin = 0.d0

    allocate(xrec(n))
    allocate(Xe(n))
    allocate(tau(n))
    allocate(dtau(n))
    allocate(logtau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(ne(n))
    allocate(ne2(n))
    allocate(g(n))
    allocate(dg(n))
    allocate(g2(n))
    allocate(g22(n))

! ---------------- x_rec grid ----------------

    xstep = (xstop-xstart)/(n-1)
    xrec(1) = xstart
    do i = 2,n
        xrec(i) = xrec(i-1) + xstep
    end do

    hstep = abs(1.d-3*(xrec(1)-xrec(2)))     ! ode step length

! ---------------- electron fraction and density ----------------

    use_saha = .true.
    do i = 1, n
        nb = omega_b * rho_c / (m_H*exp(xrec(i))**3)
        if (use_saha) then
        ! use saha
            Tb = T0 / exp(xrec(i))
            Xe0= ((m_e*kb*Tb) / (2.d0*pi*hbar**2))**1.5d0 * exp(-eps_0/(kb*Tb))/nb
            Xe(i) = (-Xe0 + sqrt(Xe0**2 +4.d0*Xe0))/2.d0
            if (Xe(i) < saha_lim) use_saha = .false.
        else
        ! use peebles
            Xe(i) = Xe(i-1)
            call odeint(Xe(i:i), xrec(i-1), xrec(i), eps, hstep, hmin, dXedx, bsstep, output)
        end if
        ne(i) = Xe(i) * nb
    end do

    ! spline logarithmic electron density
    ne = log(ne)
    call spline(xrec, ne, yp1, ypn, ne2)


! ---------------- optical depth tau ----------------

    tau(n) = 0.d0

    do i = n-1,1,-1
        tau(i) = tau(i+1)
        call odeint(tau(i:i), xrec(i+1), xrec(i), eps, hstep, hmin, dtaudx, bsstep, output)
!        logtau(i) = log(tau(i))
    end do

    ! spline optical depth and second derivative
    call spline(xrec, tau, yp1, ypn, tau2)
    call spline(xrec, tau2, yp1, ypn, tau22)

    ! dtau
    do i = 1,n
        dtau(i) = get_dtau(xrec(i))
    end do

! or
!    do i = 1,n
!        Hrec(i) = get_H(xrec(i))
!    end do
!    do i = 1,n
!        dtau(i) = -ne(i) * sigmaT * c / Hrec(i)
!    end do



! ---------------- visibility function g ----------------

    do i=1,n
        g(i) = -dtau(i) * exp(-tau(i))
    end do

    ! spline g and second derivative
    call spline(xrec,g,yp1,ypn,g2)
    call spline(xrec,g2,yp1,ypn,g22)

    ! save dg
    do i = 1,n
        dg(i) = get_dg(xrec(i))
    end do







! ---------------- write all ----------------

    open(1, file=fileplace//"Xe.dat")
    open(2, file=fileplace//"tau.dat")
    open(3, file=fileplace//"gfunc.dat")
        do i=1, n
            write(1,'(2(ES20.7E3))') xrec(i), Xe(i)
            write(2,'(4(ES20.7E3))') xrec(i), tau(i), tau2(i), dtau(i)
            write(3,'(4F20.7)') xrec(i), g(i), g2(i), dg(i)
        end do
    close(1)
    close(2)
    close(3)

end subroutine initialize_rec_mod




! -------------- derivatives for odeint --------------

subroutine dXedx(x, Xe, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: Xe
    real(dp), dimension(:), intent(out) :: dydx
    real(dp)                            :: Tb, nb, phi2, alpha2, beta, beta2, n1s, lambda_2s1s, lambda_alpha, Cr, H
    H      = get_H(x)
    Tb     = T0 / exp(x)
    nb     = omega_b * rho_c / (m_H*exp(x)**3)
    phi2   = 0.448d0 * log(eps_0/(kb*Tb))
    alpha2 = 64.d0 * pi / sqrt(27.d0*pi) * (alpha/m_e)**2 * sqrt(eps_0/(kb*Tb)) * phi2 * hbar**2 /c
    beta   = alpha2 *((m_e*kb*Tb)/(2.d0*pi*hbar**2))**1.5 * exp(-eps_0/(kb*Tb))

        ! to prevent beta2 going to infinity, set to 0
    if (Tb <= 169.d0) then
        beta2 = 0.d0
    else
        beta2 = beta * exp((3.d0*eps_0)/(4.d0*kb*Tb))
    end if

    n1s             = (1.d0-Xe(1)) * nb
    lambda_2s1s     = 8.227d0
    lambda_alpha    = H * (3.d0*eps_0)**3 / ((8.d0*pi)**2 * n1s) /(c*hbar)**3
    Cr              = (lambda_2s1s + lambda_alpha) / (lambda_2s1s + lambda_alpha + beta2)

    dydx            = Cr/H * (beta*(1.d0-Xe(1)) - nb*alpha2*Xe(1)**2)

end subroutine dXedx


subroutine dtaudx(x, tau, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: tau
    real(dp), dimension(:), intent(out) :: dydx
    real(dp)                            :: ne
    real(dp)                            :: H

    ne      = get_ne(x)
    H       = get_H(x)
    dydx    = -ne * sigmaT / H * c

end subroutine dtaudx




! -------------- functions --------------

  ! electron density
  function get_ne(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_ne

    get_ne = splint(xrec, ne, ne2, x)
    get_ne = exp(get_ne)

  end function get_ne

  ! tau
  function get_tau(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_tau

    get_tau = splint(xrec, tau, tau2, x)
 !   get_tau = exp(get_tau)

  end function get_tau

  ! derivative of tau
  function get_dtau(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_dtau
    real(dp)             :: ne, Hp
    Hp         = get_Hp(x)
    ne         = get_ne(x)
    get_dtau   = -ne * sigmaT * exp(x) * c/Hp
!    get_dtau    =  splint_deriv(xrec, tau, tau2, x)

  end function get_dtau

  ! second derivative of tau
  function get_ddtau(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_ddtau

    get_ddtau = splint(xrec,tau2, tau22, x)

  end function get_ddtau

  ! visibility function, g
  function get_g(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_g

    get_g = splint(xrec, g, g2, x)

  end function get_g

  ! derivative of g
  function get_dg(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_dg

    get_dg = splint_deriv(xrec,g,g2,x)

  end function get_dg

  ! second derivative of g
  function get_ddg(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_ddg

    get_ddg = splint(xrec,g2,g22,x)

  end function get_ddg


end module rec_mod
