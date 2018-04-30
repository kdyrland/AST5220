module time_mod
  use healpix_types
  use params
  use spline_1D_mod
  use ode_solver
  implicit none

  integer(i4b)                           :: nt                ! length x
  real(dp),    allocatable, dimension(:) :: xt                ! conformal time grid
  real(dp),    allocatable, dimension(:) :: at                ! scale factor grid

  integer(i4b)                           :: n_eta             ! length eta grid
  real(dp),    allocatable, dimension(:) :: xeta              ! eta grid
  real(dp),    allocatable, dimension(:) :: eta, eta2         ! eta and ddeta at each grid point

contains

subroutine initialize_time_mod
  implicit none
  CHARACTER(*), PARAMETER :: fileplace = "/uio/hume/student-u06/kariandy/AST5220/code/data/"

  integer(i4b) :: i, n, n1, n2, n3
  real(dp)     :: zri, zrf, z0, xri, xrf, x0, xeta1, xeta2, eta1, a1, a2, xstep, eps, xstepmin, yp1, ypn, rho_crit, rho_m, rho_b, rho_r, rho_v, z,xinit, eta_init

    ! define grid of two epochs
    n1          = 200                       ! grid values
    n2          = 300
    nt          = n1 + n2

    ! redshift z
    zri         = 1630.4d0                  ! start rec
    zrf         = 614.2d0                   ! end rec
    z0          = 0.d0                      ! today

    ! conformal time ln(a) = ln(1+z)
    xri         = -log(1.d0 + zri)
    xrf         = -log(1.d0 + zrf)
    x0          = 0.d0

    ! values for evaluating conformal time
    n_eta       = 1000                      ! grid points for spline
    a1          = 1.d-10                    ! scale factor, start
    a2          = 1.d0                      ! scale factor, end
    xeta1       = log(a1)
    xeta2       = 0.d0
    xinit       = xeta1

    yp1         = 1.d30
    ypn         = 1.d30
    xstepmin    = 0.d0
    eps         = 1.d-10
    eta_init    = a1 / (H_0*sqrt(omega_r))

  ! Task: Fill in x and a grids
  allocate(xt(nt))
  allocate(at(nt))

  xstep = (xrf-xri)/(n2)
  do i = 1,n2   ! During recombination
     xt(n1+i) = xri + i*xstep
  end do

  xstep = (x0-xrf)/(n3)
  do i = 1,n3   ! After recombination
     xt(n1+n2+i) = xrf + i*xstep
  end do

  at = exp(xt)

  ! compute conformal time at each time step
  allocate(xeta(n_eta))
  allocate(eta(n_eta))
  allocate(eta2(n_eta))

  ! uniform x grid
  xstep = (xeta2 - xeta1) / (n_eta-1)
  xeta(1) = xeta1
  do i = 2, n_eta
     xeta(i) = xeta(i-1) + xstep
  end do

 ! integrate eta
  xstep  = abs(1.d-2*(xeta(1)-xeta(2)))         ! xstep length
  eta(1) = eta1                                 ! initial value of et
  do i = 2, n_eta
     eta(i) = eta(i-1)
     call odeint(eta(i:i), xeta(i-1),xeta(i), eps, xstep, xstepmin, etaderivs, bsstep, output)
  end do

  ! write to file: eta, xeta
  open(1, file=fileplace//"eta.dat", action="write",status="replace")
  do i=1,n_eta
     write(1,*) eta(i), xeta(i)
  end do
  close(1)

  ! splining eta
  call spline(xeta, eta,yp1,ypn,eta2)

  ! write to file: spline
  open (2,file=fileplace//"etaspline.dat", action="write", status="replace")
  do i=1,nt
     write (2,*) get_eta(xt(i)), xt(i)
  end do
  close(2)

  ! calculating omegas
  open(3, file=fileplace//"omega_mb.dat", action="write", status="replace")
  open(5, file=fileplace//"omega_rv.dat", action="write", status="replace")
  do i=1, n_eta
    rho_crit = 3 * get_H(xeta(i))**2 / (8*pi*G_grav)

    rho_m = omega_m * rho_c * exp(xeta(i))**-3
    rho_b = omega_b * rho_c * exp(xeta(i))**-3
    rho_r = omega_r * rho_c * exp(xeta(i))**-4
    rho_v = omega_v * rho_c
    write(3,*) rho_m/rho_crit, rho_b/rho_crit
    write(5,*) rho_r/rho_crit, rho_v/rho_crit
  end do
  close(3)
  close(5)

  ! write to file: H(x), H(z)
  open(4, file=fileplace//"hubble.dat", action="write",status="replace")
  do i=1,n_eta
    z = 1-exp(-xeta(i))
    write(4,*) get_H(xeta(i)), z, get_H(log((1-z)**-1))
  end do
  close(4)
  end subroutine initialize_time_mod

  ! dnu/dx = c/H_p
  subroutine etaderivs(x,eta, derivative)
    use healpix_types
    implicit none
    real(dp), intent(in)                :: x
    real(dp), dimension(:), intent(in)  :: eta
    real(dp), dimension(:), intent(out) :: derivative
    derivative = c/get_Hp(x)
    end subroutine etaderivs

  ! find H at any time x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H
    get_H = H_0*sqrt((omega_b+omega_m)*exp(-3*x)+(omega_r+omega_nu)*exp(-4*x)+omega_v)
  end function get_H


  ! H' = aH : H prime
  function get_Hp(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_Hp
    get_Hp = H_0*sqrt((omega_b+omega_m)*exp(-x)+(omega_r+omega_nu)*exp(-2*x)+omega_v*exp(2*x))
  end function get_Hp



  ! dH'/dx : derivative of H prime
  function get_dHp(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dHp
    get_dHp = -H_0**2.d0*(0.5d0*(omega_b + omega_m)*exp(-x) + omega_r*exp(-2.d0*x) - omega_v*exp(2.d0*x))/(get_Hp(x))
  end function get_dHp



  ! find eta using splined values
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta
    get_eta =  splint(xeta, eta, eta2, x_in)
  end function get_eta


  subroutine output(x, y)
     implicit none
     real(dp),               intent(in)  :: x
     real(dp), dimension(:), intent(in)  :: y
  end subroutine output

end module time_mod
