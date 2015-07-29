! This program solves 1D Burger's equation
! with small space step and big space range
 
program fullburgers
  
  use file_utils, only: init_file_utils, finish_file_utils
  
  implicit none
  include 'fftw3.f'
  
  
  ! (kind=8) means value is defined with double precision
  
  ! constants whose values remain unchanged in program
  real (kind=8), parameter :: pi=3.141592653589793
  complex (kind=8), parameter :: zi=(0.0,1.0)
  
  ! nx is the number of real-space grid points/modes in x-direction
  integer :: nx
  ! lx is the box length in x-direction
  real (kind=8) :: lx

  ! nstep is the number of time steps
  ! nwrite is time step interval to write to file
  integer :: nstep, nwrite
  ! dt is the time step size
  real :: dt

  ! gamma1 and gamma2 are the forcing parameters (equivalent growth rate) for k1 and k2
  real :: gamma1, gamma2
  ! ikf1 and ikf2 are indices at which to force k1 and k2
  integer :: ikf1, ikf2

  ! visc is the viscosity coefficient
  real :: visc

  ! nk_fft is the number of non-negative wavenumbers used in FFTs to avoid aliasing
  integer :: nk_fft
  ! nk is the number of non-negative wavenumbers actually evolved
  integer :: nk

  integer :: it
  real :: time = 0.0
  ! defining some units for file i/o
  integer :: spectrum_unit=101, energy_unit=103

  ! xgrid contains location of grid points in x-direction
  real (kind=8), dimension (:), allocatable :: xgrid
  ! kgrid contains location of grid points in k-space
  ! kgrid_x contains negative k values 
  real (kind=8), dimension (:), allocatable :: kgrid, kgrid_x
  ! gamma is the forcing parameter (effective growth rate spectrum)
  complex (kind=8), dimension (:), allocatable :: gamma
  ! fk is what is solved for, fknew is value of fk at next time step
  complex(kind=8), dimension(:), allocatable :: fk, fknew
  ! spectrum is energy density at each wavevector k
  real (kind=8), dimension (:), allocatable :: spectrum, spectrum_avg, spectrum_tot
  ! total_energy is total energy at each time step
  real (kind=8) :: total_energy, energy_avg=0.0

  call init_file_utils
  call read_input_file
  call init_grids
  call init_f
  call init_diagnostics
  call init_forcing

  call write_diagnostics (step=0)
  do it = 1, nstep
     call time_advance
     if (mod(it,nwrite)==0) call write_diagnostics (it)
  end do

  call finish_diagnostics
  call finish_file_utils

  deallocate (xgrid, kgrid, kgrid_x)
  deallocate (fk, fknew)
  deallocate (gamma)
  deallocate (spectrum)

contains

  subroutine read_input_file

    use file_utils, only: input_unit_exist, input_unit

    implicit none

    integer :: in_file
    logical :: exist

    namelist / parameters / nx, lx, nstep, dt, nwrite, &
         gamma1, gamma2, ikf1, ikf2, visc

    ! default values

    ! spatial grid
    nx = 4096
    lx = 2.0*pi
    ! time grid
    nstep = 200
    dt = 0.001
    nwrite = 20
    !ntransfer = 50
    ! forcing
    gamma1 = 2.0
    gamma2 = 10.0
    ikf1 = 2
    ikf2 = 101
    ! dissipation
    visc = 40000

    in_file = input_unit_exist ("parameters", exist)
    if (exist) read (unit=input_unit("parameters"), nml=parameters)

  end subroutine read_input_file

  subroutine init_grids
    
    implicit none
    
    integer :: ix, ik

    nk_fft = nx/2+1
    nk = nx/3+1
    
    allocate (xgrid(nx))
    allocate (kgrid(nk), kgrid_x(2*nk-1))

    do ix = 1, nx
       xgrid(ix) = (lx*(ix-1))/nx
    end do

    do ik = 1, nk
       kgrid(ik) = (2.0*pi/lx)*(ik-1)
    end do

    do ik = 1, (2*nk-1)
       kgrid_x(ik) = (2.0*pi/lx)*(ik-nk)
    end do

  end subroutine init_grids

  subroutine init_f

    implicit none

    real :: fkinit = 1.0

    ! allocate k-space arrays
    allocate (fk(nk), fknew(nk))

    ! initialize padded arrays with zeros
    fk = 0.

    ! initialize a single k with amplitude fkinit
    fk(ikf1) = fkinit
    fk(ikf2) = fkinit

    fknew = fk

  end subroutine init_f

  subroutine init_diagnostics

    use file_utils, only: open_output_file

    implicit none

    call open_output_file (spectrum_unit,'.spectrum')
    write (spectrum_unit,'(4a12)') '# time', 'k', '|f(k)|**2', 'specavg'

    call open_output_file (energy_unit,'.energy')
    write (energy_unit,'(3a12)') '# time', 'energy', 'energyavg'

    allocate (spectrum(nk)) ; spectrum = 0.0

  end subroutine init_diagnostics

  subroutine init_forcing

    implicit none
    
    allocate (gamma(nk)) ; gamma = 0.0
    gamma(ikf1) = gamma1*(1.0+zi)
    gamma(ikf2) = gamma2*(1.0+zi)

  end subroutine init_forcing

  subroutine time_advance

    implicit none
    
    complex (kind=8), dimension (:), allocatable :: nlk
    integer :: ik
    allocate (nlk(nk_fft)) ; nlk = 0.0

    call get_nonlinearity (fk, nlk)

    fknew = (fk-0.5*dt*nlk(:nk)) / (1.0 + 0.5*dt &
         *(-gamma + visc*(kgrid/kgrid(nk))**2))
    ! checker
    do ik = 1, nk
       if (abs(fknew(ik)).LE.1.0D-10) then
          fknew(ik) = 0.0
       end if
       if (abs(fknew(ik)).GE.1.0D10) then
          print*,'fknew blows up at time =', time
          stop
       end if
    end do

    call get_nonlinearity (fknew, nlk)

    ! note that hypervisc is normalized by kmax, but visc is not
    fknew = (fk-dt*nlk(:nk)) / (1.0 + dt &
         *(-gamma + visc*(kgrid/kgrid(nk))**2))

    ! checker
    do ik = 1, nk
       if (abs(fknew(ik)).LE.1.0D-10) then
          fknew(ik) = 0.0
       end if
       if (abs(fknew(ik)).GE.1.0D10) then
          print*,'fknew blows up at time =', time
          stop
       end if
    end do

    ! update time
    time = time + dt
    fk = fknew

    deallocate (nlk)

  end subroutine time_advance

  subroutine get_nonlinearity (fk, nlk)
    
    implicit none
    
    complex (kind=8), dimension (:), intent (in) :: fk
    complex (kind=8), dimension (:), intent (out) :: nlk
    real(kind=8), dimension(:), allocatable :: fx, dfdx
    complex(kind=8), dimension(:), allocatable :: fk_fft, dfk_fft    
    real (kind=8), dimension (:), allocatable :: nlx
    
    allocate(fk_fft(nk_fft), dfk_fft(nk_fft))
    ! pad f with zeros for dealiasing
    fk_fft(:nk) = fk ; fk_fft(nk+1:) = 0.
    ! this is df/dx in k-space
    dfk_fft(:nk) = zi*kgrid*fk ; dfk_fft(nk+1:) = 0.
    
    allocate(fx(nx), dfdx(nx))
    ! calculate nonlinearity explicitly using pseudospectral method
    call k2x_fft (fk_fft, fx)
    call k2x_fft (dfk_fft, dfdx)
    
    allocate (nlx(nx))
    
    nlx = fx*dfdx
    
    call x2k_fft (nlx, nlk)
    
    ! take care of fft normalization
    nlk(:nk) = nlk(:nk)/nx
    
    deallocate (nlx)
    deallocate (fx, dfdx)
    deallocate (fk_fft, dfk_fft)
    
  end subroutine get_nonlinearity

  subroutine k2x_fft (fk, fx)

    implicit none
    include 'fftw3.f'

    complex(kind=8), dimension(:), intent(in) :: fk
    complex(kind=8), dimension(:), allocatable :: a
    real(kind=8), dimension(:), intent(out) :: fx
    real(kind=8), dimension(:), allocatable :: b
    integer(kind=8) :: plan

    allocate (a(nk_fft),b(nx))
    call dfftw_plan_dft_c2r_1d (plan, nx, a, b, fftw_estimate)
    a = fk
    call dfftw_execute_dft_c2r (plan, a, b)
    call dfftw_destroy_plan (plan)
    fx = b

    deallocate(a, b)

  end subroutine k2x_fft

  subroutine x2k_fft (fx, fk)

    implicit none
    include 'fftw3.f'

    real(kind=8), dimension(:), intent(in) :: fx
    real(kind=8), dimension(:), allocatable :: b
    complex(kind=8), dimension(:), intent(out) :: fk
    complex(kind=8), dimension(:), allocatable :: a
    integer(kind=8) :: plan

    allocate (a(nk_fft),b(nx))
    call dfftw_plan_dft_r2c_1d (plan, nx, b, a, fftw_estimate)
    b = fx
    call dfftw_execute_dft_r2c (plan, b, a)
    call dfftw_destroy_plan (plan)
    fk = a

    deallocate (a, b)

  end subroutine x2k_fft

  subroutine write_diagnostics (step)

    implicit none

    integer, intent (in) :: step

    integer :: ik, ik1
    real (kind=8) :: tmin, tmax

    call get_energy_transfer_new (step)

    write (energy_unit,'(3e12.4)') time, total_energy, energy_avg

    do ik = 1, nk
       write (spectrum_unit,'(4e12.4)') time, kgrid(ik), spectrum(ik), spectrum_avg(ik)
    end do
    write (spectrum_unit,*)

  end subroutine write_diagnostics

  subroutine get_energy_transfer_new (step)

    implicit none

    integer, intent (in) :: step

    integer, save :: navg = 0
    real (kind=8), save :: energy_tot = 0.0
    integer :: ik1, ik

    if (.not.allocated(spectrum_tot)) then
       allocate (spectrum_tot(nk)) ; spectrum_tot = 0.0
       allocate (spectrum_avg(nk)) ; spectrum_avg = 0.0
    end if

    total_energy = sum(real(fknew*conjg(fknew)))*0.5
    spectrum = 0.5*real(fknew*conjg(fknew))
    
    if (step > nstep/5) then
       energy_tot = energy_tot + total_energy
       spectrum_tot = spectrum_tot + spectrum

       navg = navg + 1

       energy_avg = energy_tot/navg
       spectrum_avg = spectrum_tot/navg

    end if

  end subroutine get_energy_transfer_new

   subroutine get_energy_transfer (step)

    implicit none

    integer, intent (in) :: step

    integer :: ik1, ik2
    integer, save :: navg = 0
    real (kind=8), save :: energy_tot = 0.0

 if (.not.allocated(spectrum_tot)) then
       allocate (spectrum_tot(nk)) ; spectrum_tot = 0.0
       allocate (spectrum_avg(nk)) ; spectrum_avg = 0.0
    end if

    total_energy = sum(real(fknew*conjg(fknew)))*0.5
    spectrum = 0.5*real(fknew*conjg(fknew))

    if (step > nstep/5) then
       energy_tot = energy_tot + total_energy
       spectrum_tot = spectrum_tot + spectrum

       navg = navg + 1

       energy_avg = energy_tot/navg
       spectrum_avg = spectrum_tot/navg

    end if

  end subroutine get_energy_transfer

  subroutine finish_diagnostics

    use file_utils, only: close_output_file

    implicit none

    call close_output_file (spectrum_unit)
    call close_output_file (energy_unit)

  end subroutine finish_diagnostics

end program fullburgers
