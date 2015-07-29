! This program solves 1D Burger's equation at two scales
program tsburgers
  ! file_utils module takes care of input file and file i/o
  use file_utils, only: init_file_utils, finish_file_utils
  ! fftw library is used to perform 1D c2c,c2r,r2c Fourier Transform
  use fftw3

  implicit none

  ! basic parameters
  ! (kind=8) everywhere means values are in double precision
  ! constants remain unchanged in program
  real (kind=8), parameter :: pi=3.141592653589793
  complex (kind=8), parameter :: zi=(0.0,1.0)
  ! the problem is 1D, take the axis to be x axis
  ! lo_nx is total number of real-space grid points & mode numbers in low k(kappa)-space
  ! hi_nx is that in high k-space
  integer :: lo_nx, hi_nx
  ! lo_lx is box length for low k(kappa)-grid, hi_lx is that for high k-grid
  real (kind=8) :: lo_lx, hi_lx
  
  ! nstep is total number of time steps, nwrite is time step interval to write to file
  integer :: nstep, nwrite
  ! it is counter for time step
  integer :: it
  ! dt is time step size
  real (kind=8) :: dt
  ! time is  wall time
  real (kind=8) :: time = 0.0
  
  ! lo_visc is viscosity coefficient for low k(kappa), hi_visc is that for high k
  ! visc is to be used with wavenumbers divided by the one with maximum magnitude
  real (kind=8):: lo_visc, hi_visc

  ! hi_nk is number of useful non-negative wavenumbers evolved on low k(kappa)-grid
  ! lo_nk is that on high k-grid
  ! lo_nk_brk is location of the last wavenumber that is positive
  integer :: lo_nk, lo_nk_brk, hi_nk

  ! units for file i/o
  !integer :: hi_unit = 104!, lo_unit = 105

  ! xgrid contains location of grid points on low k(kappa)-grid
  real (kind=8), dimension(:), allocatable :: xgrid
  ! lo_kgrid contains location of grid points on low k(kappa)-grid
  ! hi_kgrid contains location of grid points on high k-grid
  real (kind=8), dimension(:), allocatable :: lo_kgrid, hi_kgrid, hi_kgrid_extended

  ! lo_fk is f(k=0,kappa) vector
  complex (kind=8), dimension(:), allocatable :: lo_fk, lo_fknew
  ! lo_fk_fft is lo_fk padded with zeros to avoid aliasing
  complex (kind=8), dimension(:), allocatable:: lo_fk_fft
  ! hi_fk is f(k,x) matrix
  complex (kind=8), dimension(:,:), allocatable :: hi_fk, hi_fknew
  ! f_k_kappa is the FT of hi_fk ( which is f(k,x) ) into kappa space
  complex (kind=8), dimension(:,:), allocatable :: f_k_kappa
  complex (kind=8), dimension(:), allocatable :: temp
  
  ! forcing
  ! lo_gamma is forcing parameter (effective growth rate) in low k(kappa) space
  ! hi_gamma is that in high k space
  complex (kind=8), dimension(:), allocatable :: lo_gamma, hi_gamma
  ! lo_gamma1 is value of forcing parameters in lo_gamma
  ! hi_gamma1 is that in hi_gamma
  real (kind=8) :: lo_gamma1, hi_gamma1
  ! lo_ikf is the wavenumber at which forcing of lo_gamma is applied
  ! hi_ifk is that of hi_gamma
  integer :: lo_ikf, hi_ikf

  ! nonlinearity
  ! hi_nk_fft is number of wavenumbers required as input to FFTW on high k grid
  ! lo_nk_fft = lo_nx
  integer :: hi_nk_fft

  ! diagnostics
  ! lo_total_energy is total energy in low k(kappa) lo_fk 
  ! hi_total_energy is that in high k hi_fk 
  real (kind=8) :: lo_total_energy, hi_total_energy
  ! units for file i/o writing diagnostics
  integer :: lo_spectrum_unit = 101, hi_spectrum_unit = 103, energy_unit = 102
  ! lo_spectrum vector is energy density spectrum in low k(kappa) lo_fk 
  ! hi_spectrum matrix is that in high k hi_fk 
  real (kind=8), dimension(:), allocatable :: lo_spectrum, hi_spectrum

  call init_file_utils
  call read_input_file
  call init_grids
  call init_f
  call init_diagnostics
  call init_forcing

  call write_diagnostics
  do it = 1, nstep
     call time_advance
     if (mod(it,nwrite)==0) call write_diagnostics
  end do

  call finish_diagnostics
  call finish_file_utils

  deallocate(xgrid, lo_kgrid, hi_kgrid, hi_kgrid_extended)
  deallocate(lo_fk, lo_fknew, lo_fk_fft)
  deallocate(hi_fk, hi_fknew)
  deallocate(f_k_kappa, temp)
  deallocate(lo_gamma, hi_gamma)
  deallocate(lo_spectrum, hi_spectrum)

contains

  subroutine read_input_file

    use file_utils, only: input_unit_exist, input_unit

    implicit none

    integer :: in_file
    logical :: exist
    ! lo_nx,hi_nx,nstep,nwrite,lo_ikf,hi_ikf are assumed to be integers
    ! lo_lx,hi_lx,dt,lo_visc,hi_visc,lo_gamma,hi_gamma are assumed to be real
    namelist / parameters / lo_nx, hi_nx, lo_lx, hi_lx, nstep, dt, nwrite,&
         lo_visc, hi_visc, lo_gamma1, hi_gamma1, lo_ikf, hi_ikf
    ! default values
    ! spatial grid
    lo_nx = 64
    hi_nx = 64
    lo_lx = 2.0*pi
    hi_lx = 2.0*pi/50.0
    ! time grid
    nstep = 100
    dt = 0.001
    nwrite = 1
    ! dissipation
    lo_visc = 40000.0
    hi_visc = 40000.0
    ! forcing
    lo_gamma1 = 1.0
    hi_gamma1 = 100.0
    lo_ikf = 2
    hi_ikf = 2
    
    in_file = input_unit_exist("parameters", exist)
    if (exist) read (unit=input_unit("parameters"), nml=parameters)

  end subroutine read_input_file

  subroutine init_grids
    
    implicit none

    integer :: ix, ik, ikappa

    ! IFT of hi_fk of each column is assumed to be real function, thus hi_fk(-k)=hi_fk(k)*
    ! half of hi_nx is required as input to IFT
    hi_nk_fft = hi_nx/2+1
    ! IFT of lo_fk is complex function, thus required number of input is lo_nx
    
    ! 2/3 rule is applied to avoid aliasing
    ! useful hi_nk is a third of hi_nx
    hi_nk = hi_nx/3+1
    ! useful lo_nk is a third of lo_nx both positive and negative
    lo_nk = 2*lo_nx/3+1
    ! as convention of fftw, output from x2k correspond to fk on this k-grid:
    ! 0, kf, 2*kf,...,(nx/2-1)*kf,-nx/2*kf,...,-kf
    ! lo_nk_brk is the position of the last useful wavenumber that is positive
    lo_nk_brk = lo_nx/3

    allocate(xgrid(lo_nx))
    allocate(lo_kgrid(lo_nk), hi_kgrid(hi_nk))

    ! x step size in low k(kappa) real space is lo_lx/lo_nx and it starts from 0
    do ix = 1, lo_nx
       xgrid(ix) = (lo_lx/lo_nx)*(ix-1)
    end do
    ! kappa step size is lo_kf = 2*pi/lo_lx
    ! kappa goes as 0, lo_kf, 2*lo_kf,...,(lo_nk_brk-1)*lo_kf, 
    ! (lo_nk_brk-lo_nk)*lo_kf,..., -lo_kf
    do ik = 1, lo_nk_brk
       lo_kgrid(ik) = (2.0*pi/lo_lx)*(ik-1)
    end do
    do ik = lo_nk_brk+1, lo_nk
       lo_kgrid(ik) = (2.0*pi/lo_lx)*(ik-lo_nk-1)
    end do
    ! k step size is hi_kf=2*pi/hi_lx and it goes as 0, hi_kf,..., hi_nx/3*hi_kf
    do ik = 1, hi_nk
       hi_kgrid(ik) = (2.0*pi/hi_lx)*(ik-1)
    end do

    allocate(hi_kgrid_extended(hi_nk*lo_nk))

    do ik = 1, hi_nk
       do ikappa = 1, lo_nk
          hi_kgrid_extended((ik-1)*lo_nk+ikappa) = &
               hi_kgrid(ik)+lo_kgrid(ikappa)
       end do
    end do
    
  end subroutine init_grids

  subroutine init_f
    
    implicit none

    real (kind=8) :: fkinit = 1.0
    integer :: ik

    allocate(lo_fk(lo_nk), lo_fknew(lo_nk)) ; lo_fk = 0.0 ; lo_fknew = 0.0
    allocate(lo_fk_fft(lo_nx)) ; lo_fk_fft = 0.0
    allocate(hi_fk(hi_nk,lo_nx), hi_fknew(hi_nk,lo_nx)) ; hi_fk = 0.0 ; hi_fknew = 0.0
    allocate(f_k_kappa(hi_nk, lo_nk), temp(lo_nx))

    ! lo_fk is initialized at the wavenumber at which forcing is applied
    ! lo_fk is initialized in this way just to make lo_fx real, this is not necessary
    lo_fk(lo_ikf) = fkinit
    lo_fk(lo_nk-lo_ikf+2) = fkinit
    
    ! lo_fk and hi_fk is set up in the way that IFT of lo_fk is equal to hi_fk(1,:)
    ! no rescaling here because fftw is assumed to follow the convention of magnifying
    ! at r2c and forward transformations but not at c2r or backward ones
    lo_fk_fft(:lo_nk_brk) = lo_fk(:lo_nk_brk)
    lo_fk_fft(lo_nx-lo_nk+lo_nk_brk+1:) = lo_fk(lo_nk_brk+1:)
    lo_fk_fft(lo_nk_brk+1:lo_nx-lo_nk+lo_nk_brk) = 0.0
    call lo_k2x_fft(lo_nx, lo_fk_fft, hi_fk(1,:))
    
    ! hi_fk is initialized at the wavenumber at which forcing is applied
    hi_fk(hi_ikf,:) = fkinit
    
    lo_fknew = lo_fk
    hi_fknew = hi_fk

    ! First row of f_k_kappa is for k=0, therefore equal to the array of lo_fknew
    f_k_kappa(1,:) = lo_fknew(:)
    ! From the second row on, each row of f_k_kappa is calculated by FT of 
    ! the corresponding row of hi_fknew
    do ik = 2, hi_nk
       call lo_x2k_fft(lo_nx, hi_fknew(ik,:), temp)
       f_k_kappa(ik,:lo_nk_brk) = temp(:lo_nk_brk)/lo_nx
       f_k_kappa(ik,lo_nk_brk+1:) = temp(lo_nk_brk+1+lo_nx-lo_nk:)/lo_nx
    end do

  end subroutine init_f

  subroutine init_diagnostics

    use file_utils, only: open_output_file

    implicit none
    
    integer :: istat

    allocate (lo_spectrum(lo_nk),hi_spectrum(hi_nk*lo_nk),STAT=istat)
    ! checker: if allocation is successful
    if (istat.NE.0) then
       print *, 'lo_spectrum or hi_spectrum is not allocated successfully'
    end if
    lo_spectrum = 0.0 ; hi_spectrum = 0.0 
    
    !open file and write title for each column
    call open_output_file (lo_spectrum_unit,'.lospectrum')
    write(lo_spectrum_unit,'(3a15)') '#time', 'lo_k', 'lo|f(k)|**2'

    call open_output_file (hi_spectrum_unit,'.hispectrum')
    write(hi_spectrum_unit,'(3a15)') '#time', 'hi_k', 'hi|f(k)|**2'

    call open_output_file (energy_unit,'.energy')
    write(energy_unit,'(3a15)') '#time', 'lo_total_energy', 'hi_total_energy'
    
  end subroutine init_diagnostics

  subroutine init_forcing

    implicit none

    ! allocate forcing arrays
    allocate(lo_gamma(lo_nk), hi_gamma(hi_nk))
    lo_gamma = 0.0 ; hi_gamma = 0.0

    ! forcing gamma of lo_gamma is applied at lo_k=lo_ikf
    ! *(1.0+zi) is to have linear growth rate and oscillation frequency of same magnitude
    lo_gamma(lo_ikf) = lo_gamma1*(1.0+zi)
    ! and also lo_k = -lo_ikf to make lo_fx real
    lo_gamma(lo_nk-lo_ikf+2) = lo_gamma1*(1.0+zi)

    ! forcing gamma of hi_gamma is applied at hi_k=hi_ikf
    hi_gamma(hi_ikf) = hi_gamma1*(1.0+zi)

  end subroutine init_forcing
  
  subroutine time_advance
    
    implicit none

    integer :: ix, ik, istat
    complex (kind=8), dimension(:), allocatable :: lo_nlk
    complex (kind=8), dimension(:,:), allocatable :: hi_nlk
    complex (kind=8), dimension(:), allocatable :: lo_interaction
    complex (kind=8), dimension(:,:), allocatable :: hi_interaction
    complex (kind=8) :: c

    ! allocate arrays lo_nlk, hi_nlk which are nonlinearity f*df/dx in k-space
    allocate (lo_nlk(lo_nx), hi_nlk(hi_nk_fft,lo_nx), STAT=istat)
    if (istat.NE.0) then
       print *, 'lo_nlk or hi_nlk is not allocated successfully'
    end if
    ! call subroutines to calculate nonlinearity
    call get_lo_nonlinearity (lo_fk, lo_nlk)
    call get_hi_nonlinearity (hi_fk, hi_nlk)
    
    ! allocate arrays lo_interaction and hi_interaction which are
    ! the interaction terms in the evolution of lo_fk and hi_fk
    allocate (lo_interaction(lo_nk), hi_interaction(hi_nk,lo_nx))
    lo_interaction = 0.0 ; hi_interaction = 0.0
    !call subroutines to calculate interaction
    call get_lo_interaction (lo_interaction)
    call get_hi_interaction (lo_fk, hi_fk, hi_interaction)

    ! centered scheme is used for viscosity and forcing terms
    ! (2-stage Runge-Kutta) Midpoint method is used for nonlinearity and interaction
    ! 0.5*dt is the middle stage

    ! update of lo_fknew is split into two do loops because
    ! lo_nlk array have different indexing than lo_fk and others
    do ik = 1, lo_nk_brk 
       c = 0.5*lo_visc*(lo_kgrid(ik)/maxval(abs(hi_kgrid)))**2-0.5*lo_gamma(ik)
       lo_fknew (ik) = (lo_fk(ik)*(1.0-0.5*dt*c)&
            -0.5*dt*lo_nlk(ik)-0.5*dt*lo_interaction(ik))&
            /(1.0+0.5*dt*c)
       ! checker: if lo_fknew is small set value to zero, if big stop and report
       if (abs(lo_fknew(ik)).LE.1.0D-10) then
          lo_fknew (ik) = 0.0
       end if
       if (abs(lo_fknew(ik)).GE.1.0D10) then
          print *, 'lo_fknew blows up at time =', time
          stop
       end if
    end do
    do ik = lo_nk_brk+1, lo_nk
       c = 0.5*lo_visc*(lo_kgrid(ik)/maxval(abs(hi_kgrid)))**2-0.5*lo_gamma(ik)
       lo_fknew (ik) = (lo_fk(ik)*(1.0-0.5*dt*c)&
            -0.5*dt*lo_nlk(ik-lo_nk+lo_nx)-0.5*dt*lo_interaction(ik))&
            /(1.0+0.5*dt*c)
       ! checker
       if (abs(lo_fknew(ik)).LE.1.0D-10) then
          lo_fknew (ik) = 0.0
       end if
       if (abs(lo_fknew(ik)).GE.1.0D10) then
          print *, 'lo_fknew blows up at time =', time
          stop
       end if
    end do

    ! update hi_fknew(1,:) with IFT of lo_fknew
    lo_fk_fft(:lo_nk_brk) = lo_fknew(:lo_nk_brk)
    lo_fk_fft(lo_nx-lo_nk+lo_nk_brk+1:) = lo_fknew(lo_nk_brk+1:)
    lo_fk_fft(lo_nk_brk+1:lo_nx-lo_nk+lo_nk_brk) = 0.0
    call lo_k2x_fft(lo_nx, lo_fk_fft, hi_fknew(1,:))

    ! the rest follows from equation of hi_fk
    do ik = 2, hi_nk
       do ix = 1, lo_nx
          c = 0.5*hi_visc*(hi_kgrid(ik)/maxval(abs(hi_kgrid)))**2-0.5*hi_gamma(ik)
          hi_fknew (ik,ix) = (hi_fk(ik,ix)*(1.0-0.5*dt*c)&
               -0.5*dt*hi_nlk(ik,ix)-0.5*dt*hi_interaction(ik,ix))/&
               (1.0+0.5*dt*c)
          ! checker
          if (abs(hi_fknew(ik,ix)).LE.1.0D-10) then
             hi_fknew(ik,ix) = 0.0
          end if
          if (abs(hi_fknew(ik,ix)).GE.1.0D10) then
             print *, 'hi_fknew blows up at time =', time
             stop
          end if
       end do
    end do

    ! time step advances with lo_fknew and hi_fknew at the middle stage
    ! call subroutines to calculate nonlinearities from lo_fknew and hi_fknew
    call get_lo_nonlinearity (lo_fknew, lo_nlk)
    call get_hi_nonlinearity (hi_fknew, hi_nlk)
    
    ! call subroutines to calculate interactions from lo_fknew and hi_fknew
    call get_lo_interaction (lo_interaction)
    call get_hi_interaction (lo_fknew, hi_fknew, hi_interaction)

    ! update with full time step
    do ik = 1, lo_nk_brk
       c = 0.5*lo_visc*(lo_kgrid(ik)/maxval(abs(hi_kgrid)))**2-0.5*lo_gamma(ik)
       lo_fknew(ik) = (lo_fk(ik)*(1.0-dt*c)&
            -dt*lo_nlk(ik)-dt*lo_interaction(ik))&
            /(1.0+dt*c)
       ! checker                     
       if (abs(lo_fknew(ik)).LE.1.0D-10) then
          lo_fknew(ik) = 0.0
       end if
       if (abs(lo_fknew(ik)).GE.1.0D10) then
          print *, 'lo_fknew blows up at time =', time
          stop
       end if
    end do
    do ik = lo_nk_brk+1, lo_nk
       c = 0.5*lo_visc*(lo_kgrid(ik)/maxval(abs(hi_kgrid)))**2-0.5*lo_gamma(ik)
       lo_fknew(ik) = (lo_fk(ik)*(1.0-dt*c)&
            -dt*lo_nlk(ik-lo_nk+lo_nx)-dt*lo_interaction(ik))&
            /(1.0+dt*c)
       ! checker                                          
       if (abs(lo_fknew(ik)).LE.1.0D-10) then
          lo_fknew(ik) = 0.0
       end if
       if (abs(lo_fknew(ik)).GE.1.0D10) then
          print *, 'lo_fknew blows up at time =', time
          stop
       end if
    end do
    
    ! update hi_fknew(1,:) with IFT of lo_fknew
    lo_fk_fft(:lo_nk_brk) = lo_fknew(:lo_nk_brk)
    lo_fk_fft(lo_nx-lo_nk+lo_nk_brk+1:) = lo_fknew(lo_nk_brk+1:)
    lo_fk_fft(lo_nk_brk+1:lo_nx-lo_nk+lo_nk_brk) = 0.0
    call lo_k2x_fft(lo_nx, lo_fk_fft, hi_fknew(1,:))

    ! the rest of hi_fknew is updated according to equation of hi_fk
    do ik = 2, hi_nk
       do ix = 1, lo_nx
          c = 0.5*hi_visc*(hi_kgrid(ik)/maxval(abs(hi_kgrid)))**2-0.5*hi_gamma(ik)
          hi_fknew (ik,ix) = (hi_fk(ik,ix)*(1.0-dt*c)&
               -dt*hi_nlk(ik,ix)-0.5*dt*hi_interaction(ik,ix))/&
               (1.0+dt*c)
          ! checker                  
          if (abs(hi_fknew(ik,ix)).LE.1.0D-10) then
             hi_fknew(ik,ix) = 0.0
          end if
          if (abs(hi_fknew(ik,ix)).GE.1.0D10) then
             print *, 'hi_fknew blows up at time =', time
             stop
          end if
       end do
    end do

    time = time + dt

    hi_fk = hi_fknew
    lo_fk = lo_fknew

    ! write result at this time step to array f_k_kappa
    f_k_kappa(1,:) = lo_fknew(:)
    do ik = 2, hi_nk
       call lo_x2k_fft (lo_nx, hi_fknew(ik,:), temp)
       f_k_kappa(ik, :lo_nk_brk) = temp(:lo_nk_brk)/lo_nx
       f_k_kappa(ik, lo_nk_brk+1:) = temp(lo_nk_brk+1+lo_nx-lo_nk:)/lo_nx
    end do

    deallocate(lo_nlk, hi_nlk, lo_interaction, hi_interaction)

  end subroutine time_advance

  subroutine get_lo_nonlinearity (lo_fk, lo_nlk)

    implicit none

    complex (kind=8), dimension(:), intent(in) :: lo_fk
    complex (kind=8), dimension(:), intent(out) :: lo_nlk
    complex (kind=8), dimension(:), allocatable :: lo_dfk_fft
    complex (kind=8), dimension(:), allocatable :: lo_fx, lo_dfdx
    complex (kind=8), dimension(:), allocatable :: lo_nlx
    integer :: istat

    allocate (lo_dfk_fft(lo_nx), STAT=istat)
    if (istat.NE.0) then
       print *, 'lo_dfk_fft is not allocated successfully'
    end if
    !padding lo_fk_fft with zeros
    lo_fk_fft(:lo_nk_brk) = lo_fk(:lo_nk_brk)
    lo_fk_fft(lo_nx-lo_nk+lo_nk_brk+1:) = lo_fk(lo_nk_brk+1:)
    lo_fk_fft(lo_nk_brk+1:lo_nx-lo_nk+lo_nk_brk) = 0.0
    !lo_dfk_fft is df/dx in k-space, also padded with zeros
    lo_dfk_fft(:lo_nk_brk) = zi*lo_kgrid(:lo_nk_brk)*lo_fk(:lo_nk_brk)
    lo_dfk_fft(lo_nx-lo_nk+lo_nk_brk+1:) = zi*lo_kgrid(lo_nk_brk+1:)*lo_fk(lo_nk_brk+1:)
    lo_dfk_fft(lo_nk_brk+1:lo_nx-lo_nk+lo_nk_brk) = 0.0
    
    allocate (lo_fx(lo_nx), lo_dfdx(lo_nx), STAT=istat)
    ! checker
    if (istat.NE.0) then
       print *,'lo_dfdx is not allocated successfully'
    end if

    ! inverse FFT and get fx, df/dx in real space
    call lo_k2x_fft (lo_nx, lo_fk_fft, lo_fx)
    call lo_k2x_fft (lo_nx, lo_dfk_fft, lo_dfdx)

    ! calculate nonlinear term f*(df/dx)in real space
    allocate (lo_nlx(lo_nx))
    lo_nlx = lo_fx*lo_dfdx

    ! forward FFT and get nonlinearity in low kappa space
    call lo_x2k_fft (lo_nx, lo_nlx, lo_nlk)
    ! normalize
    lo_nlk = lo_nlk/lo_nx

    deallocate(lo_dfk_fft, lo_fx, lo_dfdx, lo_nlx)
    
  end subroutine get_lo_nonlinearity

  subroutine lo_k2x_fft (n, fk, fx)

    use fftw3
    
    implicit none
    
    integer, intent(in) :: n
    complex (kind=8), dimension(:), intent(in) :: fk
    complex (kind=8), dimension(:), allocatable :: a
    complex (kind=8), dimension(:), intent(out) :: fx
    complex (kind=8), dimension(:), allocatable :: b
    integer (kind=8) :: plan
    
    allocate (a(n),b(n))
    ! create plan for 1d forward fft
    call dfftw_plan_dft_1d(plan, n, a, b, fftw_backward, fftw_estimate)
    a = fk
    ! execute 
    call dfftw_execute_dft(plan, a, b)
    call dfftw_destroy_plan(plan)
    fx = b

    deallocate (a, b)

  end subroutine lo_k2x_fft
  
  subroutine lo_x2k_fft (n, fx, fk)

    use fftw3

    implicit none

    integer, intent(in) :: n
    complex (kind=8), dimension(:), intent(in) :: fx
    complex (kind=8), dimension(:), allocatable :: b
    complex (kind=8), dimension(:), intent(out) :: fk
    complex (kind=8), dimension(:), allocatable :: a
    integer (kind=8) :: plan

    allocate (a(n),b(n))
    ! create plan for 1d backward fft
    call dfftw_plan_dft_1d(plan, n, b, a, fftw_forward, fftw_estimate)
    b = fx
    ! execute
    call dfftw_execute_dft(plan, b, a)
    call dfftw_destroy_plan(plan)
    fk = a
    
    deallocate (a, b)

  end subroutine lo_x2k_fft

  subroutine get_hi_nonlinearity (hi_fk, hi_nlk)

    implicit none

    complex (kind=8), dimension(:,:), intent(in) :: hi_fk
    complex (kind=8), dimension(:,:), intent(out) :: hi_nlk
    complex (kind=8), dimension(:,:), allocatable :: hi_fk_fft, hi_dfk_fft
    real (kind=8), dimension(:,:), allocatable :: hi_fx, hi_dfdx
    real (kind=8), dimension(:,:), allocatable :: hi_nlx
    integer :: istat, ix, ik
    
    allocate (hi_fk_fft(hi_nk_fft,lo_nx), hi_dfk_fft(hi_nk_fft,lo_nx), STAT=istat)
    ! checker
    if ((istat.NE.0).OR.(size(hi_dfk_fft).NE.lo_nx*hi_nk_fft).OR.&
         (size(hi_fk_fft).NE.lo_nx*hi_nk_fft)) then
       print *, 'hi_fk_fft or hi_dfk_fft is not allocated successfully'
    end if
    
    ! padding hi_fk_fft with zeros, calculating df/dx in k space
    do ik = 1, hi_nk
       do ix = 1, lo_nx
          hi_fk_fft(ik,ix) = hi_fk(ik,ix)
          hi_dfk_fft(ik,ix) = zi*hi_kgrid(ik)*hi_fk(ik,ix)
       end do
    end do
    hi_fk_fft(hi_nk+1:,:) = 0.0
    hi_dfk_fft(hi_nk+1:,:) = 0.0

    allocate (hi_fx(hi_nx,lo_nx), hi_dfdx(hi_nx,lo_nx), STAT=istat)
    ! checker
    if ((istat.NE.0).OR.(size(hi_dfdx).NE.hi_nx*lo_nx)) then
       print *, 'hi_dfdx is not allocated successfully'
    end if

    ! at each location of coarse grid ix, perform inverse FFT and 
    ! find f(x) and df/dx in real-space that correspond to high k-space
    do ix = 1, lo_nx
       call hi_k2x_fft (hi_nx, hi_fk_fft(:,ix), hi_fx(:,ix))
       call hi_k2x_fft (hi_nx, hi_dfk_fft(:,ix), hi_dfdx(:,ix))
    end do

    allocate(hi_nlx(hi_nx,lo_nx), STAT=istat)
    ! checker
    if ((istat.NE.0).OR.(size(hi_nlx).NE.hi_nx*lo_nx)) then
       print *, 'hi_nlx is not allocated successfully'
    end if
    
    ! at each location of coarse grid ix, calculate nonlinearity in real-space
    ! forward FT and get nonlinearity in high k-space
    do ix = 1, lo_nx
       hi_nlx(:,ix) = hi_fx(:,ix) * hi_dfdx(:,ix)
       call hi_x2k_fft (hi_nx, hi_nlx(:,ix), hi_nlk(:,ix))
    end do
    ! normalize
    hi_nlk(:hi_nk,:) = hi_nlk(:hi_nk,:)/hi_nx

    deallocate(hi_fk_fft, hi_dfk_fft, hi_fx, hi_dfdx, hi_nlx)

  end subroutine get_hi_nonlinearity

  subroutine hi_k2x_fft (n, fk, fx)
    
    use fftw3

    implicit none

    integer, intent(in) :: n
    integer :: n_fft
    complex (kind=8), dimension(:), intent(in) :: fk
    complex (kind=8), dimension(:), allocatable :: a
    real (kind=8), dimension(:), intent(out) :: fx
    real (kind=8), dimension(:), allocatable :: b
    integer (kind=8) :: plan
   
    n_fft = n/2+1
    allocate (b(n), a(n_fft))
    ! create plan for 1d complex-to-real FFT
    call dfftw_plan_dft_c2r_1d (plan, n, a, b, fftw_estimate)
    a = fk
    ! execute
    call dfftw_execute_dft_c2r (plan, a, b)
    call dfftw_destroy_plan (plan)
    fx = b
    
    deallocate (a, b)

  end subroutine hi_k2x_fft

  subroutine hi_x2k_fft (n, fx, fk)

    use fftw3

    implicit none

    integer, intent(in) :: n
    integer :: n_fft
    real (kind=8), dimension(:), intent(in) :: fx
    real (kind=8), dimension(:), allocatable :: b
    complex (kind=8), dimension(:), intent(out) :: fk
    complex (kind=8), dimension(:), allocatable :: a
    integer (kind=8) :: plan

    n_fft = n/2+1
    allocate (b(n), a(n_fft))
    ! create plan for 1d real-to-complex FFT
    call dfftw_plan_dft_r2c_1d (plan, n, b, a, fftw_estimate)
    b = fx
    ! execute
    call dfftw_execute_dft_r2c (plan, b, a)
    call dfftw_destroy_plan (plan)
    fk = a

    deallocate (a, b)

  end subroutine hi_x2k_fft

  subroutine get_lo_interaction (lo_interaction)

    implicit none

    complex (kind=8), dimension(:), intent(out) :: lo_interaction
    integer :: ik, ikappa, ikappa1, ikappa2

    lo_interaction = 0.0
    ! ikappa is index for target low wavenumber, k is index for target high wavenumber
    ! ikappa1 and ikappa2 are indices for source wavenumbers
    ! start with ikappa = 2 because ikappa = 1 has no contribution 
    do ikappa = 2, lo_nk 
       do ik = 1, hi_nk
          do ikappa1 = 2, lo_nk
             do ikappa2 = 2, lo_nk
                if (lo_kgrid(ikappa1) == lo_kgrid(ikappa)) then
                   lo_interaction (ikappa) = lo_interaction (ikappa) + &
                        zi*lo_kgrid(ikappa)*f_k_kappa(ik, ikappa1) * &
                        conjg(f_k_kappa(ik, 1))
                else if (lo_kgrid(ikappa2) == lo_kgrid(ikappa)) then
                   lo_interaction (ikappa) = lo_interaction (ikappa) + &
                        zi*lo_kgrid(ikappa)*f_k_kappa(ik, 1) * &
                        conjg(f_k_kappa(ik, lo_nk-ikappa2+2))
                   if (lo_kgrid(ikappa1) + lo_kgrid(ikappa2) == lo_kgrid(ikappa)) then
                      lo_interaction (ikappa) = lo_interaction (ikappa) + &
                           zi*lo_kgrid(ikappa)*f_k_kappa(ik, ikappa1) * &
                           conjg(f_k_kappa(ik, lo_nk-ikappa2+2))
                   end if
                end if
             end do
          end do
       end do
    end do
    
  end subroutine get_lo_interaction

  subroutine get_hi_interaction (lo_fk, hi_fk, hi_interaction)

    implicit none

    complex (kind=8), dimension(:), intent(in) :: lo_fk
    complex (kind=8), dimension(:,:), intent(in) :: hi_fk
    complex (kind=8), dimension(:,:), intent(out) :: hi_interaction
    integer :: ix, ik, ikappa

    hi_interaction = 0.0
    do ik = 1, hi_nk
       do ix = 1, lo_nx
          do ikappa = 1, lo_nk
             hi_interaction(ik,ix) = hi_interaction(ik,ix) &
                  +hi_fk(ik,ix)*zi*(hi_kgrid(ik)+lo_kgrid(ikappa))*lo_fk(ikappa) &
                  *exp(zi*lo_kgrid(ikappa)*xgrid(ix))
          end do
       end do
    end do

  end subroutine get_hi_interaction

  subroutine write_diagnostics 
    
    implicit none
    
    integer :: ik, ikappa

    ! subroutine to calculate spectrum and total energy
    call get_diagnostics 
    ! output lo_spectrum
    do ik = 1, lo_nk
       write(lo_spectrum_unit, '(3e15.6)') time, lo_kgrid(ik), lo_spectrum(ik)
    end do
    write(lo_spectrum_unit,*)

    ! output hi_spectrum
    do ikappa = 1, lo_nk_brk
       write(hi_spectrum_unit,'(3e15.6)') time, hi_kgrid_extended(ikappa), hi_spectrum(ikappa)
    end do
    write(hi_spectrum_unit,*)
    do ik = 2, hi_nk
       do ikappa = lo_nk_brk+1, lo_nk
          write(hi_spectrum_unit,'(3e15.6)') time, hi_kgrid_extended((ik-1)*lo_nk+ikappa),&
               hi_spectrum((ik-1)*lo_nk+ikappa)
       end do
       write(hi_spectrum_unit,*)
       do ikappa = 1, lo_nk_brk
          write(hi_spectrum_unit,'(3e15.6)') time, hi_kgrid_extended((ik-1)*lo_nk+ikappa),&
               hi_spectrum((ik-1)*lo_nk+ikappa)
       end do
       write(hi_spectrum_unit,*)
    end do

    ! output energy
    write (energy_unit,'(3e15.6)') time, lo_total_energy, hi_total_energy

  end subroutine write_diagnostics

  subroutine get_diagnostics

    implicit none

    integer :: ik, ikappa
    ! fkkappa_density is energy density of f_k_kappa
    real(kind=8), dimension(:,:), allocatable :: fkkappa_density
    
    ! spectrum = 1/2*|f|**2
    lo_spectrum = 0.5*real(lo_fknew*conjg(lo_fknew))
    ! energy = sum of energy density in the spectrum
    lo_total_energy = sum(lo_spectrum)

    allocate (fkkappa_density(hi_nk, lo_nk))
    fkkappa_density = 0.5*real(f_k_kappa*conjg(f_k_kappa))

    ! hi_spectrum is fkkappa_density reorganized in a linear vector
    ! hi_kgrid step size = int(hi_kgrid(2))
    do ik = 1, hi_nk
       do ikappa = 1, lo_nk
          hi_spectrum((ik-1)*lo_nk+ikappa) = fkkappa_density(ik,ikappa)
       end do
    end do
    hi_total_energy = sum(hi_spectrum)

  end subroutine get_diagnostics

  subroutine finish_diagnostics

    use file_utils, only: close_output_file

    implicit none

    call close_output_file (hi_spectrum_unit)
    call close_output_file (lo_spectrum_unit)
    call close_output_file (energy_unit)

  end subroutine finish_diagnostics

end program tsburgers
