! This program solves 1D Burger's equation
! with small space step and big space range
 
program fullburgers
  
  use file_utils, only: init_file_utils, finish_file_utils
  
!  use fftw3
  
  include 'fftw3.f03'
  
  implicit none
  
  ! (kind=8) means value is defined with double precision
  
  ! constants whose values remain unchanged in program
  real (kind=8), parameter :: pi=3.141592653589793
  complex (kind=8), parameter :: zi=(0.0,1.0)
  
  ! nx is the number of real-space grid points/modes in x-direction
  integer :: nx
  ! Lx is the box length in x-direction
  real (kind=8) :: Lx

  ! nstep is the number of time steps
  ! nwrite is time step interval to write to file
  ! ntransfer is time step interval to evaluate transfer function
  integer :: nstep, nwrite, ntransfer
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
  integer :: spectrum_unit=101, energy_unit=103, transfer_unit=102, trans_fn_unit=104

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
  ! transfer is matrix of transfer function 
  ! transfer_x is transfer function for positive k from both positive and negative k1
  ! sub_sum has two columns, the first is the sum of transfer_x whose ik is smaller than k
  ! the second column is the sum of transfer_x whose ik is bigger than k
  ! part_sum has two columns, the first is the sum of transfer_x whose ik is smaller than ikf2
  ! the second column is the sum of transfer_x whose ik is bigger than ikf2
  real (kind=8), dimension(:,:), allocatable :: transfer, transfer_avg, transfer_tot
  real (kind=8), dimension(:,:), allocatable :: transfer_x, transf_x_avg, transf_x_tot
  real (kind=8), dimension(:,:), allocatable :: sub_sum, part_sum

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

  if (allocated(spectrum_tot)) then
     deallocate (spectrum_tot)
     deallocate (spectrum_avg)
!     deallocate (transfer_tot)
!     deallocate (transfer_avg)
     deallocate (transf_x_tot)
     deallocate (transf_x_avg)
     deallocate (sub_sum)
     deallocate (part_sum)
  end if

contains

  subroutine read_input_file

    use file_utils, only: input_unit_exist, input_unit

    implicit none

    integer :: in_file
    logical :: exist

    namelist / parameters / nx, Lx, nstep, dt, nwrite, ntransfer, &
         gamma1, gamma2, ikf1, ikf2, visc

    ! default values

    ! spatial grid
    nx = 4096
    Lx = 2.0*pi
    ! time grid
    nstep = 200
    dt = 0.001
    nwrite = 20
    ntransfer = 50
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
       xgrid(ix) = (Lx*(ix-1))/nx
    end do

    do ik = 1, nk
       kgrid(ik) = (2.0*pi/Lx)*(ik-1)
    end do

    do ik = 1, (2*nk-1)
       kgrid_x(ik) = (2.0*pi/Lx)*(ik-nk)
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
!    fk(ikf2:ikf2+10) = fkinit*.3

    fknew = fk

  end subroutine init_f

  subroutine init_diagnostics

    use file_utils, only: open_output_file

    implicit none

    call open_output_file (spectrum_unit,'.spectrum')
    write (spectrum_unit,'(4a12)') '# time', 'k', '|f(k)|**2', 'specavg'

    call open_output_file (energy_unit,'.energy')
    write (energy_unit,'(3a12)') '# time', 'energy', 'energyavg'

    call open_output_file (transfer_unit,'.transfer')
!    write (transfer_unit,'(6a12)') '# time', 'k1', 'k2', 'trans_fn', &
!         'maxval', 'tavg'
    write(transfer_unit,'(6a12)') '# time', 'k', 'loc_sum', 'nloc_sum', 'low_sum', 'high_sum'
    
    call open_output_file (trans_fn_unit,'.transfn')
    write (trans_fn_unit,'(7a12)') '# time', 'k1', 'k', 'trans_fn', 'tavg'

    allocate (spectrum(nk)) ; spectrum = 0.0
    allocate (transfer(nk,nk)) ; transfer = 0.0
    allocate (transfer_x(2*nk-1,nk)) ; transfer_x = 0.0
  end subroutine init_diagnostics

  subroutine init_forcing

    implicit none
    
    allocate (gamma(nk)) ; gamma = 0.0
    gamma(ikf1) = gamma1*(1.0+zi)
!    gamma(ikf2-5:ikf2+5) = gamma2*exp(-(kgrid(ikf2)-kgrid(ikf2-5:ikf2+5))**2/(kgrid(ikf2)-kgrid(ikf2-5))**2)
!    gamma(ikf2:ikf2+10) = gamma2*(1.0 + (kgrid(ikf2:ikf2+10)-kgrid(ikf2))/kgrid(ikf2))
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

!    use fftw3

    include 'fftw3.f03'

    implicit none

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

!    use fftw3

    include 'fftw3.f03'

    implicit none

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

    if (ntransfer < step .and.  mod(step,ntransfer)==0 .and. step>0) then
       do ik1 = 1, 2*nk-1
          do ik= 1, nk
             write (trans_fn_unit,'(7e12.4)') time, kgrid_x(ik1), kgrid(ik), &
                  transfer_x(ik1,ik)/total_energy, transf_x_avg(ik1,ik)!, &
!                  sub_sum(1,ik)/total_energy, sub_sum(2,ik)/total_energy
          end do
          write(trans_fn_unit,*)
       end do
       write (trans_fn_unit,*)
       write (trans_fn_unit,*)
    end if

    !call get_energy_transfer (step)

    write (energy_unit,'(3e12.4)') time, total_energy, energy_avg

    do ik = 1, nk
       write (spectrum_unit,'(4e12.4)') time, kgrid(ik), spectrum(ik), spectrum_avg(ik)
    end do
    write (spectrum_unit,*)

    tmax = maxval(transfer)/total_energy
    tmin = minval(transfer)/total_energy

    if (ntransfer < step .and. mod(step,ntransfer)==0 .and. step>0) then
!       do ik1 = 1, nk
          do ik = 1, nk
!             write (transfer_unit,'(6e12.4)') time, kgrid(ik1), kgrid(ik), &
!                  transfer(ik1,ik)/total_energy, max(tmax,abs(tmin)), transfer_avg(ik1,ik)
             write (transfer_unit,'(6e12.4)') time, kgrid(ik), sub_sum(1,ik)/total_energy, &
                  sub_sum(2,ik)/total_energy, part_sum(1,ik)/total_energy, part_sum(2,ik)/total_energy
          end do
          write(transfer_unit,*)
!       end do
!       write(transfer_unit,*)
!       write(transfer_unit,*)
    end if

  end subroutine write_diagnostics

  subroutine get_energy_transfer_new (step)

    implicit none

    integer, intent (in) :: step

    integer, save :: navg = 0
    real (kind=8), save :: energy_tot = 0.0
    integer :: ik1, ik

    if (.not.allocated(transf_x_tot)) then
       allocate (spectrum_tot(nk)) ; spectrum_tot = 0.0
       allocate (spectrum_avg(nk)) ; spectrum_avg = 0.0
       allocate (transf_x_tot(2*nk-1,nk)) ; transf_x_tot = 0.0
       allocate (transf_x_avg(2*nk-1,nk)) ; transf_x_avg = 0.0
       allocate (sub_sum(2,nk)) ; sub_sum = 0.0
       allocate (part_sum(2,nk)) ; part_sum = 0.0
    end if

    total_energy = sum(real(fknew*conjg(fknew)))*0.5
    spectrum = 0.5*real(fknew*conjg(fknew))

    sub_sum = 0.0
    transfer_x = 0.0
    do ik = 1, nk
       do ik1 = nk, 2*nk-1
          if (ik1-nk-ik < 0) then
             transfer_x(ik1, ik) = 0.5*kgrid(ik)*aimag(fknew(ik1-nk+1)*fknew(ik-ik1+nk)*conjg(fknew(ik)))
             sub_sum(1,ik) = sub_sum(1,ik) + transfer_x(ik1,ik)
          else
             transfer_x(ik1, ik) = 0.5*kgrid(ik)*aimag(fknew(ik1-nk+1)*conjg(fknew(ik1-nk-ik+2))*conjg(fknew(ik)))
             sub_sum(2,ik) = sub_sum(2,ik) + transfer_x(ik1,ik)
          end if
       end do
    end do
    
    do ik = 1, nk
       do ik1 = 1, nk-1
          if (ik1-ik-1 < 0) then
             transfer_x(ik1, ik) = 0
          else
             transfer_x(ik1, ik) = 0.5*kgrid(ik)*aimag(conjg(fknew(nk-ik1+1))*fknew(ik+nk-ik1)*conjg(fknew(ik)))
             !transfer_x(ik1,ik) = transfer_x(2*nk+ik-ik1-1,ik)
             sub_sum(2,ik) = sub_sum(2,ik) + transfer_x(ik1,ik)
          end if
       end do
    end do

    part_sum = 0.0
    do ik = 1, nk
       do ik1 = nk, 2*nk-1
          if (kgrid_x(ik1)<kgrid(11)) then
             part_sum(1,ik) = part_sum(1,ik) + transfer_x(ik1,ik)
          else
             part_sum(2,ik) = part_sum(2,ik) + transfer_x(ik1,ik)
          end if
       end do
    end do
    
    if (step > nstep/5) then
       energy_tot = energy_tot + total_energy
       spectrum_tot = spectrum_tot + spectrum
       transf_x_tot = transf_x_tot + transfer_x/total_energy

       navg = navg + 1

       energy_avg = energy_tot/navg
       spectrum_avg = spectrum_tot/navg
       transf_x_avg = transf_x_tot/navg

    end if

  end subroutine get_energy_transfer_new

   subroutine get_energy_transfer (step)

    implicit none

    integer, intent (in) :: step

    integer :: ik1, ik2
    integer, save :: navg = 0
    real (kind=8), save :: energy_tot = 0.0

 if (.not.allocated(transfer_tot)) then
!       allocate (spectrum_tot(nk)) ; spectrum_tot = 0.0
!       allocate (spectrum_avg(nk)) ; spectrum_avg = 0.0
       allocate (transfer_tot(nk,nk)) ; transfer_tot = 0.0
       allocate (transfer_avg(nk,nk)) ; transfer_avg = 0.0
    end if

    total_energy = sum(real(fknew*conjg(fknew)))*0.5
    spectrum = 0.5*real(fknew*conjg(fknew))

    transfer = 0.0
    do ik1 = 1, nk
       do ik2 = 1, nk
          if (ik1-ik2 >= 0) then
             transfer(ik1,ik2) = 0.5*kgrid(ik2)*aimag(fknew(ik2)*conjg(fknew(ik1))*fknew(ik1-ik2+1))
          else if (nk > 1) then
             transfer(ik1,ik2) = 0.5*kgrid(ik2)*aimag(fknew(ik2)*conjg(fknew(ik1))*conjg(fknew(ik2-ik1+1)))
          end if
       end do
    end do

    if (step > nstep/5) then
       energy_tot = energy_tot + total_energy
       spectrum_tot = spectrum_tot + spectrum
       transfer_tot = transfer_tot + transfer/total_energy

       navg = navg + 1

       energy_avg = energy_tot/navg
       spectrum_avg = spectrum_tot/navg
       transfer_avg = transfer_tot/navg

    end if

  end subroutine get_energy_transfer

  subroutine finish_diagnostics

    use file_utils, only: close_output_file

    implicit none

    call close_output_file (spectrum_unit)
    call close_output_file (energy_unit)
    call close_output_file (transfer_unit)
    call close_output_file (trans_fn_unit)

    deallocate (transfer)
    deallocate (transfer_x)

  end subroutine finish_diagnostics

end program fullburgers
