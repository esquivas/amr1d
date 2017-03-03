!=======================================================================
!   This program solves the Euler equations with the  Godunov Method
!   using the HLLC Fluxes (a modification of HLLC proposed by Toro, to
!   account for the "C"ontact discontinuity, all the difference from HLL
!   is in the fluxes (prim2hllc) subroutine, with a block based Adaptive
!   Mesh
!=======================================================================

!=======================================================================
!   This module contains global variables
module globals
  implicit none
  !
  !   This is the number of points used to discretize X
  integer, parameter :: nx=50
  !   Number of levels allowed
  integer, parameter :: nlevs = 5
  !   maximumn number of blocks
  integer, parameter :: nbmax = 16
  !   Here we set the extent of X
  real, parameter :: xmax=1.
  real, parameter :: gamma=1.4
  !   This is a vector that contains u(x)
  real :: u(3,0:nx+1,nbmax), prim(3,0:nx+1,nbmax)
  !  other arrays and variables needed in AMR
  real    :: dx(nlevs)
  integer :: minID(nlevs), maxID(nlevs),ActiveBlocks(nbmax)
  integer :: lastActive
end module globals

!=======================================================================
!   main program
program euler_amr
  use globals
  implicit none
  ! declaration of some variables needed by the main program
  real            :: time, dt= 0.         !  t, $\Delta t$
  real, parameter :: tmax= .3             ! maximumn integration time
  real, parameter :: dtprint=0.01         ! interval between outputs
  real            :: tprint               ! time of next output
  integer         :: itprint              ! number of current output

  !  Initializes AMR
  call init_mesh()

  ! This subroutine generates the initial conditions
  call initconds(time, tprint, itprint)

  !   main loop, iterate until maximum time is reached
  do while (time.lt.tmax)

    ! updates the primitives
    call u2prim(u,prim)

    ! output at tprint intervals
    if(time >= tprint) then
      write(*,*) time,tmax,dt
      call output(itprint)
      tprint=tprint+dtprint
      itprint=itprint+1
    end if

    ! Obtain the $\Delta t$ allowed by the CFL criterium
    call timestep(dt)
    !
    ! Integrate u fom t to t+dt
    call tstep(dt,time)
    ! time counter increases
    time=time+dt

  end do

  stop
end program euler_amr

!=======================================================================
!  Initializes all things mesh
subroutine init_mesh()
  use globals, only : nx, nlevs, xmax, dx, minID, maxID, ActiveBlocks, &
                      lastActive
  implicit none
  integer :: nl

  !  compute dx
  dx(1) = xmax/real(nx)
  do nl=2,nlevs
    dx(nl)= dx(nl-1)/2.
  end do

  !  min and max block ID per level
  do nl =1,nlevs
    minID(nl) = 2**(nl-1)
    maxID(nl) = 2**nl -1
  end do

  ActiveBlocks(:) = -1
  ActiveBlocks(1) =  4
  ActiveBlocks(2) = 10
  ActiveBlocks(3) = 22
  ActiveBlocks(4) = 23
  ActiveBlocks(5) = 24
  ActiveBlocks(6) = 25
  ActiveBlocks(7) = 7
  ActiveBlocks(8) = 13
  lastActive = 8

  return
end subroutine init_mesh

!=======================================================================
!   Obtains the level of a block with bID
subroutine getLevel(bID,level)
  use globals, only : nlevs, minID, maxID
  implicit none
  integer, intent(in)  :: bID
  integer, intent(out) :: level
  integer  :: nl

  if (bID == -1) then
    level = -1
  else
    do nl=1,nlevs
      if ( (bID >= minID(nl)) .and.  (bID <= maxID(nl)) ) then
        level =nl
        return
      end if
    end do
  end if

  return
end subroutine

!=======================================================================
!   Obtains the position of the first cell in the block
subroutine getBlockPosition(bID,x)
  use globals, only : nx, dx, minID
  implicit none
  integer, intent(in) :: bID
  real, intent(out)   :: x
  integer :: level

  if (bID == -1) then
    x= -1
  else
    call getLevel(bID,level)
    x = (bID-minID(level)) * nx * dx(level)
  end if

  return
end subroutine getBlockPosition

!=======================================================================
!generates initial condition
subroutine initconds(time, tprint, itprint)
  use globals
  implicit none
  real, intent(out) :: time, tprint
  integer, intent (out) :: itprint
  integer :: nb, level, i
  real    :: xb, x
  real    :: xp = 0.45

  !  sweep all ActiveBlocks
  do nb = 1, lastActive
    !   get the level and position of first cell in each block
    if (ActiveBlocks(nb) /= -1) then

      call getLevel(ActiveBlocks(nb),level)
      call getBlockPosition(ActiveBlocks(nb),xb)

      do i = 0, nx+1
        x = xb + (real(i)-0.5)*dx(level)
        if (x < xp ) then
          u(1,i,nb)=1.0
          u(2,i,nb)=0.0
          u(3,i,nb)=1.0/(gamma-1.)
        else
          u(1,i,nb)=0.125
          u(2,i,nb)=0.0
          u(3,i,nb)=0.1/(gamma-1.)
        end if
        !
        if( (x-0.5*dx(level) < xp).and.(x+0.5*dx(level) > xp) ) then
          u(1,i,nb)=1.125/2.
          u(2,i,nb)=0.0
          u(3,i,nb)=1.1/2./(gamma-1.)
        end if

      end do

    end if

  end do

  !   reset the counters and time to 0
  time=0.
  tprint=0.
  itprint=0

  return
end subroutine initconds

!=======================================================================
! computes the primitives as a function of the Us
subroutine u2prim(u,prim)
  use globals, only : nx, nbmax, gamma, lastActive, ActiveBlocks
  implicit none
  real , intent(in)   :: u   (3,0:nx+1,nbmax)
  real , intent(out)  :: prim(3,0:nx+1,nbmax)
  integer :: i,nb

  do nb=1, lastActive
    if (ActiveBlocks(nb) /= -1 ) then

      do i=0,nx+1
        prim(1,i,nb) = u(1,i,nb)
        prim(2,i,nb) = u(2,i,nb)/u(1,i,nb)
        prim(3,i,nb) = (u(3,i,nb)-0.5*prim(1,i,nb)*prim(2,i,nb)**2)*(gamma-1.)
      end do

    endif
  end do

  return
end subroutine u2prim

!=======================================================================
! output to file
subroutine output(itprint)
  use globals
  implicit none
  integer, intent(in) :: itprint
  character (len=20) file1
  integer :: i, nb, level
  real    :: x

  ! open output file
  write(file1,'(a,i2.2,a)') 'euler_amr-',itprint,'.dat'
  open(unit=10,file=file1,status='unknown')

  ! writes x and rho, u(=vx) and P
  do nb=1,lastActive
    if (ActiveBlocks(nb) /= -1) then
      call getLevel(ActiveBlocks(nb),level)
      call getBlockPosition(ActiveBlocks(nb),x)
      do i=1,nx
        write(10,*) x+real((i)-0.5)*dx(level),prim(:,i,nb)
      end do
    end if
  end do

  ! closes output file
  close(10)

  return
end subroutine output
!=======================================================================
! computes the timestep allowed by the CFL criterium
subroutine timestep(dt)
  use globals
  implicit none
  real, intent(out) ::dt
  ! Courant number =0.9
  real, parameter :: Co=0.9
  real :: del, cs
  integer :: i, nb, lev

  del=1E30
  !   sweep all active blocks
  do nb = 1, lastActive
    if (ActiveBlocks(nb) /= -1) then

      call getLevel(ActiveBlocks(nb),lev)

      do i=0,nx+1
        cs=sqrt(gamma*prim(3,i,nb)/prim(1,i,nb))
        del=min( del,dx(lev)/(abs(prim(2,i,nb))+cs) )
      end do

    end if
  end do

  dt=Co*del

  return
end subroutine timestep

!=======================================================================
! integration from t to t+dt with the method of Lax
subroutine tstep(dt,time)
  use globals
  implicit none
  real, intent(in) :: dt, time
  real :: up(3,0:nx+1,nbmax), f(3,0:nx+1,nbmax)
  real :: dtx
  integer :: i, nb, level

  !  obtain the fluxes
  call hllcfluxes(prim,f)

  do nb=1,lastActive
    if (ActiveBlocks(nb) /= -1) then
      call getLevel(ActiveBlocks(nb),level)
      dtx=dt/dx(level)
      do i=1,nx
        up(:,i,nb)=u(:,i,nb)-dtx*(f(:,i,nb)-f(:,i-1,nb))
      end do

    end if
  end do

  !   Boundary conditions to the U^n+1
  call boundaries(up)

  ! copy the up to the u
  u(:,:,1:lastActive)=up(:,:,1:lastActive)

  return
end subroutine tstep

!=======================================================================
!  computes the HLLC fluxes in the entire domain
subroutine hllcfluxes(prim,f)
  use globals, only : gamma, nx, nbmax, ActiveBlocks, lastActive
  implicit none
  real,    intent(in) :: prim(3,0:nx+1,nbmax)
  real,    intent(out)::    f(3,0:nx+1,nbmax)
  integer :: i, nb
  real :: priml(3), primr(3), ff(3)

  do nb=1,lastActive
    if(ActiveBlocks(nb) /= -1) then
      do i=0,nx

        primL(:)= prim(:,i  ,nb)
        primR(:)= prim(:,i+1,nb)

        call prim2hllc(gamma, primL, primR, ff)
        f(:,i,nb)=ff(:)

      end do
    end if
  end do
end subroutine hllcfluxes

!=======================================================================
! Obtain the HLLC fluxes
subroutine prim2hllc(gamma,primL,primR,ff)
  implicit none
  real, intent(in) :: gamma, primL(3), primR(3)
  real, intent(out):: ff(3)
  real :: sl, sr, csr, csl, sst, fl(3), fr(3), ust(3), ek, uL(3), uR(3), uu(3)

  csl= sqrt(gamma*primL(3)/primL(1))
  csr= sqrt(gamma*primR(3)/primR(1))

  sl = min(primL(2)-csl, primR(2)-csr )
  sr = max(primL(2)+csl, primR(2)+csr )

  sst=( primR(3)-primL(3)+primL(1)*primL(2)*(sl-primL(2) )   &
                         -primR(1)*primR(2)*(sr-primR(2) ) )   &
           /( primL(1)*(sl-primL(2)) - primR(1)*(sr-primR(2)) )

  if(sl > 0. ) then

    call eulerfluxes(gamma,primL,ff)

  else if (sr < 0.) then

    call eulerfluxes(gamma,primR,ff)

  else if(sst >= 0.) then

    ust(1)=primL(1)*(sl-primL(2) )/(sl-sst)
    ust(2)=ust(1)*sst

    ek=(0.5*primL(1)*primL(2)**2)+primL(3)/(gamma-1.)
    ust(3)=ust(1)*(ek/primL(1)+ (sst-primL(2))*  &
         (sst+primL(3)/( primL(1)*(sl-primL(2) ) ) ) )

    call eulerfluxes(gamma,primL(:),fl)
    call prim2u     (gamma,primL(:),uu)

    ff(:)= fl(:)+ sl * (ust(:)-uu(:) )

  else if (sst <= 0.) then

    ust(1)=primR(1)*(sr-primR(2)) /(sr-sst)
    ust(2)=ust(1)*sst

    ek=(0.5*primR(1)*primR(2)**2)+primR(3)/(gamma-1.)
    ust(3)=ust(1)*(ek/primR(1)+ (sst-primR(2))*  &
         (sst+primR(3)/( primR(1)*(sr-primR(2)) ) ) )

    call eulerfluxes(gamma,primR(:),fr)
    call prim2u     (gamma,primR(:),uu)
    ff(:)= fr(:)+sr*(ust(:)-uu(:) )

  else
    print*,'Error in HLLC', sl, sr, sst
    stop
  endif

  return
end subroutine prim2hllc

!=======================================================================
!  computes the euler fluxes, one cell
subroutine eulerfluxes(gamma,pp,ff)
  implicit none
  real, intent(in) :: gamma, pp(3)
  real, intent(out):: ff(3)

  ff(1)=pp(1)*pp(2)
  ff(2)=pp(1)*pp(2)**2+pp(3)
  ff(3)=pp(2)*(0.5*pp(1)*pp(2)**2+gamma*pp(3)/(gamma-1.) )

  return
end subroutine eulerfluxes

!=======================================================================
! computes the primitives as a function of the Us, only in one cell
subroutine prim2u(gamma,pp,uu)
  implicit none
  real, intent(in)    :: gamma
  real , intent(in)   :: pp(3)
  real , intent(out)  :: uu(3)

    uu(1) = pp(1)
    uu(2) = pp(1)*pp(2)
    uu(3) = 0.5*pp(1)*pp(2)**2 +pp(3)/(gamma-1.)

  return
end subroutine prim2u

!=======================================================================
! Set boundary conditions
subroutine boundaries(u)
  use globals, only : lastActive, ActiveBlocks, nx, nbmax
  implicit none
  real,    intent(out):: u(3,0:nx+1,nbmax)
  integer :: i,nb

  !   Periodic boundary conditions
  !u(:,0 )=u(:,nx)
  !u(:,nx+1)=u(:,1)

  do nb=1, lastActive
    if (ActiveBlocks(nb) /= -1) then
      !  open boundary conditions
      u(:,0   ,nb)=u(:,1, nb)
      u(:,nx+1,nb)=u(:,nx,nb)
    end if
  end do

  return
end subroutine boundaries
