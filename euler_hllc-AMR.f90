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
  integer, parameter :: nx=10
  !   Number of levels allowed
  integer, parameter :: nlevs = 6
  !   maximumn number of blocks
  integer, parameter :: nbmax = 64
  !   Here we set the extent of X
  real, parameter :: xmax=1.
  real, parameter :: gamma=1.4
  !   This is a vector that contains u(x)
  real :: u(3,0:nx+1,nbmax), prim(3,0:nx+1,nbmax)
  !  other arrays and variables needed in AMR
  real    :: dx(nlevs)
  integer :: minID(nlevs), maxID(nlevs),ActiveBlocks(nbmax)
  integer :: lastActive
  logical :: FlagRefine(nbmax), FlagCoarse(nbmax)
  real, parameter :: rThresh = 0.1, cThresh = 0.01
  !
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
    call update_prim(u,prim)

    ! output at tprint intervals
    if(time >= tprint) then
      write(*,*) time,tmax,dt, itprint
      call output(itprint)
      tprint=tprint+dtprint
      itprint=itprint+1
    end if

    ! Obtain the $\Delta t$ allowed by the CFL criterium
    call timestep(dt)
    !
    ! Integrate u fom t to t+dt
    call tstep(dt,time)

    ! updates mesh
    call update_mesh()

    ! time counter increases
    time=time+dt

  end do

  stop
end program euler_amr

!=======================================================================
!  Initializes all things mesh
subroutine init_mesh()
  use globals, only : nx, nlevs, xmax, dx, minID, maxID, ActiveBlocks, &
                      lastActive, FlagRefine, FlagCoarse
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

  FlagRefine(:) = .false.
  FlagCoarse(:) = .false.

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
!   Obtains nb from bID (or -1 in case is not active)
subroutine get_nb(bID, nb)
  use globals, only : lastActive, ActiveBlocks
  implicit none
  integer, intent(in)  :: bID
  integer, intent(out) :: nb

  do nb=1, lastActive
    if (ActiveBlocks(nb)==bID) return
  end do

  nb = -1

  return
end subroutine get_nb

!=======================================================================
!   Obtains the nb from the neighbor at the Left
subroutine getLeft(nb, nbLeft)
  use globals, only : ActiveBlocks, minID
  implicit none
  integer, intent(in)  :: nb
  integer, intent(out) :: nbLeft
  integer :: selfID, fatherID, level, sonID

  !   get bID of current block
  selfID = ActiveBlocks(nb)

  !   domain boundary
  call getLevel(selfID,level)
  if (selfID == minID(level)) then
    nbLeft = -1
    return
  end if

  !  same level of refinement
  call get_nb(selfID-1 ,nbLeft)
  if (nbLeft /= -1) return

  !  lower level of refinement
  fatherID = selfID/2
  call get_nb(fatherID-1 ,nbLeft)
  if (nbLeft /= -1) return

  !  higher level or refinemnt
  sonID = selfID*2
  call get_nb(sonID-1,nbLeft)

  return
end subroutine getLeft

!=======================================================================
!   Obtains the nb from the neighbor at the Right
subroutine getRight(nb, nbRight)
  use globals, only : ActiveBlocks, maxID
  implicit none
  integer, intent(in)  :: nb
  integer, intent(out) :: nbRight
  integer :: selfID, fatherID, level, sonID

  !   get bID of current block
  selfID = ActiveBlocks(nb)

  !   domain boundary
  call getLevel(selfID,level)
  if (selfID == maxID(level)) then
    nbRight = -1
    return
  end if

  !  same level of refinement
  call get_nb(selfID+1 ,nbRight)
  if (nbRight /= -1) return

  !  lower level of refinement
  fatherID = selfID/2
  call get_nb(fatherID+1 ,nbRight)
  if (nbRight /= -1) return

  !  higher level or refinemnt
  sonID = selfID*2
  call get_nb(sonID+2,nbRight)

  return
end subroutine getRight

!=======================================================================
!generates initial condition
subroutine initconds(time, tprint, itprint)
  use globals
  implicit none
  real, intent(out) :: time, tprint
  integer, intent (out) :: itprint
  integer :: nb, level, i
  real    :: xb, x
  real    :: xp = 0.5

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
! computes the primitives as a function of the Us in all  ActiveBlocks
subroutine update_prim(u,prim)
  use globals, only : nx, nbmax, gamma, lastActive, ActiveBlocks
  implicit none
  real , intent(in)   :: u   (3,0:nx+1,nbmax)
  real , intent(out)  :: prim(3,0:nx+1,nbmax)
  integer :: i,nb

  do nb=1, lastActive
    if (ActiveBlocks(nb) /= -1 ) then

      do i=0,nx+1
        call u2prim(gamma,u(:,i,nb),prim(:,i,nb))
      end do

    endif
  end do

  return
end subroutine update_prim

!=======================================================================
! Computes primitives from the conserved vars in a single cell
subroutine u2prim(gamma,uu,pp)
  implicit none
  real, intent(in)   :: gamma
  real , intent(in)  :: uu(3)
  real , intent(out) :: pp(3)

  pp(1) = uu(1)
  pp(2) = uu(2)/uu(1)
  pp(3) = (uu(3)-0.5*pp(1)*pp(2)**2)*(gamma-1.)

return
end subroutine

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
! computes the timestep allowed by the CFL
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
  integer :: nbL, nbR

  !   Periodic boundary conditions
  !u(:,0 )=u(:,nx)
  !u(:,nx+1)=u(:,1)

  do nb=1, lastActive
    if (ActiveBlocks(nb) /= -1) then

      !  open boundary conditions
      !  left
      call getLeft(nb,nbL)
      if (nbL == -1) then
        u(:,0 ,nb)=u(:,1, nb )
      else
        u(:,0 ,nb)=u(:,nx,nbL)
      endif

      !  right
      call getRight(nb,nbR)
      if (nbR == -1) then
        u(:,nx+1,nb)=u(:,nx,nb )
      else
        u(:,nx+1,nb)=u(:,1 ,nbR)
      endif

    end if
  end do

  return
end subroutine boundaries

!=======================================================================
! Mark (set flag) blocks for refinement/coarsening based on then
! gradients of density and pressure
subroutine FlagGrads(prim)
  use globals, only : FlagRefine, FlagCoarse, nbmax, rThresh, cThresh, &
                      lastActive, ActiveBlocks, dx, nlevs, nx
  implicit none
  real, intent(in):: prim(3,0:nx+1,nbmax)
  integer :: nb, level, i
  real    :: gradP, gradrho, maxgrad

  do nb = 1, lastActive
    maxgrad = 0.
    if(ActiveBlocks(nb) /= -1 ) then

      call getLevel(ActiveBlocks(nb),level)
      do i=1,nx

        gradrho = (prim(1,i+1,nb) -prim(1,i-1,nb))/(prim(1,i,nb)*2.*dx(level))
        gradP =   (prim(3,i+1,nb) -prim(3,i-1,nb))/(prim(3,i,nb)*2.*dx(level))

        maxgrad = max(maxgrad, gradrho)
        maxgrad = max(maxgrad, gradP  )

        if(maxgrad >= rThresh) then
          !  mark for refinement (only it not at max resolution already)
          !  and stop checking
          if (level < nlevs) FlagRefine(nb) = .true.
          exit
        end if
      end do

      if (maxgrad <= cThresh) then
        FlagCoarse(nb) = .true.
      end if

    end if
  end do


  return
end subroutine

!=======================================================================
! Mark blocks for refinement based on proximity criteria
subroutine FLagProx()
  use globals, only : ActiveBlocks, lastActive, nlevs, FlagRefine, &
                      FlagCoarse
  implicit none
  integer  :: nb, ilev, dadID, nbLeft, nbRight

  do ilev=nlevs,1,-1
    do nb = 1, lastActive
      if(ActiveBlocks(nb) /= -1 ) then

        !  check all blocks marked for refinement and mark the neighbor
        !  for refinement (and inhibit coarsening) if it is at a lower
        !  level, just inhibit refinement if it is at same level
        if(FlagRefine(nb)) then
          dadID = ActiveBlocks(nb)
          !  same level neighbors
          !  left
          call get_nb(dadID-1,nbLeft)
          if (nbLeft  /= -1) FlagCoarse(nbLeft ) = .false.
          !  right
          call get_nb(dadID+1,nbRight)
          if (nbRight /= -1) FlagCoarse(nbRight) = .false.

          !  lower level neighbors
          !  left
          call get_nb(dadID/2-1,nbLeft)
          if (nbLeft  /= -1) then
            FlagRefine(nbLeft ) = .true.
            FlagCoarse(nbLeft ) = .false.
          end if
          !  right
          call get_nb(dadID/2+1,nbRight)
          if (nbRight  /= -1) then
            FlagRefine(nbRight) = .true.
            FlagCoarse(nbRight) = .false.
          end if

        end if

        !  Inhibit coarsening if finer neighbors are not set for coarsening
        if(FlagCoarse(nb)) then
          !  finer neighbors
          !  left
          call get_nb(dadID*2-1,nbLeft)
          if(nbLeft /= -1) then
            if(.not.FlagCoarse(nbLeft)) FlagCoarse(nb)=.false.
          end if
          !  right
          call get_nb(dadID*2+2,nbRight)
          if(nbRight /= -1) then
            if(.not.FlagCoarse(nbRight)) FlagCoarse(nb)=.false.
          end if

        end if

      end if
    end do
  end do

 return
end subroutine FLagProx

!=======================================================================
!  Refines block nb to next level, refinement is applied both to the
!  primitives and conserved variables
subroutine refineBlock(dadNb)
  use globals, only : ActiveBlocks, lastActive, nx, nbmax, u, prim, &
                      gamma, FlagRefine
  implicit none
  integer, intent(in)  :: dadNb
  integer :: dadID, son1ID, son2ID, son1nb, son2nb, nb, i

  dadID = ActiveBlocks(dadNb)
  son1ID = dadID*2
  son2ID = son1ID + 1
  !  change the bID of the parent for that of the first son/
  ActiveBlocks(dadNb) = son1ID
  son1nb =dadNb
  !  Activate the bID of rthe second son in the first available spot
  !  increase lastActive if needed
  son2nb = -1 !  to test if nbmax exceeded what is allowed
  do nb=1,nbmax
    if (ActiveBlocks(nb) == -1) then
      ActiveBlocks(nb) = son2ID
      son2nb = nb
      if(nb > lastActive) then
        lastActive = nb
      end if
      exit
    end if
  end do
  if (son2nb == -1) stop 'nb needed ecxeeded nbmax'

  !  copy data from dad to sons and update primitives
  !  Second son
  do i = 0, nx+1
    u(:,i,son2nb) = u(:,(i+1+nx)/2,dadNb)
    call u2prim(gamma,u(:,i,son2nb),prim(:,i,son2nb))
  end do
  !  first son
  do i=nx+1,0,-1
    u(:,i,son1nb) = u(:,(i+1   )/2,dadNb)
    call u2prim(gamma,u(:,i,son1nb),prim(:,i,son1nb))
  end do

  !  reset refining flag
  FlagRefine(dadNb) = .false.

  return
end subroutine refineBlock

!=======================================================================
!  updates (refines and coarsens) mesh
subroutine update_mesh()
  use globals, only : prim, u, ActiveBlocks, lastActive, FlagRefine, FlagCoarse
  implicit none
  integer  :: nb

  !  Mark by physical criteria
  call FlagGrads(prim)
  !  Mark by proximity
  call FLagProx()

  !  Proceed w/refinement of marked blocks
  do nb=1,lastActive
    if (ActiveBlocks(nb) /= -1) then
      if(FlagRefine(nb)) then
        print*, 'refining', ActiveBlocks(nb)
        call refineBlock(nb)
        print*, lastActive
      end if
    end if
  end do


  return
end subroutine update_mesh
