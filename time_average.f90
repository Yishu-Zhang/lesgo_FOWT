!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.

!*******************************************************************************
module time_average
!*******************************************************************************
use types, only : rprec
use param, only : nx, ny, nz, lbz
#ifdef PPCGNS
use cgns
#endif

private
public :: tavg_t

real(rprec), allocatable, dimension(:,:,:) :: w_uv, u_w, v_w, dpdz_uv
real(rprec), allocatable, dimension(:,:,:) :: vortx, vorty, vortz, vorty2, vortz2
real(rprec), allocatable, dimension(:,:,:) :: pres_real ,v_omegaz2,w_omegay2
real(rprec), dimension(:,:,:), allocatable :: vdudy, vdvdx, wdudz, wdwdx, udvdy, udwdz
real(rprec), allocatable, dimension(:,:,:) :: pdetadx, prealdetadx
real(rprec), allocatable, dimension(:,:) :: pwalldetadx1, pwalldetadx2
!real(rprec), allocatable, dimension(:,:,:) :: B_uu, B_uu1, E_uu
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
real(rprec), allocatable, dimension(:,:,:) :: fza_uv
#endif

!  Sums performed over time
type tavg_t
    real(rprec), dimension(:,:,:), allocatable :: u, v, w, u_w, v_w, w_uv
    real(rprec), dimension(:,:,:), allocatable :: u2, v2, w2, uv, uw, vw
    real(rprec), dimension(:,:,:), allocatable :: txx, tyy, tzz, txy, txz, tyz
    real(rprec), dimension(:,:), allocatable :: eqmxz, wpmxz, eqmyz, wpmyz, pwalldetadx1, pwalldetadx2
    real(rprec), dimension(:,:,:), allocatable :: pdetadx, prealdetadx
    real(rprec), dimension(:,:,:), allocatable :: p, fx, fy, fz, dpdz, dpdz_uv
    real(rprec), dimension(:,:,:), allocatable :: cs_opt2, pres_real
    real(rprec), dimension(:,:,:), allocatable :: vortx, vorty, vortz
    real(rprec), dimension(:,:,:), allocatable :: dudz,dvdx,dudy,dwdx,vorty2,vortz2,v_omegaz2,w_omegay2
    real(rprec), dimension(:,:,:), allocatable :: vdudy, vdvdx, wdudz, wdwdx, udvdy, udwdz, dwdz
    real(rprec) :: total_time
    real(rprec), dimension(:,:,:), allocatable :: cc_z, dpdx, divtx, RHSx, RHSx_f, v_omegaz, u_omegay, dtxdx, dtydy, dtzdz
!    real(rprec), allocatable, dimension(:,:,:) :: B_uu, B_uu1, E_uu
    ! Time between calls of tavg_compute, built by summing dt
    real(rprec) :: dt
    ! Switch for determining if time averaging has been initialized
    logical :: initialized = .false.
contains
    procedure, public :: init
    procedure, public :: compute
    procedure, public :: finalize
    procedure, public :: checkpoint
end type tavg_t

character(:), allocatable :: checkpoint_tavg_file

contains

!*******************************************************************************
subroutine init(this)
!*******************************************************************************
use messages
use string_util
use param, only : read_endian, coord, path
implicit none

class(tavg_t), intent(inout) :: this

character(:), allocatable :: ftavg_in
#ifdef PPMPI
character (*), parameter :: MPI_suffix = '.c'
#endif
character (128) :: fname
! integer :: i, j, k

logical :: exst

! Create file name
allocate(ftavg_in, source = path // 'tavg.out')
allocate(checkpoint_tavg_file, source=path // 'tavg.out')

! Allocate and initialize
allocate( this%u(nx,ny,lbz:nz) ); this%u(:,:,:) = 0._rprec
allocate( this%v(nx,ny,lbz:nz) ); this%v(:,:,:) = 0._rprec
allocate( this%w(nx,ny,lbz:nz) ); this%w(:,:,:) = 0._rprec
allocate( this%u_w(nx,ny,lbz:nz) ); this%u_w(:,:,:) = 0._rprec
allocate( this%v_w(nx,ny,lbz:nz) ); this%v_w(:,:,:) = 0._rprec
allocate( this%w_uv(nx,ny,lbz:nz) ); this%w_uv(:,:,:) = 0._rprec
allocate( this%u2(nx,ny,lbz:nz) ); this%u2(:,:,:) = 0._rprec
allocate( this%v2(nx,ny,lbz:nz) ); this%v2(:,:,:) = 0._rprec
allocate( this%w2(nx,ny,lbz:nz) ); this%w2(:,:,:) = 0._rprec
allocate( this%uv(nx,ny,lbz:nz) ); this%uv(:,:,:) = 0._rprec
allocate( this%uw(nx,ny,lbz:nz) ); this%uw(:,:,:) = 0._rprec
allocate( this%vw(nx,ny,lbz:nz) ); this%vw(:,:,:) = 0._rprec
allocate( this%txx(nx,ny,lbz:nz) ); this%txx(:,:,:) = 0._rprec
allocate( this%tyy(nx,ny,lbz:nz) ); this%tyy(:,:,:) = 0._rprec
allocate( this%tzz(nx,ny,lbz:nz) ); this%tzz(:,:,:) = 0._rprec
allocate( this%txy(nx,ny,lbz:nz) ); this%txy(:,:,:) = 0._rprec
allocate( this%txz(nx,ny,lbz:nz) ); this%txz(:,:,:) = 0._rprec
allocate( this%tyz(nx,ny,lbz:nz) ); this%tyz(:,:,:) = 0._rprec

allocate( this%eqmxz(nx,ny) ); this%eqmxz(:,:) = 0._rprec
allocate( this%wpmxz(nx,ny) ); this%wpmxz(:,:) = 0._rprec
allocate( this%eqmyz(nx,ny) ); this%eqmyz(:,:) = 0._rprec
allocate( this%wpmyz(nx,ny) ); this%wpmyz(:,:) = 0._rprec

allocate( this%p(nx,ny,lbz:nz) ); this%p(:,:,:) = 0._rprec
allocate( this%pres_real(nx,ny,lbz:nz) ); this%pres_real(:,:,:) = 0._rprec
!allocate( this%dpdz(nx,ny,lbz:nz) ); this%dpdz(:,:,:) = 0._rprec
!allocate( this%dpdz_uv(nx,ny,lbz:nz) ); this%dpdz_uv(:,:,:) = 0._rprec
allocate( this%fx(nx,ny,lbz:nz) ); this%fx(:,:,:) = 0._rprec
allocate( this%fy(nx,ny,lbz:nz) ); this%fy(:,:,:) = 0._rprec
allocate( this%fz(nx,ny,lbz:nz) ); this%fz(:,:,:) = 0._rprec
allocate( this%cs_opt2(nx,ny,lbz:nz) ); this%cs_opt2(:,:,:) = 0._rprec
allocate( this%vortx(nx,ny,lbz:nz) ); this%vortx(:,:,:) = 0._rprec
allocate( this%vorty(nx,ny,lbz:nz) ); this%vorty(:,:,:) = 0._rprec
allocate( this%vortz(nx,ny,lbz:nz) ); this%vortz(:,:,:) = 0._rprec


allocate( this%cc_z(nx,ny,lbz:nz) ); this%cc_z(:,:,:) = 0._rprec
!allocate( this%dpdx(nx,ny,lbz:nz) ); this%dpdx(:,:,:) = 0._rprec
allocate( this%divtx(nx,ny,lbz:nz) ); this%divtx(:,:,:) = 0._rprec
allocate( this%RHSx(nx,ny,lbz:nz) ); this%RHSx(:,:,:) = 0._rprec
allocate( this%RHSx_f(nx,ny,lbz:nz) ); this%RHSx_f(:,:,:) = 0._rprec
allocate( this%v_omegaz(nx,ny,lbz:nz) ); this%v_omegaz(:,:,:) = 0._rprec
allocate( this%u_omegay(nx,ny,lbz:nz) ); this%u_omegay(:,:,:) = 0._rprec
allocate( this%dtxdx(nx,ny,lbz:nz) ); this%dtxdx(:,:,:) = 0._rprec
allocate( this%dtydy(nx,ny,lbz:nz) ); this%dtydy(:,:,:) = 0._rprec
allocate( this%dtzdz(nx,ny,lbz:nz) ); this%dtzdz(:,:,:) = 0._rprec

allocate( this%dudy(nx,ny,lbz:nz) ); this%dudy(:,:,:) = 0._rprec
allocate( this%dvdx(nx,ny,lbz:nz) ); this%dvdx(:,:,:) = 0._rprec
allocate( this%dudz(nx,ny,lbz:nz) ); this%dudz(:,:,:) = 0._rprec
allocate( this%dwdx(nx,ny,lbz:nz) ); this%dwdx(:,:,:) = 0._rprec
allocate( this%vorty2(nx,ny,lbz:nz) ); this%vorty2(:,:,:) = 0._rprec
allocate( this%vortz2(nx,ny,lbz:nz) ); this%vortz2(:,:,:) = 0._rprec
allocate( this%v_omegaz2(nx,ny,lbz:nz) ); this%v_omegaz2(:,:,:) = 0._rprec
allocate( this%w_omegay2(nx,ny,lbz:nz) ); this%w_omegay2(:,:,:) = 0._rprec

allocate( this%vdudy(nx,ny,lbz:nz) ); this%vdudy(:,:,:) = 0._rprec
allocate( this%vdvdx(nx,ny,lbz:nz) ); this%vdvdx(:,:,:) = 0._rprec
allocate( this%wdudz(nx,ny,lbz:nz) ); this%wdudz(:,:,:) = 0._rprec
allocate( this%wdwdx(nx,ny,lbz:nz) ); this%wdwdx(:,:,:) = 0._rprec
allocate( this%udvdy(nx,ny,lbz:nz) ); this%udvdy(:,:,:) = 0._rprec
allocate( this%udwdz(nx,ny,lbz:nz) ); this%udwdz(:,:,:) = 0._rprec
allocate( this%dwdz(nx,ny,lbz:nz) ); this%dwdz(:,:,:) = 0._rprec

allocate( this%pdetadx(nx,ny,lbz:nz) ); this%pdetadx(:,:,:) = 0._rprec
allocate( this%prealdetadx(nx,ny,lbz:nz) ); this%prealdetadx(:,:,:) = 0._rprec

allocate( this%pwalldetadx1(nx,ny) ); this%pwalldetadx1(:,:) = 0._rprec
allocate( this%pwalldetadx2(nx,ny) ); this%pwalldetadx2(:,:) = 0._rprec

!allocate( this%B_uu(nx,ny,lbz:nz) ); this%B_uu(:,:,:) = 0._rprec
!allocate( this%B_uu1(nx,ny,lbz:nz) ); this%B_uu1(:,:,:) = 0._rprec
!allocate( this%E_uu(nx,ny,lbz:nz) ); this%E_uu(:,:,:) = 0._rprec
fname = ftavg_in
#ifdef PPMPI
call string_concat(fname, MPI_suffix, coord)
#endif

inquire (file=fname, exist=exst)
if (.not. exst) then
    !  Nothing to read in
    if (coord == 0) then
        write(*,*) ' '
        write(*,*)'No previous time averaged data - starting from scratch.'
    end if
    this%total_time = 0._rprec
else
    open(1, file=fname, action='read', position='rewind', form='unformatted',  &
        convert=read_endian)
    read(1) this%total_time
    read(1) this%u
    read(1) this%v
    read(1) this%w
    read(1) this%u_w
    read(1) this%v_w
    read(1) this%w_uv
    read(1) this%u2
    read(1) this%v2
    read(1) this%w2
    read(1) this%uv
    read(1) this%uw
    read(1) this%vw
    read(1) this%txx
    read(1) this%tyy
    read(1) this%tzz
    read(1) this%txy
    read(1) this%txz
    read(1) this%tyz
!    read(1) this%dpdz
!    read(1) this%dpdz_uv
    read(1) this%p
    read(1) this%pres_real

    read(1) this%eqmxz
    read(1) this%wpmxz
    read(1) this%eqmyz
    read(1) this%wpmyz

    read(1) this%fx
    read(1) this%fy
    read(1) this%fz
    read(1) this%cs_opt2
    read(1) this%vortx
    read(1) this%vorty
    read(1) this%vortz

    read(1) this%cc_z
!    read(1) this%dpdx
    read(1) this%divtx
    read(1) this%RHSx
    read(1) this%RHSx_f
    read(1) this%v_omegaz
    read(1) this%u_omegay
    read(1) this%dtxdx
    read(1) this%dtydy
    read(1) this%dtzdz

    read(1) this%v_omegaz2
    read(1) this%w_omegay2
    read(1) this%dudz
    read(1) this%dvdx
    read(1) this%dudy
    read(1) this%dwdx
    read(1) this%vorty2
    read(1) this%vortz2
    
    read(1) this%wdudz
    read(1) this%vdvdx
    read(1) this%vdudy
    read(1) this%wdwdx
    read(1) this%udvdy
    read(1) this%udwdz
    read(1) this%dwdz

    read(1) this%pdetadx
    read(1) this%prealdetadx
    read(1) this%pwalldetadx1
    read(1) this%pwalldetadx2

!    read(1) this%B_uu
!    read(1) this%B_uu1
!    read(1) this%E_uu
    close(1)
end if

! Allocate computation arrays
allocate(w_uv(nx,ny,lbz:nz), u_w(nx,ny,lbz:nz), v_w(nx,ny,lbz:nz),dpdz_uv(nx,ny,lbz:nz))
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
allocate(fza_uv(nx,ny,lbz:nz))
#endif
allocate(pres_real(nx,ny,lbz:nz))
allocate(vortx(nx,ny,lbz:nz), vorty(nx,ny,lbz:nz), vortz(nx,ny,lbz:nz))
allocate(vorty2(nx,ny,lbz:nz), vortz2(nx,ny,lbz:nz))
allocate(v_omegaz2(nx,ny,lbz:nz), w_omegay2(nx,ny,lbz:nz))

allocate(vdudy(nx,ny,lbz:nz), vdvdx(nx,ny,lbz:nz), wdudz(nx,ny,lbz:nz), wdwdx(nx,ny,lbz:nz) )
allocate(udvdy(nx,ny,lbz:nz), udwdz(nx,ny,lbz:nz))

allocate(pdetadx(nx,ny,lbz:nz))
allocate(prealdetadx(nx,ny,lbz:nz))
allocate(pwalldetadx1(nx,ny))
allocate(pwalldetadx2(nx,ny))

!allocate(B_uu(nx,ny,lbz:nz))
!allocate(B_uu1(nx,ny,lbz:nz))
!allocate(E_uu(nx,ny,lbz:nz))
! Initialize dt
this%dt = 0._rprec

! Set global switch that tavg as been initialized
this%initialized = .true.

end subroutine init

!*******************************************************************************
subroutine compute(this)
!*******************************************************************************
!
!  This subroutine collects the stats for each flow
!  variable quantity
!
use param, only : ubc_mom, lbc_mom, coord, nproc, dz, nx, ny
use sgs_param, only : Cs_opt2
use sim_param, only : u, v, w, p
use sim_param, only : txx, txy, tyy, txz, tyz, tzz
use sim_param, only : eqmxz, s_wpmxz, eqmyz, s_wpmyz
use sim_param, only : detadx, u_orb, divtz
use sim_param, only : dudy, dudz, dvdx, dvdz, dwdx, dwdy, dpdz, dvdy, dwdz
use sim_param, only : cc_z, dpdx, divtx, RHSx, RHSx_f, v_omegaz, u_omegay, dtxdx, dtydy, dtzdz
#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
use sim_param, only : fxa, fya, fza
#endif
use functions, only : interp_to_uv_grid, interp_to_w_grid
!use fft
implicit none
!include'fftw3.f'

class(tavg_t), intent(inout) :: this

integer :: jz, jz_min, jz_max, jx, jy
!integer*8 plan
!real(rprec) :: u_xyAvg, uu_xyAvg, uu_xyAvg_2, u_rms
!real(rprec), dimension(:,:), allocatable :: u_inst, u_fluc, B_uu_xy
!complex(rprec), dimension(:,:), allocatable :: u_fluc_hat, B_uu_hat 

! Interpolation onto other grids
w_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz )
u_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(u(1:nx,1:ny,lbz:nz), lbz )
v_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(v(1:nx,1:ny,lbz:nz), lbz )

!dpdz_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(dpdz(1:nx,1:ny,lbz:nz), lbz )

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
fza_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(fza(1:nx,1:ny,lbz:nz), lbz )
#endif

! Calculate real pressure (not pseudo-pressure)
pres_real(1:nx,1:ny,lbz:nz) = 0._rprec
pres_real(1:nx,1:ny,lbz:nz) = p(1:nx,1:ny,lbz:nz)                              &
    - 0.5 * ( u(1:nx,1:ny,lbz:nz)**2 + w_uv(1:nx,1:ny,lbz:nz)**2               &
    + v(1:nx,1:ny,lbz:nz)**2 )

! Compute vorticity
! Use vortx as an intermediate step for performing uv-w interpolation
! Vorticity is written in w grid
vortz2(1:nx,1:ny,lbz:nz) = dvdx(1:nx,1:ny,lbz:nz) - dudy(1:nx,1:ny,lbz:nz)
vortz(1:nx,1:ny,lbz:nz) = interp_to_w_grid( vortz2(1:nx,1:ny,lbz:nz), lbz)
vortx(1:nx,1:ny,lbz:nz) = dwdy(1:nx,1:ny,lbz:nz) - dvdz(1:nx,1:ny,lbz:nz)
vorty(1:nx,1:ny,lbz:nz) = dudz(1:nx,1:ny,lbz:nz) - dwdx(1:nx,1:ny,lbz:nz)
vorty2(1:nx,1:ny,lbz:nz) = vorty(1:nx,1:ny,lbz:nz) 
 
if (coord == 0) then
    vortz(1:nx,1:ny, 1) = 0._rprec
    vorty2(1:nx,1:ny, 1) =  dudz(1:nx,1:ny,1) - 0.5_rprec*(dwdx(1:nx,1:ny,1) + dwdx(1:nx,1:ny,2))
end if
if (coord == nproc-1) then
    vorty2(1:nx,1:ny, nz) =  dudz(1:nx,1:ny,nz-1) - 0.5_rprec*(dwdx(1:nx,1:ny,nz-1) + dwdx(1:nx,1:ny,nz))
end if


! Compute each nonlinear term only for x-mom the same way its done in convec.f90 --MA
    v_omegaz2(1:nx,1:ny,lbz:nz) =  v(1:nx,1:ny,lbz:nz)*(-vortz2(1:nx,1:ny,lbz:nz))

if (coord == 0) then
    v_omegaz2(1:nx,1:ny, 1) =  v(1:nx,1:ny,1)*(-vortz2(1:nx,1:ny,1))
end if

if (coord == 0) then
     w_omegay2(1:nx,1:ny,1) = 0.5_rprec*w(1:nx,1:ny,2)*(vorty2(1:nx,1:ny,2))
     pdetadx(1:nx,1:ny,1) = p(1:nx,1:ny,1)*detadx(1:nx,1:ny)
     prealdetadx(1:nx,1:ny,1) = pres_real(1:nx,1:ny,1)*detadx(1:nx,1:ny)
     jz_min = 2
else
    jz_min = 1
end if

if (coord == nproc-1 ) then 
      w_omegay2(1:nx,1:ny,nz-1) = (0.5_rprec*w(1:nx,1:ny,nz-1)*(vorty2(1:nx,1:ny,nz-1)))
      pdetadx(1:nx,1:ny,nz-1) = p(1:nx,1:ny,nz-1)*detadx(1:nx,1:ny)
      prealdetadx(1:nx,1:ny,nz-1) = pres_real(1:nx,1:ny,nz-1)*detadx(1:nx,1:ny)
      jz_max = nz-2
else
    jz_max = nz-1
end if

do jz = jz_min, jz_max    
    w_omegay2(:,:,jz)= 0.5_rprec*(w(1:nx,1:ny,jz+1)*(vorty2(1:nx,1:ny,jz+1))&
        +w(1:nx,1:ny,jz)*(vorty2(1:nx,1:ny,jz)))
    pdetadx(1:nx,1:ny,jz) = p(1:nx,1:ny,jz)*detadx(1:nx,1:ny)
    prealdetadx(1:nx,1:ny,jz) = pres_real(1:nx,1:ny,jz)*detadx(1:nx,1:ny)
end do

! Computing the pressure at the wall by using forward finite difference 
! 1 means im using div(tau_z) at wall and 2 means i use pressure gradient (should be same)
! I am assuming dpdz is always on w grid 

pwalldetadx1(1:nx,1:ny) = (p(1:nx,1:ny,1) + (dz/2)*divtz(1:nx,1:ny,1))*detadx(1:nx,1:ny)
pwalldetadx2(1:nx,1:ny) = (pres_real(1:nx,1:ny,1) - (dz/2)*dpdz(1:nx,1:ny,1))*detadx(1:nx,1:ny)

! Computing vdudy, vdvdx, wdudz, wdwdx terms
vdudy(1:nx,1:ny,lbz:nz) = v(1:nx,1:ny,lbz:nz)*dudy(1:nx,1:ny,lbz:nz)
vdvdx(1:nx,1:ny,lbz:nz) = v(1:nx,1:ny,lbz:nz)*dvdx(1:nx,1:ny,lbz:nz)
wdwdx(1:nx,1:ny,lbz:nz) = w(1:nx,1:ny,lbz:nz)*dwdx(1:nx,1:ny,lbz:nz)
wdudz(1:nx,1:ny,lbz:nz) = w(1:nx,1:ny,lbz:nz)*dudz(1:nx,1:ny,lbz:nz)
udvdy(1:nx,1:ny,lbz:nz) = u(1:nx,1:ny,lbz:nz)*dvdy(1:nx,1:ny,lbz:nz)
udwdz(1:nx,1:ny,lbz:nz) = u(1:nx,1:ny,lbz:nz)*dwdz(1:nx,1:ny,lbz:nz)

! Since dudz at first level is at uv grid then 
if (coord == 0) then

wdudz(1:nx,1:ny,1) = 0.5_rprec*(w(1:nx,1:ny,2)+w(1:nx,1:ny,1))*dudz(1:nx,1:ny,1) 

end if

!******** Computating the correlation function at each XY plane *********!

! Getting instantaneous streamwise velocity at each z location
! and taking the spatial average

!allocate(u_fluc(nx,ny))
!allocate(u_fluc_hat(nx,ny))
!allocate(B_uu_hat(nx,ny))
!allocate(B_uu_xy(nx,ny))

!do jz = 1,nz
   
!u_xyAvg = 0._rprec
!   do jx = 1, nx 
!       do jy = 1, ny
!       u_xyAvg = u_xyAvg + u(jx,jy,jz)
!       enddo
!    enddo
!u_xyAvg =  u_xyAvg/(nx*ny)

!write(*,*) 'xyAvg velocity', u_xyAvg 

! Calculating the fluctuating velocity u'
!u_fluc(:,:) = u(:,:,jz) - u_xyAvg

! Calculating average u'u' --> B_uu(0) used for normalization

!uu_xyAvg = 0._rprec
!uu_xyAvg_2 = 0._rprec
!    do jx = 1, nx
!       do jy = 1, ny
!       uu_xyAvg = uu_xyAvg + u_fluc(jx,jy)*u_fluc(jx,jy)
!       uu_xyAvg_2 = uu_xyAvg_2 + u_inst(jx,jy)*u_inst(jx,jy)
!       enddo
!    enddo
!uu_xyAvg = uu_xyAvg/(nx*ny)
!u_rms = sqrt(uu_xyAvg_2 - u_xyAvg**2)

!write(*,*) 'xyAvg RS', uu_xyAvg
!write(*,*) 'u RMS', u_RMS

! Applying FFT 2D to the fluctuating velocity and do the convolution

!u_fluc(:,:) = u_fluc(:,:)/(nx*ny)

!write(*,*) 'Forward FFT '
!call dfftw_plan_dft_r2c_2d(plan, nx, ny, u_fluc, u_fluc_hat, FFTW_ESTIMATE)
!call dfftw_execute_dft_r2c(plan, u_fluc, u_fluc_hat)
!call dfftw_destroy_plan(plan)

!B_uu_hat = real(u_fluc_hat * conjg(u_fluc_hat)) 

!write(*,*) 'Backward FFT '
!call dfftw_plan_dft_c2r_2d(plan, nx, ny, B_uu_hat, B_uu_xy, FFTW_ESTIMATE)
!call dfftw_execute_dft_c2r(plan, B_uu_hat, B_uu_xy)
!call dfftw_destroy_plan(plan)


!write(*,*) 'Finish'
!B_uu(:,:,jz) = B_uu_xy(:,:)
!B_uu1(:,:,jz) = B_uu_xy(:,:)/uu_xyAvg
!E_uu(:,:,jz) = B_uu_hat(:,:)

!end do



! note: u_w not necessarily zero on walls, but only mult by w=0 vu u'w', so OK
! can zero u_w at BC anyway:
if(coord==0       .and. lbc_mom>0) u_w(:,:,1)  = 0._rprec
if(coord==nproc-1 .and. ubc_mom>0) u_w(:,:,nz) = 0._rprec
if(coord==0       .and. lbc_mom>0) v_w(:,:,1)  = 0._rprec
if(coord==nproc-1 .and. ubc_mom>0) v_w(:,:,nz) = 0._rprec

this%u(:,:,:) = this%u(:,:,:) + u(1:nx,1:ny,:)*this%dt          !! uv grid
this%v(:,:,:) = this%v(:,:,:) + v(1:nx,1:ny,:)*this%dt          !! uv grid
this%w(:,:,:) = this%w(:,:,:) + w(1:nx,1:ny,:)*this%dt          !! w grid
this%w_uv(:,:,:) = this%w_uv(:,:,:) + w_uv(1:nx,1:ny,:)*this%dt !! uv grid
this%u_w(:,:,:) = this%u_w(:,:,:) + u_w(1:nx,1:ny,:)*this%dt    !! w grid
this%v_w(:,:,:) = this%v_w(:,:,:) + v_w(1:nx,1:ny,:)*this%dt    !! w grid

this%u2(:,:,:) = this%u2(:,:,:) + u(1:nx,1:ny,:)*u(1:nx,1:ny,:)*this%dt     !! uv grid
this%v2(:,:,:) = this%v2(:,:,:) + v(1:nx,1:ny,:)*v(1:nx,1:ny,:)*this%dt     !! uv grid
this%w2(:,:,:) = this%w2(:,:,:) + w(1:nx,1:ny,:)*w(1:nx,1:ny,:)*this%dt     !! w grid
this%uv(:,:,:) = this%uv(:,:,:) + u(1:nx,1:ny,:)*v(1:nx,1:ny,:)*this%dt     !! uv grid
this%uw(:,:,:) = this%uw(:,:,:) + u_w(1:nx,1:ny,:)*w(1:nx,1:ny,:)*this%dt   !! w grid
this%vw(:,:,:) = this%vw(:,:,:) + v_w(1:nx,1:ny,:)*w(1:nx,1:ny,:)*this%dt   !! w grid

this%txx(:,:,:) = this%txx(:,:,:) + txx(1:nx,1:ny,:)*this%dt    !! uv grid
this%tyy(:,:,:) = this%tyy(:,:,:) + tyy(1:nx,1:ny,:)*this%dt    !! uv grid
this%tzz(:,:,:) = this%tzz(:,:,:) + tzz(1:nx,1:ny,:)*this%dt    !! uv grid
this%txy(:,:,:) = this%txy(:,:,:) + txy(1:nx,1:ny,:)*this%dt    !! uv grid
this%txz(:,:,:) = this%txz(:,:,:) + txz(1:nx,1:ny,:)*this%dt    !! w grid
this%tyz(:,:,:) = this%tyz(:,:,:) + tyz(1:nx,1:ny,:)*this%dt    !! w grid

this%eqmxz(:,:) = this%eqmxz(:,:) + eqmxz(1:nx,1:ny)*this%dt    !! 2D Data at 1st grid point (w grid)
this%wpmxz(:,:) = this%wpmxz(:,:) + s_wpmxz(1:nx,1:ny)*this%dt       !! 2D Data at 1st grid point (w grid)
this%eqmyz(:,:) = this%eqmyz(:,:) + eqmyz(1:nx,1:ny)*this%dt    !! 2D Data at 1st grid point (w grid)
this%wpmyz(:,:) = this%wpmyz(:,:) + s_wpmyz(1:nx,1:ny)*this%dt

this%pres_real(:,:,:) = this%pres_real(:,:,:) + pres_real(1:nx,1:ny,:)*this%dt
this%p(:,:,:) = this%p(:,:,:) + p(1:nx,1:ny,:)*this%dt
!this%dpdz(:,:,:) = this%dpdz(:,:,:) + dpdz(1:nx,1:ny,:)*this%dt
!this%dpdz_uv(:,:,:) = this%dpdz_uv(:,:,:) + dpdz_uv(1:nx,1:ny,:)*this%dt

this%cc_z(:,:,:) = this%cc_z(:,:,:) + cc_z(1:nx,1:ny,:)*this%dt
!this%dpdx(:,:,:) = this%dpdx(:,:,:) + dpdx(1:nx,1:ny,:)*this%dt
this%divtx(:,:,:) = this%divtx(:,:,:) + divtx(1:nx,1:ny,:)*this%dt
this%RHSx(:,:,:) = this%RHSx(:,:,:) + RHSx(1:nx,1:ny,:)*this%dt
this%RHSx_f(:,:,:) = this%RHSx_f(:,:,:) + RHSx_f(1:nx,1:ny,:)*this%dt
this%v_omegaz(:,:,:) = this%v_omegaz(:,:,:) + v_omegaz(1:nx,1:ny,:)*this%dt
this%u_omegay(:,:,:) = this%u_omegay(:,:,:) + u_omegay(1:nx,1:ny,:)*this%dt
this%dtxdx(:,:,:) = this%dtxdx(:,:,:) + dtxdx(1:nx,1:ny,:)*this%dt
this%dtydy(:,:,:) = this%dtydy(:,:,:) + dtydy(1:nx,1:ny,:)*this%dt
this%dtzdz(:,:,:) = this%dtzdz(:,:,:) + dtzdz(1:nx,1:ny,:)*this%dt


this%dudy(:,:,:) = this%dudy(:,:,:) + dudy(1:nx,1:ny,:)*this%dt
this%dvdx(:,:,:) = this%dvdx(:,:,:) + dvdx(1:nx,1:ny,:)*this%dt
this%dudz(:,:,:) = this%dudz(:,:,:) + dudz(1:nx,1:ny,:)*this%dt
this%dwdx(:,:,:) = this%dwdx(:,:,:) + dwdx(1:nx,1:ny,:)*this%dt
this%vorty2(:,:,:) = this%vorty2(:,:,:) + vorty2(1:nx,1:ny,:)*this%dt
this%vortz2(:,:,:) = this%vortz2(:,:,:) + vortz2(1:nx,1:ny,:)*this%dt
this%v_omegaz2(:,:,:) = this%v_omegaz2(:,:,:) + v_omegaz2(1:nx,1:ny,:)*this%dt
this%w_omegay2(:,:,:) = this%w_omegay2(:,:,:) + w_omegay2(1:nx,1:ny,:)*this%dt

this%vdudy(:,:,:) = this%vdudy(:,:,:) + vdudy(1:nx,1:ny,:)*this%dt
this%vdvdx(:,:,:) = this%vdvdx(:,:,:) + vdvdx(1:nx,1:ny,:)*this%dt
this%wdudz(:,:,:) = this%wdudz(:,:,:) + wdudz(1:nx,1:ny,:)*this%dt
this%wdwdx(:,:,:) = this%wdwdx(:,:,:) + wdwdx(1:nx,1:ny,:)*this%dt
this%udvdy(:,:,:) = this%udvdy(:,:,:) + udvdy(1:nx,1:ny,:)*this%dt
this%udwdz(:,:,:) = this%udwdz(:,:,:) + udwdz(1:nx,1:ny,:)*this%dt
this%dwdz(:,:,:) = this%dwdz(:,:,:) + dwdz(1:nx,1:ny,:)*this%dt

this%pdetadx(:,:,:) = this%pdetadx(:,:,:) + pdetadx(1:nx,1:ny,:)*this%dt
this%prealdetadx(:,:,:) = this%prealdetadx(:,:,:) + prealdetadx(1:nx,1:ny,:)*this%dt

this%pwalldetadx1(:,:) = this%pwalldetadx1(:,:) + pwalldetadx1(1:nx,1:ny)*this%dt
this%pwalldetadx2(:,:) = this%pwalldetadx2(:,:) + pwalldetadx2(1:nx,1:ny)*this%dt

!this%B_uu(:,:,:) = this%B_uu(:,:,:) + B_uu(1:nx,1:ny,:)*this%dt
!this%B_uu1(:,:,:) = this%B_uu1(:,:,:) + B_uu1(1:nx,1:ny,:)*this%dt
!this%E_uu(:,:,:) = this%E_uu(:,:,:) + E_uu(1:nx,1:ny,:)*this%dt

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
this%fx(:,:,1:) = this%fx(:,:,1:) + fxa(1:nx,1:ny,1:)*this%dt
this%fy(:,:,1:) = this%fy(:,:,1:) + fya(1:nx,1:ny,1:)*this%dt
this%fz(:,:,1:) = this%fz(:,:,1:) + fza_uv(1:nx,1:ny,1:)*this%dt
#endif

this%cs_opt2(:,:,1:) = this%cs_opt2(:,:,1:) + Cs_opt2(1:nx,1:ny,1:)*this%dt

this%vortx(:,:,:) = this%vortx(:,:,:) + vortx(1:nx,1:ny,:) * this%dt
this%vorty(:,:,:) = this%vorty(:,:,:) + vorty(1:nx,1:ny,:) * this%dt
this%vortz(:,:,:) = this%vortz(:,:,:) + vortz(1:nx,1:ny,:) * this%dt

! Update this%total_time for variable time stepping
this%total_time = this%total_time + this%dt

! Set this%dt back to zero for next increment
this%dt = 0._rprec

end subroutine compute

!*******************************************************************************
subroutine finalize(this)
!*******************************************************************************
use grid_m
use param, only : write_endian, lbz, path, coord, nproc
use string_util
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array,MPI_SYNC_DOWNUP
use param, only : ierr,comm
#endif
implicit none

class(tavg_t), intent(inout) :: this

#ifndef PPCGNS
character(64) :: bin_ext
#endif

character(64) :: fname_vel, fname_velw, fname_vel2, fname_tau, fname_pres

character(64) :: fname_mosd, fname_mombal, fname_Corrfunc

character(64) :: fname_f, fname_rs, fname_cs, fname_vort

! integer :: i,j,k
! Where to end with nz index.
integer :: nz_end

real(rprec), pointer, dimension(:) :: x, y, z, zw

real(rprec), allocatable, dimension(:,:,:) :: up2, vp2, wp2, upvp, upwp, vpwp

nullify(x,y,z,zw)

x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

#ifdef PPMPI
! This adds one more element to the last processor (which contains an extra one)
! Processor nproc-1 has data from 1:nz
! Rest of processors have data from 1:nz-1
if ( coord == nproc-1 ) then
    nz_end = 0
else
    nz_end = 1
end if
#else
nz_end = 0
#endif

! Common file name
fname_vel = path // 'output/veluv_avg'
fname_velw = path // 'output/velw_avg'
fname_vel2 = path // 'output/vel2_avg'
fname_tau = path // 'output/tau_avg'

fname_mosd = path // 'output/mosd_avg'
fname_mombal = path // 'output/mombal_avg'
!fname_Corrfunc = path // 'output/Corrfunc_avg'

fname_f = path // 'output/force_avg'
fname_pres = path // 'output/pres_avg'
fname_rs = path // 'output/rs'
fname_cs = path // 'output/cs_opt2'
fname_vort = path // 'output/vort_avg'

! CGNS
#ifdef PPCGNS
call string_concat(fname_vel, '.cgns')
call string_concat(fname_velw, '.cgns')
call string_concat(fname_vel2, '.cgns')
call string_concat(fname_tau, '.cgns')
call string_concat(fname_pres, '.cgns')
call string_concat(fname_f, '.cgns')
call string_concat(fname_rs, '.cgns')
call string_concat(fname_cs, '.cgns')
call string_concat(fname_vort, '.cgns')

! Binary
#else
#ifdef PPMPI
call string_splice(bin_ext, '.c', coord, '.bin')
#else
bin_ext = '.bin'
#endif
call string_concat(fname_vel, bin_ext)
call string_concat(fname_velw, bin_ext)
call string_concat(fname_vel2, bin_ext)
call string_concat(fname_tau, bin_ext)

call string_concat(fname_mosd, bin_ext)
call string_concat(fname_mombal, bin_ext)
!call string_concat(fname_Corrfunc, bin_ext)

call string_concat(fname_pres, bin_ext)
call string_concat(fname_f, bin_ext)
call string_concat(fname_rs, bin_ext)
call string_concat(fname_cs, bin_ext)
call string_concat(fname_vort, bin_ext)
#endif

! Final checkpoint all restart data
call this%checkpoint()

#ifdef PPMPI
call mpi_barrier(comm, ierr)
#endif

!  Perform time averaging operation
this%u(:,:,:) = this%u(:,:,:) /  this%total_time
this%v(:,:,:) = this%v(:,:,:) /  this%total_time
this%w(:,:,:) = this%w(:,:,:) /  this%total_time
this%u_w(:,:,:)  = this%u_w(:,:,:)  /  this%total_time
this%v_w(:,:,:)  = this%v_w(:,:,:)  /  this%total_time
this%w_uv(:,:,:) = this%w_uv(:,:,:) /  this%total_time
this%u2(:,:,:) = this%u2(:,:,:) /  this%total_time
this%v2(:,:,:) = this%v2(:,:,:) /  this%total_time
this%w2(:,:,:) = this%w2(:,:,:) /  this%total_time
this%uv(:,:,:) = this%uv(:,:,:) /  this%total_time
this%uw(:,:,:) = this%uw(:,:,:) /  this%total_time
this%vw(:,:,:) = this%vw(:,:,:) /  this%total_time
this%txx(:,:,:) = this%txx(:,:,:) /  this%total_time
this%tyy(:,:,:) = this%tyy(:,:,:) /  this%total_time
this%tzz(:,:,:) = this%tzz(:,:,:) /  this%total_time
this%txy(:,:,:) = this%txy(:,:,:) /  this%total_time
this%txz(:,:,:) = this%txz(:,:,:) /  this%total_time
this%tyz(:,:,:) = this%tyz(:,:,:) /  this%total_time

this%eqmxz(:,:) = this%eqmxz(:,:) /  this%total_time
this%wpmxz(:,:) = this%wpmxz(:,:) /  this%total_time
this%eqmyz(:,:) = this%eqmyz(:,:) /  this%total_time
this%wpmyz(:,:) = this%wpmyz(:,:) /  this%total_time

this%p(:,:,:) = this%p(:,:,:) /  this%total_time
this%pres_real(:,:,:) = this%pres_real(:,:,:) /  this%total_time
!this%dpdz(:,:,:) = this%dpdz(:,:,:) / this%total_time
!this%dpdz_uv(:,:,:) = this%dpdz_uv(:,:,:) / this%total_time
this%fx(:,:,:) = this%fx(:,:,:) /  this%total_time
this%fy(:,:,:) = this%fy(:,:,:) /  this%total_time
this%fz(:,:,:) = this%fz(:,:,:) /  this%total_time
this%cs_opt2(:,:,:) = this%cs_opt2(:,:,:) /  this%total_time
this%vortx(:,:,:) = this%vortx(:,:,:) /  this%total_time
this%vorty(:,:,:) = this%vorty(:,:,:) /  this%total_time
this%vortz(:,:,:) = this%vortz(:,:,:) /  this%total_time

this%cc_z(:,:,:) = this%cc_z(:,:,:) /  this%total_time
!this%dpdx(:,:,:) = this%dpdx(:,:,:) /  this%total_time
this%divtx(:,:,:) = this%divtx(:,:,:) / this%total_time
this%RHSx(:,:,:) = this%RHSx(:,:,:) / this%total_time
this%RHSx_f(:,:,:) = this%RHSx_f(:,:,:) /  this%total_time
this%v_omegaz(:,:,:) = this%v_omegaz(:,:,:) /  this%total_time
this%u_omegay(:,:,:) = this%u_omegay(:,:,:) /  this%total_time
this%dtxdx(:,:,:) = this%dtxdx(:,:,:) /  this%total_time
this%dtydy(:,:,:) = this%dtydy(:,:,:) /  this%total_time
this%dtzdz(:,:,:) = this%dtzdz(:,:,:) /  this%total_time

this%dudy(:,:,:) = this%dudy(:,:,:) /  this%total_time
this%dvdx(:,:,:) = this%dvdx(:,:,:) /  this%total_time
this%dudz(:,:,:) = this%dudz(:,:,:) /  this%total_time
this%dwdx(:,:,:) = this%dwdx(:,:,:) /  this%total_time
this%vorty2(:,:,:) = this%vorty2(:,:,:) /  this%total_time
this%vortz2(:,:,:) = this%vortz2(:,:,:) /  this%total_time
this%v_omegaz2(:,:,:) = this%v_omegaz2(:,:,:) /  this%total_time
this%w_omegay2(:,:,:) = this%w_omegay2(:,:,:) /  this%total_time

this%vdudy(:,:,:) = this%vdudy(:,:,:) /  this%total_time
this%vdvdx(:,:,:) = this%vdvdx(:,:,:) /  this%total_time
this%wdudz(:,:,:) = this%wdudz(:,:,:) /  this%total_time
this%wdwdx(:,:,:) = this%wdwdx(:,:,:) /  this%total_time
this%udvdy(:,:,:) = this%udvdy(:,:,:) /  this%total_time
this%udwdz(:,:,:) = this%udwdz(:,:,:) /  this%total_time

this%pdetadx(:,:,:) = this%pdetadx(:,:,:) /  this%total_time
this%prealdetadx(:,:,:) = this%prealdetadx(:,:,:) /  this%total_time
this%pwalldetadx1(:,:) = this%pwalldetadx1(:,:) /  this%total_time
this%pwalldetadx2(:,:) = this%pwalldetadx2(:,:) /  this%total_time

!this%B_uu(:,:,:) = this%B_uu(:,:,:) /  this%total_time
!this%B_uu1(:,:,:) = this%B_uu1(:,:,:) /  this%total_time
!this%E_uu(:,:,:) = this%E_uu(:,:,:) /  this%total_time

#ifdef PPMPI
call mpi_barrier( comm, ierr )
#endif

!  Sync entire tavg structure
#ifdef PPMPI
call mpi_sync_real_array( this%u(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%v(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%w(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%u2(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%v2(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%w2(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%uw(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%vw(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%uv(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%pres_real(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array( this%dpdz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%fx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%cs_opt2(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%vortx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%vorty(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%vortz(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%pdetadx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( this%prealdetadx(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array( this%B_uu(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array( this%B_uu1(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
!call mpi_sync_real_array( this%E_uu(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
#endif

! Write all the 3D data
#ifdef PPCGNS
! Write CGNS Data
call write_parallel_cgns (fname_vel ,nx, ny, nz - nz_end, nz_tot,              &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 3,                                  &
    (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                               &
    (/ this%u(1:nx,1:ny,1:nz-nz_end),                                          &
       this%v(1:nx,1:ny,1:nz-nz_end),                                          &
       this%w_uv(1:nx,1:ny,1:nz-nz_end) /) )

call write_parallel_cgns (fname_velw ,nx, ny, nz - nz_end, nz_tot,             &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ),                                    &
    1, (/ 'VelocityZ' /), (/ this%w(1:nx,1:ny,1:nz-nz_end) /) )

call write_parallel_cgns(fname_vel2,nx,ny,nz- nz_end,nz_tot,                   &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 6,                                 &
    (/ 'Mean--uu', 'Mean--vv', 'Mean--ww','Mean--uw','Mean--vw','Mean--uv'/),  &
    (/ this%u2(1:nx,1:ny,1:nz-nz_end),                                         &
       this%v2(1:nx,1:ny,1:nz-nz_end),                                         &
       this%w2(1:nx,1:ny,1:nz-nz_end),                                         &
       this%uw(1:nx,1:ny,1:nz-nz_end),                                         &
       this%vw(1:nx,1:ny,1:nz-nz_end),                                         &
       this%uv(1:nx,1:ny,1:nz-nz_end) /) )

call write_parallel_cgns(fname_tau,nx,ny,nz- nz_end,nz_tot,                    &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 6,                                 &
    (/ 'Tau--txx', 'Tau--txy', 'Tau--tyy','Tau--txz','Tau--tyz','Tau--tzz'/),  &
    (/ this%txx(1:nx,1:ny,1:nz-nz_end),                                        &
       this%txy(1:nx,1:ny,1:nz-nz_end),                                        &
       this%tyy(1:nx,1:ny,1:nz-nz_end),                                        &
       this%txz(1:nx,1:ny,1:nz-nz_end),                                        &
       this%tyz(1:nx,1:ny,1:nz-nz_end),                                        &
       this%tzz(1:nx,1:ny,1:nz-nz_end) /) )

call write_parallel_cgns(fname_pres,nx,ny,nz- nz_end,nz_tot,                   &
   (/ 1, 1,   (nz-1)*coord + 1 /),                                             &
   (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                                &
   x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 1,                                  &
   (/ 'pressure' /),                                                           &
   (/ this%p(1:nx,1:ny,1:nz-nz_end) /) )

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
call write_parallel_cgns(fname_f,nx,ny,nz- nz_end,nz_tot,                      &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 3,                                 &
    (/ 'bodyForX', 'bodyForY', 'bodyForZ' /),                                  &
    (/ this%fx(1:nx,1:ny,1:nz-nz_end),                                         &
       this%fy(1:nx,1:ny,1:nz-nz_end),                                         &
       this%fz(1:nx,1:ny,1:nz-nz_end) /) )
#endif

call write_parallel_cgns(fname_cs,nx,ny,nz- nz_end,nz_tot,                     &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 1,                                 &
    (/ 'Cs_Coeff'/),  (/ this%cs_opt2(1:nx,1:ny,1:nz- nz_end) /) )

call write_parallel_cgns(fname_vort,nx,ny,nz- nz_end,nz_tot,                   &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ), 3,                                 &
    (/ 'VorticityX', 'VorticityY', 'VorticityZ' /),                            &
    (/ this%vortx(1:nx,1:ny,1:nz-nz_end),                                      &
       this%vorty(1:nx,1:ny,1:nz-nz_end),                                      &
       this%vortz(1:nx,1:ny,1:nz-nz_end) /)

#else
! Write binary data
open(unit=13, file=fname_vel, form='unformatted', convert=write_endian,        &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%u(:nx,:ny,1:nz)
write(13,rec=2) this%v(:nx,:ny,1:nz)
write(13,rec=3) this%w_uv(:nx,:ny,1:nz)
close(13)

! Write binary data
open(unit=13, file=fname_velw, form='unformatted', convert=write_endian,       &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%w(:nx,:ny,1:nz)
close(13)

open(unit=13, file=fname_vel2, form='unformatted', convert=write_endian,       &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%u2(:nx,:ny,1:nz)
write(13,rec=2) this%v2(:nx,:ny,1:nz)
write(13,rec=3) this%w2(:nx,:ny,1:nz)
write(13,rec=4) this%uw(:nx,:ny,1:nz)
write(13,rec=5) this%vw(:nx,:ny,1:nz)
write(13,rec=6) this%uv(:nx,:ny,1:nz)
close(13)

open(unit=13, file=fname_tau, form='unformatted', convert=write_endian,        &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%txx(:nx,:ny,1:nz)
write(13,rec=2) this%txy(:nx,:ny,1:nz)
write(13,rec=3) this%tyy(:nx,:ny,1:nz)
write(13,rec=4) this%txz(:nx,:ny,1:nz)
write(13,rec=5) this%tyz(:nx,:ny,1:nz)
write(13,rec=6) this%tzz(:nx,:ny,1:nz)
close(13)

if (coord == 0) then
open(unit=13, file=fname_mosd, form='unformatted', convert=write_endian,        &
    access='direct', recl=nx*ny*1*rprec)
write(13,rec=1) this%eqmxz(:nx,:ny)
write(13,rec=2) this%wpmxz(:nx,:ny)
write(13,rec=3) this%eqmyz(:nx,:ny)
write(13,rec=4) this%wpmyz(:nx,:ny)
write(13,rec=5) this%pwalldetadx1(:nx,:ny)
write(13,rec=6) this%pwalldetadx2(:nx,:ny)
close(13)
end if

open(unit=13, file=fname_pres, form='unformatted', convert=write_endian,       &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%p(:nx,:ny,1:nz)
write(13,rec=2) this%pres_real(:nx,:ny,1:nz)
write(13,rec=3) this%pdetadx(:nx,:ny,1:nz)
write(13,rec=4) this%prealdetadx(:nx,:ny,1:nz)
close(13)

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
open(unit=13, file=fname_f, form='unformatted', convert=write_endian,          &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%fx(:nx,:ny,1:nz)
write(13,rec=2) this%fy(:nx,:ny,1:nz)
write(13,rec=3) this%fz(:nx,:ny,1:nz)
close(13)
#endif

open(unit=13, file=fname_cs, form='unformatted', convert=write_endian,         &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%cs_opt2(:nx,:ny,1:nz)
close(13)

open(unit=13, file=fname_vort, form='unformatted', convert=write_endian,       &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%vortx(:nx,:ny,1:nz)
write(13,rec=2) this%vorty(:nx,:ny,1:nz)
write(13,rec=3) this%vortz(:nx,:ny,1:nz)
close(13)

open(unit=13, file=fname_mombal, form='unformatted', convert=write_endian,       &
    access='direct', recl=nx*ny*nz*rprec)
write(13,rec=1) this%cc_z(:nx,:ny,1:nz)
!write(13,rec=2) this%dpdx(:nx,:ny,1:nz)
write(13,rec=3) this%divtx(:nx,:ny,1:nz)
write(13,rec=4) this%RHSx(:nx,:ny,1:nz)
write(13,rec=5) this%RHSx_f(:nx,:ny,1:nz)
write(13,rec=6) this%v_omegaz(:nx,:ny,1:nz)
write(13,rec=7) this%u_omegay(:nx,:ny,1:nz)
write(13,rec=8) this%dtxdx(:nx,:ny,1:nz)
write(13,rec=9) this%dtydy(:nx,:ny,1:nz)
write(13,rec=10) this%dtzdz(:nx,:ny,1:nz)

write(13,rec=11) this%dudy(:nx,:ny,1:nz)
write(13,rec=12) this%dvdx(:nx,:ny,1:nz)
write(13,rec=13) this%dudz(:nx,:ny,1:nz)
write(13,rec=14) this%dwdx(:nx,:ny,1:nz)
write(13,rec=15) this%vorty2(:nx,:ny,1:nz)
write(13,rec=16) this%vortz2(:nx,:ny,1:nz)
write(13,rec=17) this%v_omegaz2(:nx,:ny,1:nz)
write(13,rec=18) this%w_omegay2(:nx,:ny,1:nz)

write(13,rec=19) this%vdudy(:nx,:ny,1:nz)
write(13,rec=20) this%vdvdx(:nx,:ny,1:nz)
write(13,rec=21) this%wdudz(:nx,:ny,1:nz)
write(13,rec=22) this%wdwdx(:nx,:ny,1:nz)
write(13,rec=23) this%udvdy(:nx,:ny,1:nz)
write(13,rec=24) this%udwdz(:nx,:ny,1:nz)
write(13,rec=25) this%dwdz(:nx,:ny,1:nz)
close(13)

!open(unit=13, file=fname_Corrfunc, form='unformatted', convert=write_endian,         &
!    access='direct', recl=nx*ny*nz*rprec)
!write(13,rec=1) this%B_uu(:nx,:ny,1:nz)
!write(13,rec=2) this%B_uu1(:nx,:ny,1:nz)
!write(13,rec=3) this%E_uu(:nx,:ny,1:nz)
!close(13)

#endif

#ifdef PPMPI
! Ensure all writes complete before preceeding
call mpi_barrier( comm, ierr )
#endif

! Do the Reynolds stress calculations afterwards. Now we can interpolate w and
! ww to the uv grid and do the calculations. We have already written the data to
! the files so we can overwrite now
allocate( up2(nx,ny,lbz:nz) )
allocate( vp2(nx,ny,lbz:nz) )
allocate( wp2(nx,ny,lbz:nz) )
allocate( upvp(nx,ny,lbz:nz) )
allocate( upwp(nx,ny,lbz:nz) )
allocate( vpwp(nx,ny,lbz:nz) )
up2 = this%u2 - this%u * this%u
vp2 = this%v2 - this%v * this%v
wp2 = this%w2 - this%w * this%w
upvp = this%uv - this%u * this%v
!! using u_w and v_w below instead of u and v ensures that the Reynolds
!! stresses are on the same grid as the squared velocities (i.e., w-grid)
upwp = this%uw - this%u_w * this%w
vpwp = this%vw - this%v_w * this%w

#ifdef PPCGNS
! Write CGNS data
call write_parallel_cgns(fname_rs,nx,ny,nz- nz_end,nz_tot,                     &
    (/ 1, 1,   (nz-1)*coord + 1 /),                                            &
    (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                               &
    x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ), 6,                                  &
    (/ 'Meanupup', 'Meanvpvp', 'Meanwpwp','Meanupwp','Meanvpwp','Meanupvp'/),  &
    (/ up2(1:nx,1:ny,1:nz- nz_end) ,                                         &
    vp2(1:nx,1:ny,1:nz- nz_end) ,                                              &
    wp2(1:nx,1:ny,1:nz- nz_end) ,                                              &
    upwp(1:nx,1:ny,1:nz- nz_end) ,                                             &
    vpwp(1:nx,1:ny,1:nz- nz_end) ,                                             &
    upvp(1:nx,1:ny,1:nz- nz_end)  /) )
#else
! Write binary data
open(unit=13, file=fname_rs, form='unformatted', convert=write_endian,         &
    access='direct',recl=nx*ny*nz*rprec)
write(13,rec=1) up2(:nx,:ny,1:nz)
write(13,rec=2) vp2(:nx,:ny,1:nz)
write(13,rec=3) wp2(:nx,:ny,1:nz)
write(13,rec=4) upwp(:nx,:ny,1:nz)
write(13,rec=5) vpwp(:nx,:ny,1:nz)
write(13,rec=6) upvp(:nx,:ny,1:nz)
close(13)
#endif

#ifdef PPMPI
! Ensure all writes complete before preceeding
call mpi_barrier( comm, ierr )
#endif

end subroutine finalize

!*******************************************************************************
subroutine checkpoint(this)
!*******************************************************************************
!
! This subroutine writes the restart data and is to be called by 'checkpoint'
! for intermediate checkpoints and by 'tavg_finalize' at the end of the
! simulation.
!
use param, only : write_endian, coord
use string_util
implicit none

class(tavg_t), intent(inout) :: this

character(64) :: fname

fname = checkpoint_tavg_file
#ifdef PPMPI
call string_concat( fname, '.c', coord)
#endif

!  Write data to tavg.out
open(1, file=fname, action='write', position='rewind',form='unformatted',      &
    convert=write_endian)
write(1) this%total_time
write(1) this%u
write(1) this%v
write(1) this%w
write(1) this%u_w
write(1) this%v_w
write(1) this%w_uv
write(1) this%u2
write(1) this%v2
write(1) this%w2
write(1) this%uv
write(1) this%uw
write(1) this%vw
write(1) this%txx
write(1) this%tyy
write(1) this%tzz
write(1) this%txy
write(1) this%txz
write(1) this%tyz
!write(1) this%dpdz
!write(1) this%dpdz_uv
write(1) this%p
write(1) this%pres_real

write(1) this%eqmxz
write(1) this%wpmxz
write(1) this%eqmyz
write(1) this%wpmyz

write(1) this%fx
write(1) this%fy
write(1) this%fz
write(1) this%cs_opt2
write(1) this%vortx
write(1) this%vorty
write(1) this%vortz

    write(1) this%cc_z
!    write(1) this%dpdx
    write(1) this%divtx
    write(1) this%RHSx
    write(1) this%RHSx_f
    write(1) this%v_omegaz
    write(1) this%u_omegay
    write(1) this%dtxdx
    write(1) this%dtydy
    write(1) this%dtzdz

    write(1) this%v_omegaz2
    write(1) this%w_omegay2
    write(1) this%dudz
    write(1) this%dvdx
    write(1) this%dudy
    write(1) this%dwdx
    write(1) this%vorty2
    write(1) this%vortz2

    write(1) this%wdudz
    write(1) this%vdvdx
    write(1) this%vdudy
    write(1) this%wdwdx
    write(1) this%udvdy
    write(1) this%udwdz
    write(1) this%dwdz

    write(1) this%pdetadx
    write(1) this%prealdetadx
    write(1) this%pwalldetadx1
    write(1) this%pwalldetadx2
  
!    write(1) this%B_uu
!    write(1) this%B_uu1
!    write(1) this%E_uu

close(1)

end subroutine checkpoint

#ifdef PPCGNS
#ifdef PPMPI
!*******************************************************************************
subroutine write_parallel_cgns (file_name, nx, ny, nz, nz_tot, start_n_in,     &
    end_n_in, xin, yin, zin, num_fields, fieldNames, input )
!*******************************************************************************
use param, only : coord
implicit none

integer, intent(in) :: nx, ny, nz, nz_tot, num_fields
! Name of file to be written
character(*), intent(in) :: file_name
! Name of fields we are writing
character(*), intent(in), dimension(:) :: fieldNames
! Data to be written
real(rprec), intent(in), dimension(:) :: input
! Coordinates to write
real(rprec), intent(in), dimension(:) :: xin, yin, zin
! Where the total node counter starts nodes
integer, intent(in) :: start_n_in(3)
! Where the total node counter ends nodes
integer, intent(in) :: end_n_in(3)

integer :: fn=1        ! CGNS file index number
integer :: ier         ! CGNS error status
integer :: base=1      ! base number
integer :: zone=1      ! zone number
integer :: nnodes      ! Number of nodes in this processor
integer :: sol =1      ! solution number
integer :: field       ! section number
integer(cgsize_t) :: sizes(3,3)  ! Sizes

! Convert input to right data type
integer(cgsize_t) :: start_n(3)  ! Where the total node counter starts nodes
integer(cgsize_t) :: end_n(3)  ! Where the total node counter ends nodes

! Building the lcoal mesh
integer :: i,j,k
real(rprec), dimension(nx,ny,nz) :: xyz

!  ! Set the parallel communicator
!  call cgp_mpi_comm_f(cgnsParallelComm, ierr)

! Convert types such that CGNS libraries can handle the input
start_n(1) = int(start_n_in(1), cgsize_t)
start_n(2) = int(start_n_in(2), cgsize_t)
start_n(3) = int(start_n_in(3), cgsize_t)
end_n(1) = int(end_n_in(1), cgsize_t)
end_n(2) = int(end_n_in(2), cgsize_t)
end_n(3) = int(end_n_in(3), cgsize_t)

! The total number of nodes in this processor
nnodes = nx*ny*nz

! Sizes, used to create zone
sizes(:,1) = (/int(nx, cgsize_t),int(ny, cgsize_t),int(nz_tot, cgsize_t)/)
sizes(:,2) = (/int(nx-1, cgsize_t),int(ny-1, cgsize_t),int(nz_tot-1, cgsize_t)/)
sizes(:,3) = (/int(0, cgsize_t) , int(0, cgsize_t), int(0, cgsize_t)/)

! Open CGNS file
call cgp_open_f(file_name, CG_MODE_WRITE, fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write base
call cg_base_write_f(fn, 'Base', 3, 3, base, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write zone
call cg_zone_write_f(fn, base, 'Zone', sizes, Structured, zone, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write print info to screen
if (coord .eq. 0) then
    write(*,*) 'Writing, ', file_name
end if

! Create data nodes for coordinates
call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateX', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateY', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateZ', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

! This is done for the 3 dimensions x,y and z
! It writes the coordinates
! Create grid points
do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = xin(i)
end do
end do
end do

call cgp_coord_write_data_f(fn, base, zone, 1,                                 &
    start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = yin(j)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 2,   &
    start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = zin(k)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 3,   &
                            start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Create a centered solution
call cg_sol_write_f(fn, base, zone, 'Solution', Vertex, sol, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the solution
do i=1,num_fields
    call cgp_field_write_f(fn, base, zone, sol, RealDouble, fieldNames(i),     &
        field, ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

    call cgp_field_write_data_f(fn, base, zone, sol, field, start_n, end_n,    &
        input((i-1)*nnodes+1:(i)*nnodes), ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

end do

! Close the file
call cgp_close_f(fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

end subroutine write_parallel_cgns

!*******************************************************************************
subroutine write_null_cgns (file_name, nx, ny, nz, nz_tot, start_n_in,         &
    end_n_in, xin, yin, zin, num_fields, fieldNames )
!*******************************************************************************
use param, only : coord
implicit none

integer, intent(in) :: nx, ny, nz, nz_tot, num_fields
! Name of file to be written
character(*), intent(in) :: file_name
! Name of fields we are writing
character(*), intent(in), dimension(:) :: fieldNames
! Coordinates to write
real(rprec), intent(in), dimension(:) :: xin, yin, zin
! Where the total node counter starts nodes
integer, intent(in) :: start_n_in(3)
! Where the total node counter ends nodes
integer, intent(in) :: end_n_in(3)

integer :: fn=1        ! CGNS file index number
integer :: ier         ! CGNS error status
integer :: base=1      ! base number
integer :: zone=1      ! zone number
integer :: nnodes      ! Number of nodes in this processor
integer :: sol =1      ! solution number
integer :: field       ! section number
integer(cgsize_t) :: sizes(3,3)  ! Sizes

! Convert input to right data type
integer(cgsize_t) :: start_n(3)  ! Where the total node counter starts nodes
integer(cgsize_t) :: end_n(3)  ! Where the total node counter ends nodes

! Building the lcoal mesh
integer :: i,j,k
real(rprec), dimension(nx,ny,nz) :: xyz

!  ! Set the parallel communicator
!  call cgp_mpi_comm_f(cgnsParallelComm, ierr)

! Convert types such that CGNS libraries can handle the input
start_n(1) = int(start_n_in(1), cgsize_t)
start_n(2) = int(start_n_in(2), cgsize_t)
start_n(3) = int(start_n_in(3), cgsize_t)
end_n(1) = int(end_n_in(1), cgsize_t)
end_n(2) = int(end_n_in(2), cgsize_t)
end_n(3) = int(end_n_in(3), cgsize_t)

! The total number of nodes in this processor
nnodes = nx*ny*nz

! Sizes, used to create zone
sizes(:,1) = (/int(nx, cgsize_t),int(ny, cgsize_t),int(nz_tot, cgsize_t)/)
sizes(:,2) = (/int(nx-1, cgsize_t),int(ny-1, cgsize_t),int(nz_tot-1, cgsize_t)/)
sizes(:,3) = (/int(0, cgsize_t) , int(0, cgsize_t), int(0, cgsize_t)/)

! Open CGNS file
call cgp_open_f(file_name, CG_MODE_WRITE, fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write base
call cg_base_write_f(fn, 'Base', 3, 3, base, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write zone
call cg_zone_write_f(fn, base, 'Zone', sizes, Structured, zone, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write print info to screen
if (coord .eq. 0) then
    write(*,*) 'Writing, ', file_name
end if

! Create data nodes for coordinates
call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateX', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateY', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateZ', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! This is done for the 3 dimensions x,y and z
! It writes the coordinates
! Create grid points
do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = xin(i)
end do
end do
end do

call cgp_coord_write_data_f(fn, base, zone, 1, start_n, end_n, %VAL(0), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = yin(j)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 2, start_n, end_n, %VAL(0), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = zin(k)
end do
end do
end do

call cgp_coord_write_data_f(fn, base, zone, 3, start_n, end_n, %VAL(0), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Create a centered solution
call cg_sol_write_f(fn, base, zone, 'Solution', Vertex, sol, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the solution
do i = 1, num_fields
    call cgp_field_write_f(fn, base, zone, sol, RealDouble, fieldNames(i),     &
                           field, ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

    call cgp_field_write_data_f(fn, base, zone, sol, field, start_n, end_n,    &
                                %VAL(0), ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

end do

! Close the file
call cgp_close_f(fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

write(*,*) "end of write_null_cgns"

end subroutine write_null_cgns
#endif
#endif

end module time_average
