!!
!!  Copyright (C) 2010-2016  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!*******************************************************************************
module phase_average
!*******************************************************************************
use types, only : rprec
use param, only : nx, ny, nz, lbz

private
public :: pavg_t


type pavg_t
    real(rprec), dimension(:,:,:,:), allocatable :: u, v, w, u_w, v_w, w_uv, p, &
            pres_real
    real(rprec), dimension(:,:,:,:), allocatable :: uu, vv, ww, uv, uw, vw, uw_uv, vw_uv, ww_uv
    real(rprec), dimension(:,:,:,:), allocatable :: txx, tyy, tzz, txy, txz, tyz
    real(rprec), dimension(:,:,:,:), allocatable :: dudx, dudy, dudz, &
       dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, dpdz, dpdz_uv
    real(rprec), dimension(:,:,:,:), allocatable ::  fx, fy, fz
    real(rprec), dimension(:), allocatable :: np
    real(rprec) :: time
    logical :: initialized = .false.
contains
    procedure, public :: init
    procedure, public :: compute
    procedure, public :: finalize
    procedure, public :: checkpoint
end type pavg_t

character(:), allocatable :: checkpoint_pavg_file

contains

!*******************************************************************************
subroutine init(this)
!*******************************************************************************
use messages
use string_util
use param, only : read_endian, coord, path, pavg_nbins, lbz
implicit none

class(pavg_t), intent(inout) :: this

character (128) :: fname

logical :: exst
integer :: fid

! Create file name
allocate(checkpoint_pavg_file, source=path // 'pavg.out')

! Allocate and initialize

allocate( this%u(1:nx,1:ny,1:nz,pavg_nbins) ); this%u(:,:,:,:) = 0._rprec
allocate( this%v(1:nx,1:ny,1:nz,pavg_nbins) ); this%v(:,:,:,:) = 0._rprec
allocate( this%w(1:nx,1:ny,1:nz,pavg_nbins) ); this%w(:,:,:,:) = 0._rprec
allocate( this%u_w(1:nx,1:ny,1:nz,pavg_nbins) ); this%u_w(:,:,:,:) = 0._rprec
allocate( this%v_w(1:nx,1:ny,1:nz,pavg_nbins) ); this%v_w(:,:,:,:) = 0._rprec
allocate( this%w_uv(1:nx,1:ny,1:nz,pavg_nbins) ); this%w_uv(:,:,:,:) = 0._rprec

allocate( this%uw_uv(1:nx,1:ny,1:nz,pavg_nbins) ); this%uw_uv(:,:,:,:) = 0._rprec
allocate( this%vw_uv(1:nx,1:ny,1:nz,pavg_nbins) ); this%vw_uv(:,:,:,:) = 0._rprec
allocate( this%ww_uv(1:nx,1:ny,1:nz,pavg_nbins) ); this%ww_uv(:,:,:,:) = 0._rprec

allocate( this%uu(1:nx,1:ny,1:nz,pavg_nbins) ); this%uu(:,:,:,:) = 0._rprec
allocate( this%vv(1:nx,1:ny,1:nz,pavg_nbins) ); this%vv(:,:,:,:) = 0._rprec
allocate( this%ww(1:nx,1:ny,1:nz,pavg_nbins) ); this%ww(:,:,:,:) = 0._rprec
allocate( this%uv(1:nx,1:ny,1:nz,pavg_nbins) ); this%uv(:,:,:,:) = 0._rprec
allocate( this%uw(1:nx,1:ny,1:nz,pavg_nbins) ); this%uw(:,:,:,:) = 0._rprec
allocate( this%vw(1:nx,1:ny,1:nz,pavg_nbins) ); this%vw(:,:,:,:) = 0._rprec
allocate( this%txx(1:nx,1:ny,1:nz,pavg_nbins) ); this%txx(:,:,:,:) = 0._rprec
allocate( this%tyy(1:nx,1:ny,1:nz,pavg_nbins) ); this%tyy(:,:,:,:) = 0._rprec
allocate( this%tzz(1:nx,1:ny,1:nz,pavg_nbins) ); this%tzz(:,:,:,:) = 0._rprec
allocate( this%txy(1:nx,1:ny,1:nz,pavg_nbins) ); this%txy(:,:,:,:) = 0._rprec
allocate( this%txz(1:nx,1:ny,1:nz,pavg_nbins) ); this%txz(:,:,:,:) = 0._rprec
allocate( this%tyz(1:nx,1:ny,1:nz,pavg_nbins) ); this%tyz(:,:,:,:) = 0._rprec

allocate( this%dudx(1:nx,1:ny,1:nz,pavg_nbins) ); this%dudx(:,:,:,:) = 0._rprec
allocate( this%dudy(1:nx,1:ny,1:nz,pavg_nbins) ); this%dudy(:,:,:,:) = 0._rprec
allocate( this%dudz(1:nx,1:ny,1:nz,pavg_nbins) ); this%dudz(:,:,:,:) = 0._rprec
allocate( this%dvdx(1:nx,1:ny,1:nz,pavg_nbins) ); this%dvdx(:,:,:,:) = 0._rprec
allocate( this%dvdy(1:nx,1:ny,1:nz,pavg_nbins) ); this%dvdy(:,:,:,:) = 0._rprec
allocate( this%dvdz(1:nx,1:ny,1:nz,pavg_nbins) ); this%dvdz(:,:,:,:) = 0._rprec
allocate( this%dwdx(1:nx,1:ny,1:nz,pavg_nbins) ); this%dwdx(:,:,:,:) = 0._rprec
allocate( this%dwdy(1:nx,1:ny,1:nz,pavg_nbins) ); this%dwdy(:,:,:,:) = 0._rprec
allocate( this%dwdz(1:nx,1:ny,1:nz,pavg_nbins) ); this%dwdz(:,:,:,:) = 0._rprec
!allocate( this%dpdz(1:nx,1:ny,1:nz,pavg_nbins) ); this%dpdz(:,:,:,:) = 0._rprec
!allocate( this%dpdz_uv(1:nx,1:ny,1:nz,pavg_nbins) ); this%dpdz_uv(:,:,:,:) = 0._rprec

allocate( this%p(1:nx,1:ny,1:nz,pavg_nbins) ); this%p(:,:,:,:) = 0._rprec
allocate( this%pres_real(1:nx,1:ny,1:nz,pavg_nbins) ); this%pres_real(:,:,:,:) = 0._rprec

allocate( this%fx(1:nx,1:ny,1:nz,pavg_nbins) ); this%fx(:,:,:,:) = 0._rprec
allocate( this%fy(1:nx,1:ny,1:nz,pavg_nbins) ); this%fy(:,:,:,:) = 0._rprec
allocate( this%fz(1:nx,1:ny,1:nz,pavg_nbins) ); this%fz(:,:,:,:) = 0._rprec

allocate( this%np(pavg_nbins)); this%np(:) = 0._rprec

fname = checkpoint_pavg_file
call string_concat( fname, '.c', coord)

inquire (file=fname, exist=exst)
if (.not. exst) then
    !  Nothing to read in
    if (coord == 0) then
        write(*,*) ' '
        write(*,*)'No previous phase averaged data - starting from scratch.'
    end if
else
    if (coord == 0) then
        write(*,*)'Reading phase averaging data.'
    end if
    open(newunit=fid, file=fname, action='read', position='rewind', form='unformatted',  &
        convert=read_endian)
    read(fid) this%np
    read(fid) this%u
    read(fid) this%v
    read(fid) this%w
    read(fid) this%u_w
    read(fid) this%v_w
    read(fid) this%w_uv
    
    read(fid) this%uu
    read(fid) this%vv
    read(fid) this%ww
    read(fid) this%uv
    read(fid) this%uw
    read(fid) this%vw
   
    read(fid) this%uw_uv
    read(fid) this%vw_uv
    read(fid) this%ww_uv

    read(fid) this%txx
    read(fid) this%tyy
    read(fid) this%tzz
    read(fid) this%txy
    read(fid) this%txz
    read(fid) this%tyz
    
    read(fid) this%dudx
    read(fid) this%dudy
    read(fid) this%dudz
    read(fid) this%dvdx
    read(fid) this%dvdy
    read(fid) this%dvdz
    read(fid) this%dwdx
    read(fid) this%dwdy
    read(fid) this%dwdz
 !   read(fid) this%dpdz     
 !   read(fid) this%dpdz_uv   
    
    read(fid) this%p
    read(fid) this%pres_real
    
    read(fid) this%fx
    read(fid) this%fy
    read(fid) this%fz


    close(fid)
end if

!allocate(pres_real(nx,ny,lbz:nz))

! Set global switch that pavg as been initialized
this%initialized = .true.

end subroutine init

!*******************************************************************************
subroutine compute(this)
!*******************************************************************************
!
!  This subroutine collects the stats for each flow
!  variable quantity
!
use param, only : dt, ubc_mom, lbc_mom, coord, nproc
use param, only : pi, total_time, jt_total, lbz
use param, only : wave_freq,  pavg_tstart, pavg_nbins
use sim_param, only : u, v, w, p
use sim_param, only : dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, dpdz
use sim_param, only : txx, txy, tyy, txz, tyz, tzz
use sim_param, only : fxa, fya, fza
use functions, only : interp_to_uv_grid, interp_to_w_grid

implicit none

class(pavg_t), intent(inout) :: this
real(rprec), dimension(nx,ny,lbz:nz) :: u_w, v_w, w_uv, pres_real, dpdz_uv
real(rprec) :: Tf
integer :: bin

! Interpolation onto other grids
w_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz )
u_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(u(1:nx,1:ny,lbz:nz), lbz )
v_w(1:nx,1:ny,lbz:nz) = interp_to_w_grid(v(1:nx,1:ny,lbz:nz), lbz )
!dpdz_uv(1:nx,1:ny,lbz:nz) = interp_to_uv_grid(dpdz(1:nx,1:ny,lbz:nz), lbz )


!! note: u_w not necessarily zero on walls, but only mult by w=0 vu u'w', so OK
!! can zero u_w at BC anyway:
!if(coord==0       .and. lbc_mom>0) u_w(:,:,1)  = 0._rprec
!if(coord==nproc-1 .and. ubc_mom>0) u_w(:,:,nz) = 0._rprec
!if(coord==0       .and. lbc_mom>0) v_w(:,:,1)  = 0._rprec
!if(coord==nproc-1 .and. ubc_mom>0) v_w(:,:,nz) = 0._rprec

! Calculate real pressure (not pseudo-pressure)
pres_real(1:nx,1:ny,lbz:nz) = 0._rprec
pres_real(1:nx,1:ny,lbz:nz) = p(1:nx,1:ny,lbz:nz)                              &
    - 0.5 * ( u(1:nx,1:ny,lbz:nz)**2 + w_uv(1:nx,1:ny,lbz:nz)**2               &
    + v(1:nx,1:ny,lbz:nz)**2 )

this%time = total_time - pavg_tstart  
Tf = 2._rprec*pi/wave_freq
bin = ceiling(pavg_nbins/Tf*modulo(this%time,Tf))


call pavg_compute(this%time,u(1:nx,1:ny,1:nz),this%u)
call pavg_compute(this%time,v(1:nx,1:ny,1:nz),this%v)
call pavg_compute(this%time,w(1:nx,1:ny,1:nz),this%w)
call pavg_compute(this%time,u_w(1:nx,1:ny,1:nz),this%u_w)
call pavg_compute(this%time,v_w(1:nx,1:ny,1:nz),this%v_w)
call pavg_compute(this%time,w_uv(1:nx,1:ny,1:nz),this%w_uv)

call pavg_compute(this%time,u(1:nx,1:ny,1:nz)*u(1:nx,1:ny,1:nz),this%uu)
call pavg_compute(this%time,v(1:nx,1:ny,1:nz)*v(1:nx,1:ny,1:nz),this%vv)
call pavg_compute(this%time,w(1:nx,1:ny,1:nz)*w(1:nx,1:ny,1:nz),this%ww)
call pavg_compute(this%time,u(1:nx,1:ny,1:nz)*v(1:nx,1:ny,1:nz),this%uv)
call pavg_compute(this%time,u_w(1:nx,1:ny,1:nz)*w(1:nx,1:ny,1:nz),this%uw)
call pavg_compute(this%time,v_w(1:nx,1:ny,1:nz)*w(1:nx,1:ny,1:nz),this%vw)

call pavg_compute(this%time,u(1:nx,1:ny,1:nz)*w_uv(1:nx,1:ny,1:nz),this%uw_uv)
call pavg_compute(this%time,v(1:nx,1:ny,1:nz)*w_uv(1:nx,1:ny,1:nz),this%vw_uv)
call pavg_compute(this%time,w_uv(1:nx,1:ny,1:nz)*w_uv(1:nx,1:ny,1:nz),this%ww_uv)

call pavg_compute(this%time,txx(1:nx,1:ny,1:nz),this%txx)
call pavg_compute(this%time,tyy(1:nx,1:ny,1:nz),this%tyy)
call pavg_compute(this%time,tzz(1:nx,1:ny,1:nz),this%tzz)
call pavg_compute(this%time,txy(1:nx,1:ny,1:nz),this%txy)
call pavg_compute(this%time,txz(1:nx,1:ny,1:nz),this%txz)
call pavg_compute(this%time,tyz(1:nx,1:ny,1:nz),this%tyz)

call pavg_compute(this%time,dudx(1:nx,1:ny,1:nz),this%dudx)
call pavg_compute(this%time,dudy(1:nx,1:ny,1:nz),this%dudy)
call pavg_compute(this%time,dudz(1:nx,1:ny,1:nz),this%dudz)
call pavg_compute(this%time,dvdx(1:nx,1:ny,1:nz),this%dvdx)
call pavg_compute(this%time,dvdy(1:nx,1:ny,1:nz),this%dvdy)
call pavg_compute(this%time,dvdz(1:nx,1:ny,1:nz),this%dvdz)
call pavg_compute(this%time,dwdx(1:nx,1:ny,1:nz),this%dwdx)
call pavg_compute(this%time,dwdy(1:nx,1:ny,1:nz),this%dwdy)
call pavg_compute(this%time,dwdz(1:nx,1:ny,1:nz),this%dwdz)
!call pavg_compute(this%time,dpdz(1:nx,1:ny,1:nz),this%dpdz)
!call pavg_compute(this%time,dpdz_uv(1:nx,1:ny,1:nz),this%dpdz_uv)

call pavg_compute(this%time,pres_real(1:nx,1:ny,1:nz),this%pres_real)
call pavg_compute(this%time,p(1:nx,1:ny,1:nz),this%p)

call pavg_compute(this%time,fxa(1:nx,1:ny,1:nz),this%fx)
call pavg_compute(this%time,fya(1:nx,1:ny,1:nz),this%fy)
call pavg_compute(this%time,fza(1:nx,1:ny,1:nz),this%fz)

!if (coord==0) then
!call pavg_compute_ws(this%time,-txz(1:nx,1:ny,1),this%twx)
!call pavg_compute_ws(this%time,-tyz(1:nx,1:ny,1),this%twy)
!endif

this%np(bin) = this%np(bin)+1._rprec

!if (coord==0) then
!write(*,*) jt_total,total_time,this%time,bin,this%np(bin),&
!    sum(twxpp)/nx/ny,this%twxpp(bin)/this%np(bin)
!endif

end subroutine compute

!*******************************************************************************
subroutine pavg_compute(time,f,f_pavg)
!*******************************************************************************

use types, only : rprec
use param, only : pi, wave_freq, pavg_nbins
use param, only : jt_total, coord
implicit none

real(rprec), intent(in) :: time
real(rprec), dimension(nx,ny,nz), intent(in) :: f
real(rprec), dimension(nx,ny,nz,pavg_nbins), intent(inout) :: f_pavg

real(rprec), dimension(pavg_nbins) :: np
real(rprec), dimension(nx,nz) :: f_avg
real(rprec) :: time_shift
real(rprec) :: Tf
integer :: bin

!f_avg(:,:) = sum(f(:,:,:),2)/ny

Tf = 2._rprec*pi/wave_freq

bin = ceiling(pavg_nbins/Tf*modulo(time,Tf))
f_pavg(:,:,:,bin) = f_pavg(:,:,:,bin) + f(:,:,:)

end subroutine pavg_compute

!*******************************************************************************
subroutine pavg_compute2d(time,f2,f_pavg2)
!*******************************************************************************

use types, only : rprec
use param, only : pi, wave_freq, pavg_nbins
use param, only : jt_total, coord
implicit none

real(rprec), intent(in) :: time
real(rprec), dimension(nx,ny), intent(in) :: f2
real(rprec), dimension(nx,ny,pavg_nbins), intent(inout) :: f_pavg2

real(rprec), dimension(pavg_nbins) :: np
real(rprec) :: time_shift
real(rprec) :: Tf
integer :: bin

Tf = 2._rprec*pi/wave_freq

bin = ceiling(pavg_nbins/Tf*modulo(time,Tf))
f_pavg2(:,:,bin) = f_pavg2(:,:,bin) + f2(:,:)

end subroutine pavg_compute2d

!!*******************************************************************************
!subroutine pavg_compute_ws(time,f,f_pavg)
!!*******************************************************************************
!
!use types, only : rprec
!use param, only : pi, nu_molec, pulse_freq, pavg_nbins
!
!use param, only : jt_total, coord
!implicit none
!
!real(rprec), intent(in) :: time
!real(rprec), dimension(nx,ny), intent(in) :: f
!real(rprec), dimension(pavg_nbins), intent(inout) :: f_pavg
!
!real(rprec), dimension(pavg_nbins) :: np
!real(rprec) :: f_avg
!real(rprec) :: time_shift
!real(rprec) :: Tf
!integer :: bin
!
!f_avg = sum(f)/nx/ny
!Tf = 2._rprec*pi*nu_molec/pulse_freq
!! shift time so t=0 corresponds with the peak for cosine forcing
!time_shift = time - Tf/4._rprec
!bin = ceiling(pavg_nbins/Tf*modulo(time_shift,Tf))
!f_pavg(bin) = f_pavg(bin) + f_avg
!
!end subroutine pavg_compute_ws



!*******************************************************************************
subroutine finalize(this)
!*******************************************************************************
use grid_m
use param, only : write_endian, lbz, path, coord, nproc, pavg_nbins
use string_util
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array,MPI_SYNC_DOWNUP
use param, only : ierr,comm
#endif
implicit none

class(pavg_t), intent(inout) :: this

character(64) :: bin_ext

character(64) :: fname_vel,fname_vel2,fname_tau,fname_sgs,fname_vgrad,fname_ws,&
        fname_rs, fname_f

real(rprec), allocatable, dimension(:,:,:,:) :: up2, vp2, wp2, upvp, upwp, vpwp, &
        upwp_uv, vpwp_uv,wpwp_uv

integer :: jx, jy, jz
integer :: fid

! Common file name
fname_vel = path // 'output/vel_pavg'
fname_vel2 = path // 'output/vel2_pavg'
fname_tau = path // 'output/tau_pavg'
fname_vgrad = path // 'output/vgrad_pavg'
fname_rs = path // 'output/rs_pavg'
fname_f = path // 'output/f_pavg'
!fname_ws = path // 'output/ws_pavg.bin'

call string_splice(bin_ext, '.c', coord, '.bin')

call string_concat(fname_vel, bin_ext)
call string_concat(fname_vel2, bin_ext)
call string_concat(fname_tau, bin_ext)
call string_concat(fname_vgrad, bin_ext)
call string_concat(fname_rs, bin_ext)
call string_concat(fname_f, bin_ext)

! Final checkpoint all restart data
call this%checkpoint()
do jz = 1,nz
 do jy = 1,ny
  do jx = 1,nx
    this%u(jx,jy,jz,:) = this%u(jx,jy,jz,:) / this%np(:)
    this%v(jx,jy,jz,:) = this%v(jx,jy,jz,:) / this%np(:)
    this%w(jx,jy,jz,:) = this%w(jx,jy,jz,:) / this%np(:)
    this%u_w(jx,jy,jz,:) = this%u_w(jx,jy,jz,:) / this%np(:)
    this%v_w(jx,jy,jz,:) = this%v_w(jx,jy,jz,:) / this%np(:)
    this%w_uv(jx,jy,jz,:) = this%w_uv(jx,jy,jz,:) / this%np(:)
   
    this%uu(jx,jy,jz,:) = this%uu(jx,jy,jz,:) / this%np(:)
    this%vv(jx,jy,jz,:) = this%vv(jx,jy,jz,:) / this%np(:)
    this%ww(jx,jy,jz,:) = this%ww(jx,jy,jz,:) / this%np(:)
    this%uv(jx,jy,jz,:) = this%uv(jx,jy,jz,:) / this%np(:)
    this%uw(jx,jy,jz,:) = this%uw(jx,jy,jz,:) / this%np(:)
    this%vw(jx,jy,jz,:) = this%vw(jx,jy,jz,:) / this%np(:)
  
    this%uw_uv(jx,jy,jz,:) = this%uw_uv(jx,jy,jz,:) / this%np(:)
    this%vw_uv(jx,jy,jz,:) = this%vw_uv(jx,jy,jz,:) / this%np(:)
    this%ww_uv(jx,jy,jz,:) = this%ww_uv(jx,jy,jz,:) / this%np(:)

    this%txx(jx,jy,jz,:) = this%txx(jx,jy,jz,:) / this%np(:)
    this%tyy(jx,jy,jz,:) = this%tyy(jx,jy,jz,:) / this%np(:)
    this%tzz(jx,jy,jz,:) = this%tzz(jx,jy,jz,:) / this%np(:)
    this%txy(jx,jy,jz,:) = this%txy(jx,jy,jz,:) / this%np(:)
    this%txz(jx,jy,jz,:) = this%txz(jx,jy,jz,:) / this%np(:)
    this%tyz(jx,jy,jz,:) = this%tyz(jx,jy,jz,:) / this%np(:)
   
    this%dudx(jx,jy,jz,:) = this%dudx(jx,jy,jz,:) / this%np(:)
    this%dudy(jx,jy,jz,:) = this%dudy(jx,jy,jz,:) / this%np(:)
    this%dudz(jx,jy,jz,:) = this%dudz(jx,jy,jz,:) / this%np(:)
    this%dvdx(jx,jy,jz,:) = this%dvdx(jx,jy,jz,:) / this%np(:)
    this%dvdy(jx,jy,jz,:) = this%dvdy(jx,jy,jz,:) / this%np(:)
    this%dvdz(jx,jy,jz,:) = this%dvdz(jx,jy,jz,:) / this%np(:)
    this%dwdx(jx,jy,jz,:) = this%dwdx(jx,jy,jz,:) / this%np(:)
    this%dwdy(jx,jy,jz,:) = this%dwdy(jx,jy,jz,:) / this%np(:)
    this%dwdz(jx,jy,jz,:) = this%dwdz(jx,jy,jz,:) / this%np(:)
 !   this%dpdz(jx,jy,jz,:) = this%dpdz(jx,jy,jz,:) / this%np(:)
 !   this%dpdz_uv(jx,jy,jz,:) = this%dpdz_uv(jx,jy,jz,:) / this%np(:)

    this%p(jx,jy,jz,:) = this%p(jx,jy,jz,:) / this%np(:)
    this%pres_real(jx,jy,jz,:) = this%pres_real(jx,jy,jz,:) / this%np(:)
 
    this%fx(jx,jy,jz,:) = this%fx(jx,jy,jz,:) / this%np(:)
    this%fy(jx,jy,jz,:) = this%fy(jx,jy,jz,:) / this%np(:)
    this%fz(jx,jy,jz,:) = this%fz(jx,jy,jz,:) / this%np(:)

   enddo
 enddo
enddo

! Write binary data

open(newunit=fid, file=fname_vel, form='unformatted', convert=write_endian,        &
    access='direct', recl=nx*ny*nz*pavg_nbins*rprec)
write(fid,rec=1) this%u(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=2) this%v(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=3) this%w(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=4) this%u_w(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=5) this%v_w(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=6) this%w_uv(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=7) this%p(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=8) this%pres_real(1:nx,1:ny,1:nz,1:pavg_nbins)
!write(fid,rec=9) this%dpdz(1:nx,1:ny,1:nz,1:pavg_nbins)
!write(fid,rec=10) this%dpdz_uv(1:nx,1:ny,1:nz,1:pavg_nbins)
close(fid)

open(newunit=fid, file=fname_vel2, form='unformatted', convert=write_endian,        &
    access='direct', recl=nx*ny*nz*pavg_nbins*rprec)
write(fid,rec=1) this%uu(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=2) this%vv(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=3) this%ww(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=4) this%uv(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=5) this%uw(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=6) this%vw(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=7) this%uw_uv(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=8) this%vw_uv(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=9) this%ww_uv(1:nx,1:ny,1:nz,1:pavg_nbins)
close(fid)

open(newunit=fid, file=fname_tau, form='unformatted', convert=write_endian,        &
    access='direct', recl=nx*ny*nz*pavg_nbins*rprec)
write(fid,rec=1) this%txx(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=2) this%tyy(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=3) this%tzz(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=4) this%txy(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=5) this%txz(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=6) this%tyz(1:nx,1:ny,1:nz,1:pavg_nbins)
close(fid)


open(newunit=fid, file=fname_vgrad, form='unformatted', convert=write_endian,        &
   access='direct', recl=nx*ny*nz*pavg_nbins*rprec)
write(fid,rec=1) this%dudx(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=2) this%dudy(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=3) this%dudz(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=4) this%dvdx(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=5) this%dvdy(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=6) this%dvdz(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=7) this%dwdx(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=8) this%dwdy(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=9) this%dwdz(1:nx,1:ny,1:nz,1:pavg_nbins)
close(fid)

! Do the Reynolds stress calculations afterwards. Now we can interpolate w and
! ww to the uv grid and do the calculations. We have already written the data to
! the files so we can overwrite now
allocate( up2(1:nx,1:ny,1:nz,pavg_nbins) )
allocate( vp2(1:nx,1:ny,1:nz,pavg_nbins) )
allocate( wp2(1:nx,1:ny,1:nz,pavg_nbins) )
allocate( upvp(1:nx,1:ny,1:nz,pavg_nbins) )
allocate( upwp(1:nx,1:ny,1:nz,pavg_nbins) )
allocate( vpwp(1:nx,1:ny,1:nz,pavg_nbins) )

allocate( upwp_uv(1:nx,1:ny,1:nz,pavg_nbins) )
allocate( vpwp_uv(1:nx,1:ny,1:nz,pavg_nbins) )
allocate( wpwp_uv(1:nx,1:ny,1:nz,pavg_nbins) )

up2 = this%uu - this%u * this%u
vp2 = this%vv - this%v * this%v
wp2 = this%ww - this%w * this%w
upvp = this%uv - this%u * this%v

!! using u_w and v_w below instead of u and v ensures that the Reynolds
!! stresses are on the same grid as the squared velocities (i.e., w-grid)
upwp = this%uw - this%u_w * this%w
vpwp = this%vw - this%v_w * this%w

!! using w_uv instead of w to insure that Reynolds stresses are in the same
!! grid as the squared velocities (i.e., uv-grid)
upwp_uv = this%uw_uv - this%u * this%w_uv
vpwp_uv = this%vw_uv - this%v * this%w_uv
wpwp_uv = this%ww_uv - this%w_uv * this%w_uv

open(newunit=fid, file=fname_rs, form='unformatted', convert=write_endian,        &
   access='direct', recl=nx*ny*nz*pavg_nbins*rprec)
write(fid,rec=1) up2(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=2) vp2(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=3) wp2(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=4) upvp(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=5) upwp(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=6) vpwp(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=7) upwp_uv(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=8) vpwp_uv(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=9) wpwp_uv(1:nx,1:ny,1:nz,1:pavg_nbins)
close(fid)


open(newunit=fid, file=fname_f, form='unformatted', convert=write_endian,        &
    access='direct', recl=nx*ny*nz*pavg_nbins*rprec)
write(fid,rec=1) this%fx(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=2) this%fy(1:nx,1:ny,1:nz,1:pavg_nbins)
write(fid,rec=3) this%fz(1:nx,1:ny,1:nz,1:pavg_nbins)
close(fid)
!
!open(newunit=fid, file=fname_ws, form='unformatted', convert=write_endian,        &
!    access='direct', recl=pavg_nbins*rprec)
!write(fid,rec=1) this%twx(1:pavg_nbins)
!write(fid,rec=2) this%twxbar(1:pavg_nbins)
!write(fid,rec=3) this%twxpp(1:pavg_nbins)
!write(fid,rec=4) this%twy(1:pavg_nbins)
!write(fid,rec=5) this%twybar(1:pavg_nbins)
!write(fid,rec=6) this%twypp(1:pavg_nbins)
!close(fid)
!
!endif

end subroutine finalize

!*******************************************************************************
subroutine checkpoint(this)
!*******************************************************************************
!
! This subroutine writes the restart data and is to be called by 'checkpoint'
! for intermediate checkpoints and by 'pavg_finalize' at the end of the
! simulation.
!
use param, only : write_endian, coord
use string_util
implicit none

class(pavg_t), intent(inout) :: this

character(64) :: fname
integer :: fid

fname = checkpoint_pavg_file
call string_concat( fname, '.c', coord)

!  Write data to pavg.out
open(newunit=fid, file=fname, action='write', position='rewind',form='unformatted',      &
    convert=write_endian)
write(fid) this%np
write(fid) this%u
write(fid) this%v
write(fid) this%w
write(fid) this%u_w
write(fid) this%v_w
write(fid) this%w_uv
write(fid) this%uu
write(fid) this%vv
write(fid) this%ww
write(fid) this%uv
write(fid) this%uw
write(fid) this%vw

write(fid) this%uw_uv
write(fid) this%vw_uv
write(fid) this%ww_uv

write(fid) this%txx
write(fid) this%tyy
write(fid) this%tzz
write(fid) this%txy
write(fid) this%txz
write(fid) this%tyz

write(fid) this%dudx
write(fid) this%dudy
write(fid) this%dudz
write(fid) this%dvdx
write(fid) this%dvdy
write(fid) this%dvdz
write(fid) this%dwdx
write(fid) this%dwdy
write(fid) this%dwdz
!write(fid) this%dpdz
!write(fid) this%dpdz_uv

write(fid) this%p
write(fid) this%pres_real

write(fid) this%fx
write(fid) this%fy
write(fid) this%fz
close(fid)

end subroutine checkpoint

end module phase_average

