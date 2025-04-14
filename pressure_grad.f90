!!
!!  Copyright (C) 2019  Johns Hopkins University
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
module pressure_grad 
!*******************************************************************************
! This module contains all of the subroutines associated with scalar transport
use types, only : rprec
use pid_m
implicit none

save
private

public :: pressure_grad_init, pressure_grad_finalize, pressure_grad_calc

! u_set -> planar-averaged velocity averaged at top 3 layers
! height_set -> height of u2 set point (dimensional)
! Kp, Ki, Kd -> PID controller gains (dimensional)
! integer, public :: pid_time_P = 0
real(rprec), public :: u_set != 15.0
real(rprec), public :: Kp_P != 1e-4
real(rprec), public :: Ki_P != 3.8e-8
real(rprec), public :: Kd_P != 0.0
logical :: nonlinear_error = .true.

! PID controller
type(pid_t) :: pid

! Planar-averaged velocity averaged at top 3 layers
real(rprec) :: u_actual
! nonlinear error (u_set-u_actural)*|u_set-u_actural| 
real(rprec) :: u_control
! controlled pressure gradient which is added to RHS
real(rprec), public :: dp_control
! Length scale for conversion
real(rprec) :: L

contains

!*******************************************************************************
subroutine pressure_grad_init
!*******************************************************************************
! This subroutine initializes the variables for the pressure gradient 
use param, only : z_i, u_star, L_z, nz, coord, nproc, dz, read_endian
use grid_m
logical :: exst
real(rprec) :: e_int
integer :: num_t, fid, i

! Non-dimensionalize
L = L_z/z_i
Ki_P = Ki_P*z_i/u_star
Kd_P = Kd_P*u_star/z_i

! Create nonlinear PID controller
pid = pid_t(Kp_P, Ki_P, Kd_P, u_set, nonlinear_error)
! Create linear PID controller
pid = pid_t(Kp_P, Ki_P, Kd_P, u_set)

inquire (file='pressure_pid.out', exist=exst)
if (exst) then
    open(12, file='pressure_pid.out', form='unformatted', convert=read_endian)
    read(12) e_int
    close(12)
    pid%e_int = e_int
end if
end subroutine pressure_grad_init

!*******************************************************************************
subroutine pressure_grad_finalize
!*******************************************************************************
use param, only : read_endian

open(12, file='pressure_grad_pid.out', form='unformatted', convert=read_endian)
write(12) pid%e_int
close(12)


end subroutine pressure_grad_finalize

!*******************************************************************************
subroutine pressure_grad_calc
!*******************************************************************************
use param, only : MPI_RPREC, comm, ierr, dt, total_time_dim, u_star, jt_total, coord, nproc
use sim_param, only : u, v, RHSx, RHSy, nx, ny, nz
use functions, only : linear_interp
use mpi
use messages


! Ensure there are at least three layers
if (nproc*nz < 3) then
    call error ('pressure_grad/pressure_grad_calc', 'Not enough grid points in the z-direction.')
end if

! Compute the average x-velocity at the top three layers
if (coord == nproc-1) then
    u_actual = sum(u(:,:,nz-2:nz)) / (nx * ny * 3)
!end if
   ! Use PID to get new velocity square

   u_control = pid%advance(u_actual, dt)     
   ! Match the dimension
   dp_control = u_control / L
end if

! Broadcast dp_control from the last process to all others
call MPI_Bcast(dp_control, 1, MPI_DOUBLE_PRECISION, nproc-1, MPI_COMM_WORLD, ierr)

! write(*,*) "dp_control:", dp_control, coord
! Pressure gradient: add forcing to RHS
RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + dp_control  
  
end subroutine pressure_grad_calc

end module pressure_grad
