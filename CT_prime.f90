!!  Copyright (C) 2010-2026  Johns Hopkins University
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
!*******************************************************************************

module CT_prime
!*******************************************************************************
! This module computes:
! - Induction factor, a
! - Tip Speed Ratio (TSR)
! - Thrust coefficient, C_T
! - Local thrust coefficient C'_T 
! from Ud and pitch angle theta
! Reference: "Implementation of Turbine Control for LES" (Yishu Zhang)
!*******************************************************************************

use types, only : rprec
use param
use grid_m
use messages
use string_util
! use turbines, only: nturbs, disk_avg_vel, theta2
use functions, only : count_lines
use stat_defs, only : wind_farm
#ifdef PPMPI
use mpi_defs, only : MPI_SYNC_DOWNUP, mpi_sync_real_array
#endif

implicit none
private

! Public subroutines and functions
public :: compute_induction, compute_TSR, compute_CT_from_TSR, compute_CT_local_from_CT
public :: U_inf_from_Ud, Udtr_from_Utr
public :: ctprime_init, Filt_step

character(*), parameter :: mod_name = 'ct_prime'

!-----------------------------------------------------------
! Input parameters
real(rprec), public :: CT_op          ! Thrust coefficient at transition II to III
real(rprec), public :: CT_co          ! Thrust coefficient at cut-out
real(rprec), public :: Utr            ! Transitional free-stream velocity II to III
real(rprec), public :: U_co           ! Cut-out free-stream velocity
real(rprec), public :: TSR_op         ! TSR at transition II to III
real(rprec), public :: TSR_co         ! TSR at cut-out


!-----------------------------------------------------------
! Internal fit parameters
real(rprec), public :: alpha_Ct = 0.7_rprec         ! This fit parameter is an Input
real(rprec), private :: beta_Ct  = 1.0_rprec        ! This fit parameter is not user-defined

!-----------------------------------------------------------
! TSR fitting parameters (computed in ctprime_init)
real(rprec), public :: P, Q
real(rprec), public :: eps_safe

contains

!*******************************************************************************
subroutine ctprime_init()
!*******************************************************************************
! Initialize the CT prime module and compute TSR fitting parameters
  implicit none
  P = Utr / (U_co - Utr) * (TSR_op - TSR_co)
  Q = U_co / Utr * P - TSR_op
  eps_safe = 1.0e-12_rprec
end subroutine ctprime_init

!*******************************************************************************
subroutine Filt_step(f_old, f_new, dt, T, f_filtered)
!*******************************************************************************
! Apply one step of exponential moving average (EMA) filtering.
!
! Inputs:
!   f_old      - previous filtered value (f_tilde[n-1])
!   f_new      - new unfiltered value at current time step (f[n])
!   dt         - physical time step (s)
!   T          - filter time constant (s)
!
! Output:
!   f_filtered - updated filtered value (f_tilde[n])
!
! Formula:
!   eps = (dt/T) / (1 + dt/T)
!   f_tilde[n] = (1 - eps)*f_tilde[n-1] + eps*f[n]
!
! Notes:
!   - If T <= 0, no filtering is applied (f_tilde = f_new).
!   - This subroutine is numerically stable for large dt/T.
!*******************************************************************************

  implicit none
  real(rprec), intent(in)  :: f_old     ! previous filtered value
  real(rprec), intent(in)  :: f_new     ! current unfiltered value
  real(rprec), intent(in)  :: dt        ! time step
  real(rprec), intent(in)  :: T         ! filter time constant
  real(rprec), intent(out) :: f_filtered ! updated filtered value

  real(rprec) :: eps

  if (T > 0.0_rprec) then
     eps = (dt / T) / (1.0_rprec + dt / T)
  else
     eps = 1.0_rprec
  end if

  f_filtered = (1.0_rprec - eps) * f_old + eps * f_new

end subroutine Filt_step


!*******************************************************************************
function compute_induction(Ud, theta, CT_op, U_trans, alpha, beta) result(a)
!*******************************************************************************
! Compute induction factor a based on disk-averaged velocity and pitch
  implicit none
  real(rprec), intent(in) :: Ud, theta, CT_op, U_trans, alpha
  real(rprec), intent(in), optional :: beta
  real(rprec) :: a
  real(rprec) :: beta_val, Ud_tr, num, term1, term2, denom
  real(rprec) :: theta_rad

  ! Default beta
  if (present(beta)) then
     beta_val = beta
  else
     beta_val = 1.0_rprec
  end if

  ! thera from degree to radius
  theta_rad = theta * pi / 180.0_rprec

  ! Transitional disk-averaged velocity
  Ud_tr = 0.5_rprec * U_trans * (1.0_rprec + sqrt(max(1.0_rprec - CT_op*cos(theta_rad), eps_safe)))
  ! print*, 'Transicent Ud', Ud_tr
  if (Ud <= Ud_tr) then
     ! Region II
     a = (1.0_rprec - sqrt(1.0_rprec - CT_op*cos(theta_rad))) / 2.0_rprec
  !   print*, 'Region II', 'Ud =',Ud, 'pitching angle =', theta_rad, 'a=',a
  else
     ! Region III
     num = 4.0_rprec*alpha*Ud - CT_op*cos(theta_rad)
     term1 = alpha*U_trans + alpha*Ud - 1.0_rprec
     term2 = term1**2 - (4.0_rprec*alpha*Ud - CT_op*cos(theta_rad))*(alpha*U_trans - 1.0_rprec)
     denom = 2.0_rprec*term1 + 2.0_rprec*sqrt(max(term2, eps_safe))
     a = 1.0_rprec - num / denom
  !   print*, 'Region III', 'Ud =',Ud, 'pitching angle =', theta_rad, 'a=',a
  end if

end function compute_induction

!*******************************************************************************
function U_inf_from_Ud(Ud, a) result(Uinf)
!*******************************************************************************
! Compute free-stream velocity from disk-averaged velocity and induction
  implicit none
  real(rprec), intent(in) :: Ud, a
  real(rprec) :: Uinf

  if (1.0_rprec - a > 1.0e-12_rprec) then
     Uinf = Ud / (1.0_rprec - a)
  else
     Uinf = Ud
  end if

end function U_inf_from_Ud

!*******************************************************************************
function Udtr_from_Utr(Utr, theta, CT_op) result(Udtr)
!*******************************************************************************
! Compute transitional disk-averaged velocity
  implicit none
  real(rprec), intent(in) :: Utr, theta, CT_op
  real(rprec) :: Udtr
  real(rprec) :: theta_rad
    ! thera from degree to radius
  theta_rad = theta * pi / 180.0_rprec

  Udtr = 0.5_rprec * Utr * (1.0_rprec + sqrt(max(1.0_rprec - CT_op*cos(theta_rad), eps_safe)))

end function Udtr_from_Utr

!*******************************************************************************
function compute_TSR(Uinf, TSR_op, U_trans, U_co, P, Q) result(TSR_val)
!*******************************************************************************
! Compute TSR based on free-stream velocity for a single turbine
  implicit none
  real(rprec), intent(in) :: Uinf
  real(rprec), intent(in) :: TSR_op, U_trans, U_co, P, Q
  real(rprec) :: TSR_val

  if (Uinf <= U_trans) then
     ! Region II
     TSR_val = TSR_op
  else
     ! Region III
     TSR_val = max((U_co / Uinf) * P - Q, 0.0_rprec)
  end if

end function compute_TSR

!*******************************************************************************
function compute_CT_from_TSR(TSR, CT_op, TSR_op, U_co, U_trans, P, Q, alpha, beta) result(CT_val)
!*******************************************************************************
! Compute thrust coefficient based on scalar TSR
  implicit none
  real(rprec), intent(in) :: TSR
  real(rprec), intent(in) :: CT_op, TSR_op, U_co, U_trans, P, Q, alpha
  real(rprec), intent(in), optional :: beta
  real(rprec) :: CT_val
  real(rprec) :: beta_val


  ! Default beta = 1
  if (present(beta)) then
     beta_val = beta
  else
     beta_val = 1.0_rprec
  end if

  ! Only support beta = 1 for now
  if (beta_val /= 1.0_rprec) then
     print *, "Warning: Only beta=1 is supported; ignoring other values."
  end if

  ! Region II / Region III logic
  if (TSR == TSR_op) then
     ! Region II: constant CT
     CT_val = CT_op
  else if (TSR < TSR_op) then
     ! Region III: above-rated region
     CT_val = (CT_op * (TSR + Q)) / ( (1.0_rprec - alpha*U_trans)*(TSR + Q) + alpha*U_co*P )
  else 
     print *, "Warning: TSR cannot be greater than TSR_op! Region 1.5 is not considered"     
  end if

end function compute_CT_from_TSR

!*******************************************************************************
function compute_CT_local_from_CT(theta, CT) result(CT_local)
!*******************************************************************************
! Compute local thrust coefficient C'_T from scalar CT and azimuth angle theta
  implicit none
  real(rprec), intent(in) :: theta
  real(rprec), intent(in) :: CT
  real(rprec) :: CT_local
  real(rprec) :: sqrt_term
  real(rprec) :: theta_rad

  ! thera from degree to radius
  theta_rad = theta * pi / 180.0_rprec

  ! Guard sqrt argument against negative due to rounding
  sqrt_term = sqrt(max(1.0_rprec - CT*cos(theta_rad), eps_safe))

  ! Local thrust coefficient formula
  CT_local = (4.0_rprec - 4.0_rprec*sqrt_term) / (1.0_rprec + sqrt_term)

  ! Ensure non-negative
  if (CT_local < 0.0_rprec) CT_local = 0.0_rprec

end function compute_CT_local_from_CT

end module CT_prime
