!!
!!  Copyright (C) 2011-2017  Johns Hopkins University
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
module eqmfit_wm
!*******************************************************************************
! This module has all the calculations of the equilibrium wall model using the
! fit from Meneveau JoT (2020). 

use types, only : rprec
use param, only : nx, ny, ld, dx, dz,  u_star, z_i, pi, zgrid_match, nu_molec,  &
                  vonk, zo, smooth_eqm, dy
use sim_param, only : eqmxz, eqmyz, ur_mag_eqm, u_LES, u, v
use grid_m
use test_filtermodule
use functions, only : cell_indx

implicit none

private  
public :: eqmfit_calc, eqmfit_finalize

contains
!*******************************************************************************
subroutine eqmfit_finalize()
!*******************************************************************************
use param, only : nx, ny
use sim_param, only : txz, tyz
implicit none

call eqmfit_calc ()

txz(1:nx,1:ny,1) = eqmxz(:,:)
tyz(1:nx,1:ny,1) = eqmyz(:,:) 


end subroutine eqmfit_finalize

!*******************************************************************************
subroutine eqmfit_calc ()
!*******************************************************************************
use param, only: lbc_mom
use sim_param, only: dudz, dvdz, u_orb

implicit none
integer :: dz_match_ind
real(rprec) :: dz_match
real(rprec), dimension(ld, ny) :: u_eqm, v_eqm, u1, v1
real(rprec), dimension(nx, ny) :: ur_eqm, vr_eqm, u_avg,Re_delta, beta_1,   &
                                  beta_2, Re_fit, Re_delta2, beta_12,       &
                                  beta_22, Re_fit2, A, A2

! Calculating the velocity input, this is set up to manage any matching location 
! specified in lesgo.conf. (e.g. zgrid_match=1.5 is at 2nd grid point and 
! 2.5 is at 3rd)
dz_match = zgrid_match*dz
dz_match_ind = ceiling(zgrid_match)
u_eqm = u(1:ld,1:ny,dz_match_ind)
v_eqm = v(1:ld,1:ny,dz_match_ind)
call test_filter(u_eqm)
call test_filter(v_eqm)

if (lbc_mom == 4) then
ur_eqm = u_eqm(1:nx,1:ny) - u_orb(1:nx,1:ny)
else
ur_eqm = u_eqm(1:nx,1:ny)
end if

vr_eqm = v_eqm(1:nx,1:ny)
ur_mag_eqm = sqrt((ur_eqm(:,:))**2 + (vr_eqm(:,:))**2)

! Calculating the velocity input for dudz, dvdz at 1st grid point. Velocity is
! filtered and velocity magnitude is calculated.
u1 = u(1:ld,1:ny,1)
v1 = v(1:ld,1:ny,1)
call test_filter(u1)
call test_filter(v1)
u_avg = sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)

! Calculating the parameters for the EQM model. We calcualte 2 sets because 
! the second is used for the calculation of dudz,dvdz which is at 1st grid point
u_LES  = sqrt((ur_eqm(:,:))**2 + (vr_eqm(:,:))**2)
Re_delta = ((u_LES(:,:)*u_star)*((dz_match)*z_i))/nu_molec
beta_1 = (1_rprec + 0.155_rprec*Re_delta(:,:)**(-0.03_rprec))**(-1_rprec)
beta_2 = 1.7_rprec - (1_rprec + 36_rprec*Re_delta(:,:)**(-0.75_rprec))**(-1_rprec)
Re_fit = (0.005_rprec)**(beta_1(:,:)-0.5_rprec)*Re_delta(:,:)**(beta_1(:,:))* &
                    (1_rprec + (0.005_rprec*Re_delta(:,:))**(-beta_2(:,:)))** &
                    ((beta_1(:,:) -0.5_rprec)/beta_2(:,:))

Re_delta2 = ((u_avg(:,:)*u_star)*(0.5_rprec*dz*z_i))/nu_molec
beta_12 = (1_rprec + 0.155_rprec*Re_delta2(:,:)**(-0.03_rprec))**(-1_rprec)
beta_22 = 1.7_rprec - (1_rprec + 36_rprec*Re_delta2(:,:)**(-0.75_rprec))**(-1_rprec)
Re_fit2 = (0.005_rprec)**(beta_12(:,:)-0.5_rprec)*Re_delta2(:,:)**(beta_12(:,:))* &
                    (1_rprec + (0.005_rprec*Re_delta2(:,:))**(-beta_22(:,:)))** &
                    ((beta_12(:,:) -0.5_rprec)/beta_22(:,:))

if (smooth_eqm) then
        eqmxz = -(Re_fit(:,:)/Re_delta(:,:))**2 * ur_eqm(:,:)*u_LES(:,:)
        eqmyz = -(Re_fit(:,:)/Re_delta(:,:))**2 * vr_eqm(:,:)*u_LES(:,:)
        dudz(1:nx,1:ny,1) = 1/(0.5_rprec*dz*vonk)*(Re_fit2(:,:)/ &
                            Re_delta2(:,:))*u(1:nx,1:ny,1)
        dvdz(1:nx,1:ny,1) = 1/(0.5_rprec*dz*vonk)*(Re_fit2(:,:)/ &
                            Re_delta2(:,:))*v(1:nx,1:ny,1)
else
        A = ((Re_fit(:,:)/Re_delta(:,:))**(6.0_rprec) + (2.5_rprec* &
            log((dz_match)/zo))**(-6.0_rprec))**(0.33333_rprec)
        eqmxz = - A(:,:) * ur_eqm(:,:)*u_LES(:,:)
        eqmyz = - A(:,:) * vr_eqm(:,:)*u_LES(:,:)

        A2 = ((Re_fit2(:,:)/Re_delta2(:,:))**(6.0_rprec) + (2.5_rprec* &
             log((0.5_rprec*dz)/zo))**(-6.0_rprec))**(0.16666_rprec)
        dudz(1:nx,1:ny,1) = A2(:,:)/(vonk*0.5_rprec*dz)*u(1:nx,1:ny,1)
        dvdz(1:nx,1:ny,1) = A2(:,:)/(vonk*0.5_rprec*dz)*v(1:nx,1:ny,1)
end if

dudz(1:nx,1:ny,1) = merge(0._rprec,dudz(1:nx,1:ny,1),u(1:nx,1:ny,1).eq.0._rprec)
dvdz(1:nx,1:ny,1) = merge(0._rprec,dvdz(1:nx,1:ny,1),v(1:nx,1:ny,1).eq.0._rprec)

end subroutine eqmfit_calc

end module eqmfit_wm
