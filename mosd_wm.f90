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
!! Manuel trial
!*******************************************************************************
module mosd_wm
!*******************************************************************************
! This module has all the calculations of the MOSD model, the equilibirum 
! wall model (EQM) is called from eqmfit_wm.f90 and the windward potential
! model (WPM) is calculated here.

use types, only : rprec
use param, only : wave_type, nx, ny, ld, dx, dz, dt, total_time_dim, total_time,&
                  amp, wave_n, wave_freq, u_star, z_i, pi, kp_spec, zgrid_match,&
                  nu_molec, vonk, zo, smooth_eqm, dy, uns_mosd, coord
use sim_param, only : eqmxz, eqmyz, s_wpmxz, s_wpmyz,uns_wpmxz, uns_wpmyz, eta, &
                      detadx, detady, detadt, u_orb, w_orb, deta2dx2, deta2dy2, &
                      detadxdy, detadydx, eta_filter, grad_eta_mag, dgrad_etadt,&
                      Cx_wave,Cy_wave, ur_mag_wpm, CdotGradGradeta, ur_mag_eqm, &
                      u_LES, u, v, detadx_dt, detady_dt
use sim_param, only : corr_1,corr_2, corr_3, corr_4,corr_5, corr_6
use sim_param, only : eta_hat_o, eta_hat_o_filter, omega_wave
use grid_m
use test_filtermodule
use functions, only : cell_indx
use fft
use derivatives, only : ddx_fd, ddy_fd
use wave_spectrum, only : frequency_shift
use messages, only : error
use eqmfit_wm, only : eqmfit_calc
implicit none

private  
public :: wpm_calc,  mosd_finalize, mono_wave, spectrum_wave

contains
!*******************************************************************************
subroutine mosd_finalize()
!*******************************************************************************
use param, only : nx, ny
use sim_param, only : txz, tyz
implicit none

call wpm_calc ()
call eqmfit_calc ()

if (uns_mosd) then
txz(1:nx,1:ny,1) = eqmxz(:,:) +  s_wpmxz(:,:) + uns_wpmxz(:,:)
tyz(1:nx,1:ny,1) = eqmyz(:,:) +  s_wpmyz(:,:) + uns_wpmyz(:,:)
else
txz(1:nx,1:ny,1) = eqmxz(:,:) +  s_wpmxz(:,:) 
tyz(1:nx,1:ny,1) = eqmyz(:,:) +  s_wpmyz(:,:) 
end if

end subroutine mosd_finalize

!*******************************************************************************
subroutine wave_selec() 
!*******************************************************************************
! This subroutine checks what type of wave will be modeled and builds the 
! surface distribution, gradients, etc.

implicit none
character(*), parameter :: sub_name = 'wave type'

select case(wave_type)

        ! Monochromatic wave
        case (0)
        call mono_wave()
        
        ! Spectrum wave
        case (1)
        call spectrum_wave()

        ! Otherwise, invalid
        case default
        call error (sub_name, 'invalid')

end select

end subroutine wave_selec

!*******************************************************************************
subroutine mono_wave() 
!*******************************************************************************
implicit none
integer :: jx, jy
real(rprec), dimension(ld, ny) :: x_grid
real(rprec), dimension(nx, ny) :: detadx_new

do jy = 1,ny
   do jx= 1,ld
        
      x_grid(jx,jy) = (jx-1)*dx
        
   end do
end do

eta = amp*cos(wave_n*x_grid(1:nx,1:ny) - wave_freq*total_time) 
detadx_new = amp*wave_n*sin(wave_freq*total_time - wave_n*x_grid(1:nx,1:ny))
detady = 0.0_rprec
detadt = amp*wave_freq*sin(wave_n*x_grid(1:nx,1:ny)   - wave_freq*total_time)
u_orb = amp*wave_freq*cos(wave_n*x_grid(1:nx,1:ny)  - wave_freq*total_time)
w_orb = amp*wave_freq*sin(wave_n*x_grid(1:nx,1:ny)  - wave_freq*total_time)
deta2dx2 = -amp*wave_n**2*cos(-wave_freq*total_time + wave_n*x_grid(1:nx,1:ny))
deta2dy2 = 0.0_rprec
detadxdy = 0.0_rprec
detadydx = 0.0_rprec

detadx_dt = (detadx_new(:,:) - detadx(:,:))/dt
detadx = detadx_new(:,:)

end subroutine mono_wave

!*******************************************************************************
subroutine spectrum_wave() 
!*******************************************************************************
implicit none
integer :: jx, jy
real(rprec), dimension(nx, ny) :: u_orb_1, eta_new, eta_1, eta_filter_1,        &
                                  detadx_new, detady_new
complex(rprec), dimension(nx,ny) :: eta_hat, eta_shift, u_orb_hat, u_orb_shift, &
                                    eta_hat_filter, eta_filter_shift
complex(rprec), dimension(nx/2+1,ny) :: eta_hat_2, u_orb_hat_2, eta_hat_filter_2
complex(kind=8) :: ii

ii = cmplx(0.0d0,1.0d0)

do jy = 1,ny
   do jx = 1,nx
    
      eta_hat(jx,jy) = eta_hat_o(jx,jy)*exp(-ii*omega_wave(jx,jy)*total_time_dim)
      eta_hat_filter(jx,jy) = eta_hat_o_filter(jx,jy)*exp(-ii*omega_wave(jx,jy)*&
                              total_time_dim)
      u_orb_hat(jx,jy) = eta_hat(jx,jy)*omega_wave(jx,jy)

   end do
end do

! The spectrum needs to be shifted before going through fft and then only grabing
! half of the spectrum for this type of fft
call frequency_shift(eta_hat,eta_shift)
call frequency_shift(eta_hat_filter,eta_filter_shift)
call frequency_shift(u_orb_hat,u_orb_shift)

eta_hat_2 = eta_shift(:nx/2+1,1:ny)
eta_hat_filter_2 = eta_filter_shift(:nx/2+1,1:ny)
u_orb_hat_2 = u_orb_shift(:nx/2+1,1:ny)

call dfftw_execute_dft_c2r(back_wave, eta_hat_2, eta_1)
call dfftw_execute_dft_c2r(back_wave, eta_hat_filter_2, eta_filter_1)
call dfftw_execute_dft_c2r(back_wave, u_orb_hat_2, u_orb_1)

! Since the surface spectrum and orbital velocities are calculated using dimensions,
! a normalization must be applied. Because of the FFT, the surface spectrum has to be 
! devided by 2
eta_new = eta_1(:,:)/2.0_rprec/z_i
eta_filter = eta_filter_1(:,:)/2.0_rprec/z_i
u_orb = u_orb_1(:,:)/2.0_rprec/u_star
w_orb = 0.0_rprec !This needs to be changed at soem point

! Calculating the time derivative and spatial derivatives of the surface
detadt = (eta_new(:,:) - eta(:,:))/dt
eta = eta_new(:,:)

call ddx_fd(eta, detadx_new)
call ddy_fd(eta, detady_new)
detadx_new(nx,:) = -eta(nx,:)/dx
detady_new(:,ny) = -eta(:,ny)/dy

call ddx_fd(detadx,deta2dx2)
call ddy_fd(detady,deta2dy2)
deta2dx2(nx,:) = -detadx(nx,:)/dx
deta2dy2(:,ny) = -detady(:,ny)/dy

call ddy_fd(detadx,detadxdy)
call ddx_fd(detady,detadydx)
detadxdy(:,ny) = -detadx(:,ny)/dy
detadydx(nx,:) = -detady(nx,:)/dx

detadx_dt = (detadx_new(:,:) - detadx(:,:))/dt
detady_dt = (detady_new(:,:) - detady(:,:))/dt
detadx = detadx_new(:,:)
detady = detady_new(:,:)


end subroutine spectrum_wave

!*******************************************************************************
subroutine wpm_calc ()
!*******************************************************************************
implicit none
integer :: jx, jy
real(rprec), pointer, dimension(:) :: z
real(rprec) :: threshold_wave_speed, k_min, L_x_wave, L_y_wave,   &
               z_dx_low, z_dx_up, z_diff
real(rprec), dimension(ld, ny) :: u_dx, v_dx
real(rprec), dimension(nx, ny) :: grad_eta_mag_new, ur_wpm, vr_wpm, nx_wpm,   &
                                  ny_wpm, H_arg, H_wpm, s_wpm, uns_wpm
nullify(z)
z => grid % z

! Calculating the velocity input for WPM at a height dx. Not using the linear_interp 
! because it does not work. The velocity is also filtered at delta scale.
! NOTE: ensure that height dx lies in the coord = 0.
z_dx_low = cell_indx('k',dz,dx)
z_dx_up = z_dx_low + 1
z_diff = dx - z(z_dx_low)
u_dx = u(1:ld,1:ny,z_dx_low) + (u(1:ld,1:ny,z_dx_up) - u(1:ld,1:ny,z_dx_low))*z_diff/dz
v_dx = v(1:ld,1:ny,z_dx_low) + (v(1:ld,1:ny,z_dx_up) - v(1:ld,1:ny,z_dx_low))*z_diff/dz
call test_filter(u_dx)
call test_filter(v_dx)

call wave_selec()

! Calculating the magnitude of the gradient, the temporal derivative of the gradient
! and the velocity components of the wave. The phase velocity is clipped to the
! velocity of the wave of size 0.25*kp
grad_eta_mag_new = sqrt((detadx(:,:))**2 + (detady(:,:))**2)
dgrad_etadt = (grad_eta_mag_new(:,:) - grad_eta_mag(:,:))/dt
grad_eta_mag = grad_eta_mag_new(:,:)

k_min = 0.25_rprec*(kp_spec/z_i)
threshold_wave_speed = sqrt(9.81_rprec/k_min)/u_star
Cx_wave = -detadt(:,:)*detadx(:,:)*(1/grad_eta_mag(:,:)**2)
Cx_wave = min(max(Cx_wave(:,:),-threshold_wave_speed),threshold_wave_speed)
Cy_wave = -detadt(:,:)*detady(:,:)*(1/grad_eta_mag(:,:)**2)
Cy_wave = min(max(Cy_wave(:,:),-threshold_wave_speed),threshold_wave_speed)

! Calcualting the realtive velocities. u_dx and v_dx are already filtered 
! and we CANNOT filter Cx_wave, Cy_wave because they won't cancel out with 
! the surface gradients specially for spectrum wave
ur_wpm = u_dx(1:nx,1:ny) - Cx_wave(:,:)
vr_wpm = v_dx(1:nx,1:ny) - Cy_wave(:,:)

! Calculating the remainder of the steady WPM component
nx_wpm = detadx(:,:)/grad_eta_mag(:,:)
ny_wpm = detady(:,:)/grad_eta_mag(:,:)
ur_mag_wpm = sqrt((ur_wpm(:,:)*nx_wpm(:,:))**2 + (vr_wpm(:,:)*ny_wpm(:,:))**2)
H_arg = (ur_wpm(:,:)*detadx(:,:) + vr_wpm(:,:)*detady(:,:))
H_wpm = (H_arg(:,:) + abs(H_arg(:,:)))/(2*H_arg(:,:))
s_wpm = (1/pi)*(ur_mag_wpm(:,:)**2)*(grad_eta_mag(:,:)**2)

! Calculating the unsteady WPM component
CdotGradGradeta = Cx_wave(:,:)*(detadx(:,:)*deta2dx2(:,:) + detady(:,:)* &
                  detadxdy(:,:))/grad_eta_mag(:,:) + Cy_wave(:,:)*(detadx(:,:)* &
                  detadydx(:,:) + detady(:,:)*deta2dy2(:,:))/grad_eta_mag(:,:)

uns_wpm = 1/(4*pi)*ur_mag_wpm(1:nx,1:ny)*sqrt(dx**2 + dy**2)*grad_eta_mag(:,:)* &
          (dgrad_etadt(:,:) + CdotGradGradeta(:,:))

! Calculating the components of both the steady and unsteady WPM
uns_wpmxz =  uns_wpm(:,:)*nx_wpm(:,:)*H_wpm(:,:)
uns_wpmyz =  uns_wpm(:,:)*ny_wpm(:,:)*H_wpm(:,:)

s_wpmxz  = -s_wpm(:,:)*nx_wpm(:,:)*H_wpm(:,:)
s_wpmyz  = -s_wpm(:,:)*ny_wpm(:,:)*H_wpm(:,:)

!*******************This section is for debuging of variables********************!

corr_1 = ur_mag_wpm(:,:)**2*ny_wpm(:,:)
corr_2 = grad_eta_mag(:,:)**2*ny_wpm(:,:)
corr_3 = H_wpm(:,:)*ny_wpm(:,:)
corr_4 = ur_mag_wpm(:,:)**2*H_wpm(:,:)*ny_wpm(:,:)
corr_5 = grad_eta_mag(:,:)**2*H_wpm(:,:)*ny_wpm(:,:)
corr_6 = (ur_mag_wpm(:,:)**2)*(grad_eta_mag(:,:)**2)*ny_wpm(:,:)

end subroutine wpm_calc

end module mosd_wm
