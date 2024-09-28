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
module wave_spectrum
!*******************************************************************************
use types, only : rprec

save
private rprec
public  

contains

!*******************************************************************************
subroutine spectrum_calc ()
!*******************************************************************************
! This subroutine will initialize all variables for the ocean spectrum, this will
! only be called once at the beginning of the simulation in initialize.f90.

use param, only : L_x, L_y, pi, nx, ny, dx, dy, alpha_spec, wp_spec, kp_spec,   &
                  ld, z_i, u_star, nu_molec, path, theta_main, L_platform
use sim_param, only: K_wave, D_spread, S_kx_ky, phi, omega_wave, eta,           &
                     eta_hat_o, eta_hat_o_filter, S_filter
implicit none

complex(kind=8) :: ii
integer :: jx, jy
real(rprec) :: grav, delta_kx, delta_ky, k_filter, omega_filter, sigma_spec, &
               L_x_wave, L_y_wave, K_platform
real(rprec), dimension(nx,ny) :: Theta_wave,  S_k_D, S_k, S_omega, kx_Wave,  &
                                 ky_wave
ii = cmplx(0.0d0,1.0d0)
grav = 9.81_rprec

! Since all variables are normalized by z_i,u_star and ocean spectrum is built
! for dimensional values, we must have everything with correct dimensions
! This is assuming z_i is equal to lambda_p

L_x_wave = L_x*z_i
L_y_wave = L_y*z_i
wp_spec = wp_spec*u_star/z_i

do jy = -ny/2+1, ny/2
        ky_wave(:,jy+ny/2) = (2_rprec*pi/L_y_wave)*(jy)
end do

do jx = -nx/2, nx/2-1
        kx_wave(jx+nx/2+1,:) = (2_rprec*pi/L_x_wave)*(jx)
end do

do jy = 1,ny
 do jx = 1,nx 
        K_wave(jx,jy) = sqrt(kx_wave(jx,jy)**2 + ky_wave(jx,jy)**2)
        omega_wave(jx,jy) = sqrt(grav*K_wave(jx,jy))
        Theta_wave(jx,jy) =  atan2(ky_wave(jx,jy),kx_wave(jx,jy))
 end do
end do
delta_kx = abs(kx_wave(2,1) - kx_wave(1,1))
delta_ky = abs(ky_wave(1,2) - ky_wave(1,1))

call wave_spectrum_read ()

! Creating the 3D spectrum using the JONSWAP specturm with directional 
! spreading (Same as Yang, et al 2014). The directional spreading has 
! the capability to introduce a main wave direction.
theta_main = theta_main*pi/180_rprec
K_platform = pi/L_platform   ! This has to be with dimensions

do jy = 1, ny
   do jx = 1, nx

      if (omega_wave(jx,jy) <= wp_spec) then
          sigma_spec = 0.07_rprec
      else 
          sigma_spec = 0.09_rprec
      end if

      if (abs(Theta_wave(jx,jy)) <= pi/2.0_rprec) then
          D_spread(jx,jy) = (2.0_rprec/pi)*cos(Theta_wave(jx,jy) - &
                            theta_main)**2
      else
          D_spread(jx,jy) = 0.0
      end if
 
      S_omega(jx,jy) = (alpha_spec*grav**2) / omega_wave(jx,jy)**5 * &
             exp(-1.25_rprec*(wp_spec/omega_wave(jx,jy))**4) * (3.3_rprec)** &
             exp(-((omega_wave(jx,jy) - wp_spec)**2 / (2.0_rprec*sigma_spec** &
             2*wp_spec**2)))
      S_k(jx,jy) = S_omega(jx,jy) * 0.5_rprec * sqrt(grav/K_wave(jx,jy)) 
      S_k_D(jx,jy) = S_k(jx,jy)*D_spread(jx,jy)
      S_kx_ky(jx,jy) = S_k_D(jx,jy)* (1.0_rprec/K_wave(jx,jy))
      S_kx_ky(nx/2+1,ny/2) = 0.0_rprec
      eta_hat_o(jx,jy) = sqrt(2*S_kx_ky(jx,jy)*delta_kx*delta_ky)*exp(ii*phi(jx,jy))

      ! This is for offshore turbine angles. We need to calcualte the filtered ocean 
      ! surface using the size of a typical platform (from input file) such that the 
      ! angles are not to small
      S_filter(jx,jy) = S_kx_ky(jx,jy)

      if (K_wave(jx,jy) >= K_platform) then        
        
      S_filter(jx,jy) = 0.0_rprec
      eta_hat_o_filter(jx,jy) = sqrt(2*S_filter(jx,jy)*delta_kx*delta_ky)*exp(ii*phi(jx,jy))
      end if 
         
   end do
end do

! We need to calculate eta(-k) and then take the complex conjugate. So we need to grab 
! eta_hat_o rotate it 180 degrees and apply conj. Therefore, eta(-k) is the complex 
! conjugate mirror of eta(k)

!do jx = 1,nx/2-1
!   do jy = 1,ny
!        eta_hat_o_minusK(jx, jy) = conjg(eta_hat_o(nx-jx+1,ny-jy+1))
!   end do
!end do


end subroutine spectrum_calc

!*******************************************************************************
subroutine wave_spectrum_checkpoint ()
!*******************************************************************************
! This subroutine will store the random phases and previous wave surface (eta) 
! for restarting simulation

use param, only : write_endian, path, checkpoint_wave_spectrum
use sim_param, only :  phi, eta
implicit none
character(64) :: fname_wave

fname_wave = checkpoint_wave_spectrum
open(11, file=fname_wave, form='unformatted', convert=write_endian, &
     status='unknown', position='rewind')
write (11) phi(:,:) , eta(:,:)
close(11)

end subroutine wave_spectrum_checkpoint

!*******************************************************************************
subroutine wave_spectrum_read ()
!*******************************************************************************
! This subroutine checks if there is a restart file to read the random phase and
! the surface distribution. If there is not, it creates it. The surface 
! distribution is need it available for the calculation of deta_dt in 
!wallstress.f90
use param, only : write_endian, path, checkpoint_wave_spectrum, read_endian,   &
                  nx, ny, pi
use sim_param, only :  phi, eta
implicit none

character(64) :: fname_wave
logical :: wave_spectrum_file_flag
real(rprec), dimension(nx,ny) :: rand
allocate(checkpoint_wave_spectrum , source = path // 'wave_spectrum.out')

! Creating the random numbers
call init_random_seed
call random_number(rand)

inquire (file='ocean_spectrum.out', exist=wave_spectrum_file_flag)
if  (.not. wave_spectrum_file_flag)  then
        write(*,*) 'Creating random phases for ocean spectrum'
        phi(:nx,:ny) = rand(:nx,:ny)*2*pi
else
        write(*,*) 'Reading the random phases for ocean spectrum and previous wave surface'
        open(12, file=checkpoint_wave_spectrum, form='unformatted', convert=read_endian)
        read(12) phi(:,:), eta(:,:)
        close(12)
end if


end subroutine wave_spectrum_read


!*******************************************************************************
subroutine frequency_shift(f,f_shift)
!*******************************************************************************
!
! This subroutine is used to shift and array of size nx,ny, such that the zero
! frequency is located at the corner instead of the center of the domain.
! This subroutine mimics fftshift from matlab.
!
use types, only : rprec
use param, only : nx, ny
implicit none

complex(rprec), dimension(:,:), intent(in) :: f
complex(rprec), dimension(:,:), intent(inout) :: f_shift

f_shift(nx/2+1:, ny/2+1:) = f(:nx/2, :ny/2)
f_shift(nx/2+1:, :ny/2) =  f(:nx/2, ny/2+1:)
f_shift(:nx/2, ny/2+1:) =  f(nx/2+1:, :ny/2)
f_shift(:nx/2, :ny/2) =  f(nx/2+1:, ny/2+1:)


end subroutine frequency_shift

end module wave_spectrum
