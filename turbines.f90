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
module turbines
!*******************************************************************************
! This module contains all of the subroutines associated with drag-disk turbines

use types, only : rprec
use param
use grid_m
use messages
use string_util
use turbine_indicator
use functions, only : count_lines
use stat_defs, only : wind_farm
#ifdef PPMPI
use mpi_defs, only : MPI_SYNC_DOWNUP, mpi_sync_real_array
#endif

implicit none

save
private

public :: turbines_init, turbines_forcing, turbine_vel_init, turbines_finalize

character (*), parameter :: mod_name = 'turbines'

! The following values are read from the input file
! number of turbines in the x-direction
integer, public :: num_x
! number of turbines in the y-direction
integer, public :: num_y
! baseline diameter in meters
real(rprec), public :: dia_all
! baseline height in meters
real(rprec), public :: height_all
! baseline thickness in meters
real(rprec), public :: thk_all
! orientation of turbines
integer, public :: orientation
! stagger percentage from baseline
real(rprec), public :: stag_perc
! angle from upstream (CCW from above, -x dir is zero)
real(rprec), public :: theta1_all
! angle above horizontal
real(rprec), public :: theta2_all
! thrust coefficient (default 1.33)
real(rprec), public :: Ct_prime
! Read parameters from input_turbines/param.dat
logical, public :: read_param
! Dynamically change theta1 from input_turbines/theta1.dat
logical, public :: dyn_theta1
! Dynamically change theta2 from input_turbines/theta2.dat
logical, public :: dyn_theta2
! Dynamically change Ct_prime from input_turbines/Ct_prime.dat
logical, public :: dyn_Ct_prime
! Use ADM with rotation
logical, public :: use_rotation = .false.
! Tip speed ratio for ADM with rotation
real(rprec), public :: tip_speed_ratio = 7
! disk-avg time scale in seconds (default 600)
real(rprec), public :: T_avg_dim
! filter size as multiple of grid spacing
real(rprec), public :: alpha1
real(rprec), public :: alpha2
! indicator function only includes values above this threshold
real(rprec), public :: filter_cutoff
! Correct ADM for filtered indicator function
logical, public :: adm_correction
! Number of timesteps between the output
integer, public :: tbase
integer, public :: angle_type
! The following are derived from the values above
integer :: nloc             ! total number of turbines
real(rprec) :: sx           ! spacing in the x-direction, multiple of diameter
real(rprec) :: sy           ! spacing in the y-direction

! Arrays for interpolating dynamic controls
real(rprec), dimension(:,:), allocatable :: theta1_arr
real(rprec), dimension(:), allocatable :: theta1_time
real(rprec), dimension(:,:), allocatable :: theta2_arr
real(rprec), dimension(:), allocatable :: theta2_time
real(rprec), dimension(:,:), allocatable :: Ct_prime_arr
real(rprec), dimension(:), allocatable :: Ct_prime_time

! Input files
character(:), allocatable :: input_folder
character(:), allocatable :: param_dat, theta1_dat, theta2_dat, Ct_prime_dat

! Output files
character(:), allocatable :: output_folder
character(:), allocatable :: vel_top_dat , u_d_T_dat
integer, dimension(:), allocatable :: forcing_fid

! epsilon used for disk velocity time-averaging
real(rprec) :: eps

! Commonly used indices
integer :: i, j, k, i2, j2, k2, l, s, n, m
integer :: k_start, k_end

!Store local rx, ry, rz values to an array
real(rprec), allocatable, dimension(:,:,:) ::  rx_l, ry_l, rz_l

contains

!*******************************************************************************
subroutine turbines_init()
!*******************************************************************************
!
! This subroutine creates the 'turbine' folder and starts the turbine forcing
! output files. It also creates the indicator function (Gaussian-filtered from
! binary locations - in or out) and sets values for turbine type
! (node locations, etc)
!
implicit none

real(rprec), pointer, dimension(:) :: x,y,z
character (*), parameter :: sub_name = mod_name // '.turbines_init'
integer :: fid
real(rprec) :: T_avg_dim_file, delta1, delta2
logical :: test_logical, exst
character (100) :: string1

! Set pointers
nullify(x,y,z)
x => grid % x
y => grid % y
z => grid % z

! Input/Output file names
allocate(input_folder, source = 'input_turbines/')
allocate(param_dat, source = path // input_folder // 'param.dat')
allocate(theta1_dat, source = path // input_folder // 'theta1.dat')
allocate(theta2_dat, source = path // input_folder // 'theta2.dat')
allocate(Ct_prime_dat, source = path // input_folder // 'Ct_prime.dat')
allocate(output_folder, source = 'turbine/')
allocate(vel_top_dat, source = path // output_folder // 'vel_top.dat')
allocate(u_d_T_dat, source = path // output_folder // 'u_d_T.dat')

! Allocate and initialize
nloc = num_x*num_y
nullify(wind_farm%turbine)
allocate(wind_farm%turbine(nloc))
allocate(forcing_fid(nloc))

! allocate local rx ry rz
 allocate(rx_l(ld,ny,lbz:nz), ry_l(ld,ny,lbz:nz), rz_l(ld,ny,lbz:nz))

! Create turbine directory
call system('mkdir -vp ' // path // output_folder)

! Non-dimensionalize length values by z_i
height_all = height_all / z_i
dia_all = dia_all / z_i
thk_all = thk_all / z_i

! Spacing between turbines (as multiple of mean diameter)
sx = L_x / (num_x * dia_all )
sy = L_y / (num_y * dia_all )

! Place the turbines and specify some parameters
call place_turbines

! Resize thickness to capture at least on plane of gridpoints
! and set baseline values for size
do k = 1, nloc
    wind_farm%turbine(k)%thk = max(wind_farm%turbine(k)%thk, dx * 1.01)
end do

! Specify starting and ending indices for the processor
#ifdef PPMPI
k_start = 1+coord*(nz-1)
k_end = nz-1+coord*(nz-1)
#else
k_start = 1
k_end = nz
#endif

! Read dynamic control input files
call read_control_files

!Compute a lookup table object for the indicator function
delta1 = alpha1 * sqrt(dx**2 + dy**2 + dz**2)
delta2 = alpha2 * sqrt(dx**2 + dy**2 + dz**2)
do s = 1, nloc
    call  wind_farm%turbine(s)%turb_ind_func%init(delta1, delta2,              &
            wind_farm%turbine(s)%thk, wind_farm%turbine(s)%dia)
end do

! Find turbine nodes - including filtered ind, n_hat, num_nodes, and nodes for
! each turbine. Each processor finds turbines in its domain
call turbines_nodes

! Read the time-averaged disk velocities from file if available
if (coord == 0) then
    inquire (file=u_d_T_dat, exist=exst)
    if (exst) then
        write(*,*) 'Reading from file ', trim(u_d_T_dat)
        open(newunit=fid, file=u_d_T_dat, status='unknown', form='formatted',  &
            position='rewind')
        do i=1,nloc
            read(fid,*) wind_farm%turbine(i)%u_d_T
        end do
        read(fid,*) T_avg_dim_file
        if (T_avg_dim_file /= T_avg_dim) then
            write(*,*) 'Time-averaging window does not match value in ',       &
                       trim(u_d_T_dat)
        end if
        close (fid)
    else
        write (*, *) 'File ', trim(u_d_T_dat), ' not found'
        write (*, *) 'Assuming u_d_T = -1. for all turbines'
        do k=1,nloc
            wind_farm%turbine(k)%u_d_T = -1.
        end do
    end if
end if

! Generate top of domain file
if (coord .eq. nproc-1) then
    open(newunit=fid, file=vel_top_dat, status='unknown', form='formatted',    &
        position='rewind')
    close(fid)
end if

! Generate the files for the turbine forcing output
if(coord==0) then
    do s=1,nloc
        call string_splice( string1, path // 'turbine/turbine_', s, '.dat' )
        open(newunit=forcing_fid(s), file=string1, status='unknown',           &
            form='formatted', position='append')
    end do
end if

nullify(x,y,z)

end subroutine turbines_init

!*******************************************************************************
subroutine turbines_nodes
!*******************************************************************************
!
! This subroutine locates nodes for each turbine and builds the arrays: ind,
! n_hat, num_nodes, and nodes.
!
use functions, only : cross_product, bilinear_interp
use param, only : dx, dy
use sim_param, only : detadx, detadx_dt
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_nodes'

real(rprec) :: rx,ry,rz,rxx, rzz,r,r_norm,r_disk
real(rprec), pointer :: p_xloc => null(), p_yloc => null(), p_height => null()
real(rprec), pointer :: p_dia => null(), p_thk => null()
real(rprec), pointer :: p_theta1 => null(), p_theta2 => null()
real(rprec), pointer :: p_nhat1 => null(), p_nhat2=> null(), p_nhat3 => null()
integer :: icp, jcp, kcp
integer :: imax, jmax, kmax
integer :: min_i, max_i, min_j, max_j, min_k, max_k
integer :: count_i, count_n
real(rprec), dimension(nz_tot) :: z_tot
real(rprec), dimension(3) :: temp_vec

#ifdef PPMPI
real(rprec), dimension(nloc) :: buffer_array
#endif
real(rprec), pointer, dimension(:) :: x, y, z

real(rprec) :: filt, filt_t, search_rad, filt_max
real(rprec), dimension(nloc) :: sumA, turbine_vol
logical :: test_logical

nullify(x,y,z)

x => grid % x
y => grid % y
z => grid % z

sumA = 0._rprec
rx_l = 0.0_rprec
ry_l = 0.0_rprec
rz_l = 0.0_rprec


! z_tot for total domain (since z is local to the processor)
do k = 1,nz_tot
    z_tot(k) = (k - 0.5_rprec) * dz
end do

do s=1,nloc
    count_n = 0    !used for counting nodes for each turbine
    count_i = 1    !index count - used for writing to array "nodes"

    !set pointers
    p_xloc => wind_farm%turbine(s)%xloc
    p_yloc => wind_farm%turbine(s)%yloc
    p_height => wind_farm%turbine(s)%height
    p_dia => wind_farm%turbine(s)%dia
    p_thk => wind_farm%turbine(s)%thk
    p_theta1 => wind_farm%turbine(s)%theta1
    p_theta2 => wind_farm%turbine(s)%theta2
    p_nhat1 => wind_farm%turbine(s)%nhat(1)
    p_nhat2 => wind_farm%turbine(s)%nhat(2)
    p_nhat3 => wind_farm%turbine(s)%nhat(3)

    ! Find the center of each turbine
    do k = 1,nloc
    wind_farm%turbine(k)%icp = nint(wind_farm%turbine(k)%xloc/dx)
    wind_farm%turbine(k)%jcp = nint(wind_farm%turbine(k)%yloc/dy)
    wind_farm%turbine(k)%kcp = nint(wind_farm%turbine(k)%height/dz + 0.5)

    ! Check if turbine is the current processor. MA: this part was in 
    ! turbines_init before, but was moved since kcp can be moving
    test_logical = wind_farm%turbine(k)%kcp >= k_start .and.                   &
           wind_farm%turbine(k)%kcp<=k_end
    if (test_logical) then
        wind_farm%turbine(k)%center_in_proc = .true.
    else
        wind_farm%turbine(k)%center_in_proc = .false.
    end if

    ! Make kcp the local index
    wind_farm%turbine(k)%kcp = wind_farm%turbine(k)%kcp - k_start + 1
    end do

    
    !identify "search area"
    search_rad = 0.5_rprec*p_dia + 3*max(alpha1, alpha2) * sqrt(dx**2 + dy**2 + dz**2)
    imax = min(int(search_rad/dx + 2), Nx/2)
    jmax = min(int(search_rad/dy + 2), Ny/2)
    kmax = int(search_rad/dz + 2)

    !determine unit normal vector for each turbine
    p_nhat1 = -cos(pi*p_theta1/180.)*cos(pi*p_theta2/180.)
    p_nhat2 = -sin(pi*p_theta1/180.)*cos(pi*p_theta2/180.)
    p_nhat3 = sin(pi*p_theta2/180.)

    !determine nearest (i,j,k) to turbine center
    icp = nint(p_xloc/dx)
    jcp = nint(p_yloc/dy)
    kcp = nint(p_height/dz + 0.5)

    !determine limits for checking i,j,k
    !due to spectral BCs, i and j may be < 1 or > nx,ny
    !the mod function accounts for this when these values are used
    min_i = icp-imax+1
    max_i = icp+imax
    min_j = jcp-jmax+1
    max_j = jcp+jmax
    min_k = max((kcp-kmax),1)
    max_k = min((kcp+kmax),nz_tot)
    wind_farm%turbine(s)%nodes_max(1) = min_i
    wind_farm%turbine(s)%nodes_max(2) = max_i
    wind_farm%turbine(s)%nodes_max(3) = min_j
    wind_farm%turbine(s)%nodes_max(4) = max_j
    wind_farm%turbine(s)%nodes_max(5) = min_k
    wind_farm%turbine(s)%nodes_max(6) = max_k

    ! check neighboring grid points
    ! update num_nodes, nodes, and ind for this turbine
    ! split domain between processors
    ! z(nz) and z(1) of neighboring coords match so each coord gets
    ! (local) 1 to nz-1
    wind_farm%turbine(s)%ind = 0._rprec
    wind_farm%turbine(s)%nodes = 0
    wind_farm%turbine(s)%num_nodes = 0
    count_n = 0
    count_i = 1

    ! Maximum value the filter takes (should be 1/volume)
    call wind_farm%turbine(s)%turb_ind_func%val(0._rprec, 0._rprec, filt_max, r_disk)

    ! Here the r vector is calculated for the indicator function and also the "arm"
    ! for the disk velocity calculation for dyanmic tilting and yawing. Note that
    ! the r vector for the ind is from the center and the "arm" is from the 
    ! platform base. (rx,ry,rz) are the components for r in ind. (rxx,ryy,rzz) is 
    ! for the "arm". 

    do k=k_start,k_end  !global k
        do j=min_j,max_j
            do i=min_i,max_i
                if (i<1) then
                    i2 = mod(i+nx-1,nx)+1
                    rx = (x(i2)-L_x) - p_xloc
                    rxx = (x(i2)-L_x) - wind_farm%turbine(s)%xloc_og
                elseif (i>nx) then
                    i2 = mod(i+nx-1,nx)+1
                    rx = (L_x+x(i2)) - p_xloc
                    rxx = (L_x+x(i2)) - wind_farm%turbine(s)%xloc_og
                else
                    i2 = i
                    rx = x(i) - p_xloc
                    rxx = x(i) - wind_farm%turbine(s)%xloc_og
                end if
                if (j<1) then
                    j2 = mod(j+ny-1,ny)+1
                    ry = (y(j2)-L_y) - p_yloc
                elseif (j>ny) then
                    j2 = mod(j+ny-1,ny)+1
                    ry = (L_y+y(j2)) - p_yloc
                else
                    j2 = j
                    ry = y(j) - p_yloc
                end if
                    rz = z_tot(k) - p_height
                    rzz = z_tot(k) 

                ! Local rx_l
                    rx_l(i2,j2,k-coord*(nz-1)) = rxx
                    ry_l(i2,j2,k-coord*(nz-1)) = ry
                    rz_l(i2,j2,k-coord*(nz-1)) = rzz

                r = sqrt(rx*rx + ry*ry + rz*rz)
                !length projected onto unit normal for this turbine
                r_norm = abs(rx*p_nhat1 + ry*p_nhat2 + rz*p_nhat3)
                !(remaining) length projected onto turbine disk
                r_disk = sqrt(r*r - r_norm*r_norm)
                ! get the filter value
                call wind_farm%turbine(s)%turb_ind_func%val(r_disk, r_norm, filt, filt_t)

                if ( filt > filter_cutoff * filt_max ) then
                    wind_farm%turbine(s)%ind(count_i) = filt
                    wind_farm%turbine(s)%ind_t(count_i) = filt_t
                    temp_vec(1) = rx-r_norm*p_nhat1
                    temp_vec(2) = ry-r_norm*p_nhat2
                    temp_vec(3) = rz-r_norm*p_nhat3
                    wind_farm%turbine(s)%e_theta(count_i,:) =                  &
                        cross_product(wind_farm%turbine(s)%nhat, temp_vec)
                    wind_farm%turbine(s)%e_theta(count_i,:) =                  &
                        wind_farm%turbine(s)%e_theta(count_i,:)                &
                        / sqrt(wind_farm%turbine(s)%e_theta(count_i,1)**2      &
                        + wind_farm%turbine(s)%e_theta(count_i,2)**2           &
                        + wind_farm%turbine(s)%e_theta(count_i,3)**2)
                    wind_farm%turbine(s)%nodes(count_i,1) = i2
                    wind_farm%turbine(s)%nodes(count_i,2) = j2
                    wind_farm%turbine(s)%nodes(count_i,3) = k-coord*(nz-1)!local
                    count_n = count_n + 1
                    count_i = count_i + 1
                    sumA(s) = sumA(s) + filt * dx * dy * dz
                end if
           end do
       end do
    end do
    wind_farm%turbine(s)%num_nodes = count_n

    ! Calculate turbine volume
    turbine_vol(s) = 0.25 * pi* p_dia**2 * p_thk

end do

! Sum the indicator function across all processors if using MPI
#ifdef PPMPI
buffer_array = sumA
call MPI_Allreduce(buffer_array, sumA, nloc, MPI_rprec, MPI_SUM, comm, ierr)
#endif

! Normalize the indicator function integrate to 1
do s = 1, nloc
    wind_farm%turbine(s)%ind=wind_farm%turbine(s)%ind(:)/sumA(s)
    wind_farm%turbine(s)%ind_t=wind_farm%turbine(s)%ind_t(:)/sumA(s)
end do

! Cleanup
nullify(x,y,z)

end subroutine turbines_nodes

!*******************************************************************************
subroutine turbines_forcing()
!*******************************************************************************
!
! This subroutine applies the drag-disk forcing
!
use param, only : pi, lbz, total_time, coord, wave_type
use sim_param, only : u, v, w, fxa, fya, fza, detadx, detadx_dt, eta
use functions, only : linear_interp, interp_to_uv_grid, interp_to_w_grid
use functions, only : bilinear_interp
use mpi
use mosd_wm, only : mono_wave, spectrum_wave
implicit none

character(*), parameter :: sub_name = mod_name // '.turbines_forcing'
character(*), parameter :: sub_name2 = 'angle type'

real(rprec), pointer :: p_u_d => null(), p_u_d_T => null(), p_f_n => null()
real(rprec), pointer :: p_Ct_prime => null()
integer, pointer :: p_icp => null(), p_jcp => null(), p_kcp => null()
real(rprec), pointer :: p_omegax => null(), p_omegay => null(), p_omegaz => null()
! real(rprec), pointer :: p_u1 => null()

integer :: fid

real(rprec) :: ind2
real(rprec), dimension(nloc) :: disk_avg_vel
real(rprec), dimension(nloc) :: u_vel_center, v_vel_center, w_vel_center
real(rprec), allocatable, dimension(:,:,:) ::  u_vel_disk, v_vel_disk, w_vel_disk
!real(rprec), allocatable, dimension(:,:,:) ::  u1, u2, u3
real(rprec), dimension(nloc) :: u1_center, u2_center, u3_center
real(rprec), allocatable, dimension(:,:,:) :: w_uv
real(rprec), pointer, dimension(:) :: y, z,x
real(rprec), dimension(nloc) :: buffer_array
real(rprec) :: eta_val

nullify(x,y,z)

x => grid % x
y => grid % y
z => grid % z

allocate(w_uv(ld,ny,lbz:nz))
allocate(u_vel_disk(ld,ny,lbz:nz), v_vel_disk(ld,ny,lbz:nz), w_vel_disk(ld,ny,lbz:nz))


#ifdef PPMPI
!syncing intermediate w-velocities
call mpi_sync_real_array(w, 0, MPI_SYNC_DOWNUP)
#endif

w_uv = interp_to_uv_grid(w, lbz)


! Here we select which type of angle we want for the turbines. Make sure
! theta2 is in degrees and the omega's are in radians 
select case(angle_type)

        case (0)        ! Offshore angles. Angles coming from the wave

! This is done to ensure that if the turbine is divided into
! processors, then each processor has information of the wave
#ifdef PPMPI 
        if (coord>0) then
                if (wave_type==0) then
                call mono_wave
                else 
                call spectrum_wave
                end if
        end if 
#endif
!write(*,*) 'size x', size(x)
!write(*,*) 'size y', size(y)
!write(*,*) 'detadx', size(detadx(:,:),1)
!write(*,*) 'detadx2', size(detadx(:,:),2)
        do s = 1,nloc
             wind_farm%turbine(s)%theta1 = 0.0_rprec
             wind_farm%turbine(s)%theta2 = -(bilinear_interp(x(1:nx),y(1:ny),detadx(:,:),    &
                                           wind_farm%turbine(s)%xloc_og,         &
                                           wind_farm%turbine(s)%yloc_og))*180/pi
              wind_farm%turbine(s)%omegax = 0.0_rprec
              wind_farm%turbine(s)%omegay = (bilinear_interp(x(1:nx),y(1:ny),detadx_dt(:,:), &
                                           wind_farm%turbine(s)%xloc_og,         &
                                           wind_farm%turbine(s)%yloc_og))
              wind_farm%turbine(s)%omegaz = 0.0_rprec
              wind_farm%turbine(s)%xloc = wind_farm%turbine(s)%xloc_og +         &  
                                          wind_farm%turbine(s)%height_og*        &
                                          sin(wind_farm%turbine(s)%theta2*pi/180)
              wind_farm%turbine(s)%yloc = wind_farm%turbine(s)%yloc_og
              wind_farm%turbine(s)%height = wind_farm%turbine(s)%height_og*      &
                                            cos(wind_farm%turbine(s)%theta2*pi/180)
        end do
        call turbines_nodes
        
        case (1)       ! Forced angle. Angle coming from equation
        do s = 1,nloc
             wind_farm%turbine(s)%theta1 = 0.0_rprec ! Yaw
             wind_farm%turbine(s)%theta2 = theta2_amp*sin(theta2_freq*2*pi*      & 
                                           total_time) + phi2
             wind_farm%turbine(s)%omegax = 0.0_rprec
             wind_farm%turbine(s)%omegay = (theta2_freq*2*pi*theta2_amp*          &
                                           cos(theta2_freq*2*pi*total_time))*pi/180
             wind_farm%turbine(s)%omegaz = 0.0_rprec
             wind_farm%turbine(s)%xloc = wind_farm%turbine(s)%xloc_og +            &  
                                         x_amp*sin(x_freq*2*pi*total_time)+        &
                                         wind_farm%turbine(s)%height_og*           &
                                         sin(wind_farm%turbine(s)%theta2*pi/180)
             wind_farm%turbine(s)%yloc = wind_farm%turbine(s)%yloc_og
             if (.NOT. ALLOCATED(wind_farm%turbine(s)%u1)) THEN
             !     PRINT *, 'Allocating u for turbine ', s
                 allocate(wind_farm%turbine(s)%u1(ld,ny,lbz:nz))
             end if
             !  wind_farm%turbine(s)%u2(ld,ny,lbz:nz),                              &
             !  wind_farm%turbine(s)%u3(ld,ny,lbz:nz))

             wind_farm%turbine(s)%u1(:,:,:) = x_amp*x_freq*sin(x_freq*2*pi*total_time)
             wind_farm%turbine(s)%height = wind_farm%turbine(s)%height_og*        &
                                           cos(wind_farm%turbine(s)%theta2*pi/180)
        end do
        call turbines_nodes

        case (2)       ! Dyanmic angles. Angles coming from table
        do s = 1, nloc
             if (dyn_theta1) wind_farm%turbine(s)%theta1 =                              &
                linear_interp(theta1_time, theta1_arr(s,:), total_time_dim)
             if (dyn_theta2) wind_farm%turbine(s)%theta2 =                              &
                linear_interp(theta2_time, theta2_arr(s,:), total_time_dim)
             if (dyn_Ct_prime) wind_farm%turbine(s)%Ct_prime =                          &
                linear_interp(Ct_prime_time, Ct_prime_arr(s,:), total_time_dim)
        end do
        call turbines_nodes
        
        case (3)        ! No angles
        do s = 1,nloc
             wind_farm%turbine(s)%omegax = 0.0_rprec
             wind_farm%turbine(s)%omegay = 0.0_rprec
             wind_farm%turbine(s)%omegaz = 0.0_rprec
        end do
        case default
        call error (sub_name2, 'invalid')
end select


! Here we assign the surge velocity of each turbine from input u1_freq and u1_amp  
do s = 1,nloc
    if (.NOT. ALLOCATED(wind_farm%turbine(s)%u1)) THEN
  !     PRINT *, 'Allocating u1 for turbine ', s
        allocate(wind_farm%turbine(s)%u1(ld,ny,lbz:nz))
    end if        
  !  wind_farm%turbine(s)%u2(ld,ny,lbz:nz),                              &
  !  wind_farm%turbine(s)%u3(ld,ny,lbz:nz))
    
    wind_farm%turbine(s)%u1(:,:,:) = x_amp*x_freq*sin(x_freq*2*pi*total_time)
end do
!Each processor calculates the weighted disk-averaged velocity
disk_avg_vel = 0._rprec
u_vel_center = 0._rprec
v_vel_center = 0._rprec
w_vel_center = 0._rprec
u1_center = 0._rprec
u_vel_disk = 0._rprec
v_vel_disk = 0._rprec
w_vel_disk = 0._rprec

do s=1,nloc
    !set pointers for omega
    p_omegax => wind_farm%turbine(s)%omegax
    p_omegay => wind_farm%turbine(s)%omegay
    p_omegaz => wind_farm%turbine(s)%omegaz

    !set pointers for velocity
   ! p_u1 => wind_farm%turbine(s)%u1
   ! p_u2 => wind_farm%turbine(s)%u2
   ! p_u3 => wind_farm%turbine(s)%u3
   
    do l=1,wind_farm%turbine(s)%num_nodes
        i2 = wind_farm%turbine(s)%nodes(l,1)
        j2 = wind_farm%turbine(s)%nodes(l,2)
        k2 = wind_farm%turbine(s)%nodes(l,3)

        u_vel_disk(i2,j2,k2) =  wind_farm%turbine(s)%u1(i2,j2,k2)                &
            + (p_omegay*rz_l(i2,j2,k2) - p_omegaz*ry_l(i2,j2,k2))
        v_vel_disk(i2,j2,k2) = -(p_omegax*rz_l(i2,j2,k2) - p_omegaz*rx_l(i2,j2,k2))
        w_vel_disk(i2,j2,k2) = (p_omegax*ry_l(i2,j2,k2) - p_omegay*rx_l(i2,j2,k2))

        disk_avg_vel(s) = disk_avg_vel(s)                                        &
            + dx*dy*dz*wind_farm%turbine(s)%ind(l)                               &
            * ( wind_farm%turbine(s)%nhat(1)*(u(i2,j2,k2)-u_vel_disk(i2,j2,k2))  &
            + wind_farm%turbine(s)%nhat(2)*(v(i2,j2,k2)-v_vel_disk(i2,j2,k2))    &
            + wind_farm%turbine(s)%nhat(3)*(w_uv(i2,j2,k2)-w_vel_disk(i2,j2,k2)) )
    end do

    ! Set pointers
    p_icp => wind_farm%turbine(s)%icp
    p_jcp => wind_farm%turbine(s)%jcp
    p_kcp => wind_farm%turbine(s)%kcp

    ! Calculate disk center velocity
    if (wind_farm%turbine(s)%center_in_proc) then
        u_vel_center(s) = u(p_icp, p_jcp, p_kcp)
        v_vel_center(s) = v(p_icp, p_jcp, p_kcp)
        w_vel_center(s) = w_uv(p_icp, p_jcp, p_kcp)
        u1_center(s) = wind_farm%turbine(s)%u1(p_icp, p_jcp, p_kcp)
    end if
end do

! Calculate disk velocities by summing all processors and multiplying by disk volume
#ifdef PPMPI
call MPI_Allreduce(disk_avg_vel, buffer_array, nloc, MPI_rprec, MPI_SUM, comm, ierr)
disk_avg_vel = buffer_array
call MPI_Allreduce(u_vel_center, buffer_array, nloc, MPI_rprec, MPI_SUM, comm, ierr)
u_vel_center = buffer_array
call MPI_Allreduce(v_vel_center, buffer_array, nloc, MPI_rprec, MPI_SUM, comm, ierr)
v_vel_center = buffer_array
call MPI_Allreduce(w_vel_center, buffer_array, nloc, MPI_rprec, MPI_SUM, comm, ierr)
w_vel_center = buffer_array
call MPI_Allreduce(u1_center, buffer_array, nloc, MPI_rprec, MPI_SUM, comm, ierr)
u1_center = buffer_array
#endif

! Update epsilon for the new timestep (for cfl_dt)
if (T_avg_dim > 0.) then
    eps = (dt_dim / T_avg_dim) / (1. + dt_dim / T_avg_dim)
else
    eps = 1.
end if

! Calculate and apply disk force
do s=1,nloc
    !set pointers
    p_u_d => wind_farm%turbine(s)%u_d
    p_u_d_T => wind_farm%turbine(s)%u_d_T
    p_f_n => wind_farm%turbine(s)%f_n
    p_Ct_prime => wind_farm%turbine(s)%Ct_prime

    !add this current value to the "running average" (first order filter)
    p_u_d = disk_avg_vel(s)
    if (adm_correction) then
        p_u_d = p_u_d /(1 + 0.25_rprec                                         &
            * (1-wind_farm%turbine(s)%turb_ind_func%M)*p_Ct_prime)
    end if
    p_u_d_T = (1.-eps)*p_u_d_T + eps*p_u_d

    !calculate total thrust force for each turbine  (per unit mass)
    !force is normal to the surface (calc from u_d_T, normal to surface)
    !write force to array that will be transferred via MPI
    p_f_n = -0.5*p_Ct_prime*abs(p_u_d_T)*p_u_d_T*0.25*pi*wind_farm%turbine(s)%dia**2

    !write values to file
    if (modulo (jt_total, tbase) == 0 .and. coord == 0) then
!        if (angle_type==0) then
!                eta_val = (bilinear_interp(x(1:nx),y(1:ny),eta(:,:),                    &
!                                           wind_farm%turbine(s)%xloc_og,    &
!                                         wind_farm%turbine(s)%yloc_og))
!        else
!                eta_val = 0.0 
!        end if
        write( forcing_fid(s), *) total_time_dim, wind_farm%turbine(s)%xloc,          & 
            u_vel_center(s), v_vel_center(s), w_vel_center(s), -p_u_d,                &
            -p_u_d_T, wind_farm%turbine(s)%theta1, wind_farm%turbine(s)%theta2,       &
            u1_center(s), p_Ct_prime, jt_total                   
!    eta_val
    end if


    do l=1,wind_farm%turbine(s)%num_nodes
        i2 = wind_farm%turbine(s)%nodes(l,1)
        j2 = wind_farm%turbine(s)%nodes(l,2)
        k2 = wind_farm%turbine(s)%nodes(l,3)
        ind2 = wind_farm%turbine(s)%ind(l)
        fxa(i2,j2,k2) = p_f_n*wind_farm%turbine(s)%nhat(1)*ind2
        fya(i2,j2,k2) = p_f_n*wind_farm%turbine(s)%nhat(2)*ind2
        fza(i2,j2,k2) = p_f_n*wind_farm%turbine(s)%nhat(3)*ind2
        if (use_rotation) then
            ind2 = wind_farm%turbine(s)%ind_t(l)
            fxa(i2,j2,k2) = fxa(i2,j2,k2)                                      &
                + p_f_n*wind_farm%turbine(s)%e_theta(l,1)*ind2/tip_speed_ratio
            fya(i2,j2,k2) = fya(i2,j2,k2)                                      &
                + p_f_n*wind_farm%turbine(s)%e_theta(l,2)*ind2/tip_speed_ratio
            fza(i2,j2,k2) = fza(i2,j2,k2)                                      &
                + p_f_n*wind_farm%turbine(s)%e_theta(l,3)*ind2/tip_speed_ratio
        end if
    end do
end do

! Interpolate force onto the w grid
call mpi_sync_real_array( fxa(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( fya(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( fza(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
fza = interp_to_w_grid(fza,lbz)

!spatially average velocity at the top of the domain and write to file
if (coord .eq. nproc-1) then
    open(newunit=fid, file=vel_top_dat, status='unknown', form='formatted',    &
        action='write', position='append')
    write(fid,*) total_time, sum(u(:,:,nz-1))/(nx*ny)
    close(fid)
end if

! Cleanup
deallocate(w_uv)
nullify(x,y,z)
nullify(p_icp, p_jcp, p_kcp)
nullify(p_omegax,p_omegay, p_omegaz)

end subroutine turbines_forcing

!*******************************************************************************
subroutine turbines_finalize ()
!*******************************************************************************
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_finalize'

!write disk-averaged velocity to file along with T_avg_dim
!useful if simulation has multiple runs   >> may not make a large difference
call turbines_checkpoint

!deallocate
deallocate(wind_farm%turbine)

end subroutine turbines_finalize

!*******************************************************************************
subroutine turbines_checkpoint ()
!*******************************************************************************
!
!
!
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_checkpoint'
integer :: fid

!write disk-averaged velocity to file along with T_avg_dim
!useful if simulation has multiple runs   >> may not make a large difference
if (coord == 0) then
    open(newunit=fid, file=u_d_T_dat, status='unknown', form='formatted',      &
        position='rewind')
    do i=1,nloc
        write(fid,*) wind_farm%turbine(i)%u_d_T
    end do
    write(fid,*) T_avg_dim
    close (fid)
end if

end subroutine turbines_checkpoint

!*******************************************************************************
subroutine turbine_vel_init(zo_high)
!*******************************************************************************
!
! called from ic.f90 if initu, lbc_mom==1, S_FLAG are all false.
! this accounts for the turbines when creating the initial velocity profile.
!
use param, only: zo
implicit none
character (*), parameter :: sub_name = mod_name // '.turbine_vel_init'

real(rprec), intent(inout) :: zo_high
real(rprec) :: cft, nu_w, exp_KE, induction_factor, Ct_noprime

! Convert Ct' to Ct
! a = Ct'/(4+Ct'), Ct = 4a(1-a)
induction_factor = Ct_prime / (4._rprec + Ct_prime)
Ct_noprime = 4*(induction_factor) * (1 - induction_factor)

! friction coefficient, cft
cft = pi*Ct_noprime/(4.*sx*sy)

!wake viscosity
nu_w = 28.*sqrt(0.5*cft)

!turbine friction height, Calaf, Phys. Fluids 22, 2010
zo_high = height_all*(1.+0.5*dia_all/height_all)**(nu_w/(1.+nu_w))* &
  exp(-1.*(0.5*cft/(vonk**2) + (log(height_all/zo* &
  (1.-0.5*dia_all/height_all)**(nu_w/(1.+nu_w))) )**(-2) )**(-0.5) )

exp_KE =  0.5*(log(0.45/zo_high)/0.4)**2

if(.false.) then
    write(*,*) 'sx,sy,cft: ',sx,sy,cft
    write(*,*) 'nu_w: ',nu_w
    write(*,*) 'zo_high: ',zo_high
    write(*,*) 'approx expected KE: ', exp_KE
end if
end subroutine turbine_vel_init

!*******************************************************************************
subroutine place_turbines
!*******************************************************************************
!
! This subroutine places the turbines on the domain. It also sets the values for
! each individual turbine. After the subroutine is called, the following values
! are set for each turbine in wind_farm: xloc, yloc, height, dia, thk, theta1,
! theta2, and Ct_prime.
!
use param, only: pi, z_i
use messages
implicit none

character(*), parameter :: sub_name = mod_name // '.place_turbines'

real(rprec) :: sxx, syy, shift_base, const
real(rprec) :: dummy, dummy2
logical :: exst
integer :: fid

! Read parameters from file if needed
if (read_param) then
    ! Check if file exists and open
    inquire (file = param_dat, exist = exst)
    if (.not. exst) then
        call error (sub_name, 'file ' // param_dat // 'does not exist')
    end if

    ! Check that there are enough lines from which to read data
    nloc = count_lines(param_dat)
    if (nloc < num_x*num_y) then
        nloc = num_x*num_y
        call error(sub_name, param_dat // 'must have num_x*num_y lines')
    else if (nloc > num_x*num_y) then
    call warn(sub_name, param_dat // ' has more than num_x*num_y lines. '  &
                  // 'Only reading first num_x*num_y lines')
    end if

    ! Read from parameters file, which should be in this format:
    ! xloc [meters], yloc [meters], height [meters], dia [meters], thk [meters],
    ! theta1 [degrees], theta2 [degrees], Ct_prime [-]
    write(*,*) "Reading from", param_dat
    open(newunit=fid, file=param_dat, status='unknown', form='formatted',      &
        position='rewind')
    do k = 1, nloc
        read(fid,*) wind_farm%turbine(k)%xloc, wind_farm%turbine(k)%yloc,      &
            wind_farm%turbine(k)%height, wind_farm%turbine(k)%dia,             &
            wind_farm%turbine(k)%thk, wind_farm%turbine(k)%theta1,             &
            wind_farm%turbine(k)%theta2, wind_farm%turbine(k)%Ct_prime
    end do
    close(fid)

    ! Make lengths dimensionless
    do k = 1, nloc
        wind_farm%turbine(k)%xloc = wind_farm%turbine(k)%xloc / z_i
        wind_farm%turbine(k)%yloc = wind_farm%turbine(k)%yloc / z_i
        wind_farm%turbine(k)%height = wind_farm%turbine(k)%height / z_i
        wind_farm%turbine(k)%dia = wind_farm%turbine(k)%dia / z_i
        wind_farm%turbine(k)%thk = wind_farm%turbine(k)%thk / z_i
    end do

    ! flag for setting the prescribed turbine location as the orginal (OG)
    ! location. this is for turbines dynamically tilting, yawing, etc
    if (angle_type < 3) then

    wind_farm%turbine(:)%xloc_og = wind_farm%turbine(:)%xloc
    wind_farm%turbine(:)%yloc_og = wind_farm%turbine(:)%yloc
    wind_farm%turbine(:)%height_og = wind_farm%turbine(:)%height

    end if

else
    ! Set values for each turbine based on values in input file
    wind_farm%turbine(:)%height = height_all
    wind_farm%turbine(:)%dia = dia_all
    wind_farm%turbine(:)%thk = thk_all
    wind_farm%turbine(:)%theta1 = theta1_all
    wind_farm%turbine(:)%theta2 = theta2_all
    wind_farm%turbine(:)%Ct_prime = Ct_prime

    ! Set baseline locations (evenly spaced, not staggered aka aligned)
    k = 1
    sxx = sx * dia_all  ! x-spacing with units to match those of L_x
    syy = sy * dia_all  ! y-spacing
    do i = 1,num_x
        do j = 1,num_y
            wind_farm%turbine(k)%xloc = sxx*real(2*i-1)/2
            wind_farm%turbine(k)%yloc = syy*real(2*j-1)/2
            k = k + 1
        end do
    end do

    ! Place turbines based on orientation flag
    ! This will shift the placement relative to the baseline locations abive
    select case (orientation)
        ! Evenly-spaced, not staggered
        case (1)

        ! Evenly-spaced, horizontally staggered only
        ! Shift each row according to stag_perc
        case (2)
            do i = 2, num_x
                do k = 1+num_y*(i-1), num_y*i
                    shift_base = syy * stag_perc/100.
                    wind_farm%turbine(k)%yloc = mod( wind_farm%turbine(k)%yloc &
                                                    + (i-1)*shift_base , L_y )
                end do
            end do

        ! Evenly-spaced, only vertically staggered (by rows)
        case (3)
            ! Make even rows taller
            do i = 2, num_x, 2
                do k = 1+num_y*(i-1), num_y*i
                    wind_farm%turbine(k)%height = height_all*(1.+stag_perc/100.)
                end do
            end do
            ! Make odd rows shorter
            do i = 1, num_x, 2
                do k = 1+num_y*(i-1), num_y*i
                    wind_farm%turbine(k)%height = height_all*(1.-stag_perc/100.)
                end do
            end do

        ! Evenly-spaced, only vertically staggered, checkerboard pattern
        case (4)
            k = 1
            do i = 1, num_x
                do j = 1, num_y
                    ! this should alternate between 1, -1
                    const = 2.*mod(real(i+j),2.)-1.
                    wind_farm%turbine(k)%height = height_all                   &
                                                  *(1.+const*stag_perc/100.)
                    k = k + 1
                end do
            end do

        ! Aligned, but shifted forward for efficient use of simulation space
        ! during CPS runs
        case (5)
        ! Shift in spanwise direction: Note that stag_perc is now used
            k=1
            dummy=stag_perc                                                    &
                  *(wind_farm%turbine(2)%yloc - wind_farm%turbine(1)%yloc)
            do i = 1, num_x
                do j = 1, num_y
                    dummy2=dummy*(i-1)
                    wind_farm%turbine(k)%yloc=mod( wind_farm%turbine(k)%yloc   &
                                                  + dummy2,L_y)
                    k=k+1
                end do
            end do

        case default
            call error (sub_name, 'invalid orientation')
    end select

    ! flag for setting the prescribed turbine location as the orginal (OG)
    ! location. this is for turbines dynamically tilting, yawing, etc
    if (angle_type < 3) then

    wind_farm%turbine(:)%xloc_og = wind_farm%turbine(:)%xloc
    wind_farm%turbine(:)%yloc_og = wind_farm%turbine(:)%yloc
    wind_farm%turbine(:)%height_og = wind_farm%turbine(:)%height

    end if

end if

end subroutine place_turbines

!*******************************************************************************
subroutine read_control_files
!*******************************************************************************
!
! This subroutine reads the input files for dynamic controls with theta1,
! theta2, and Ct_prime. This is calles from turbines_init.
!
use param, only: pi
use messages
implicit none

character(*), parameter :: sub_name = mod_name // '.place_turbines'

integer :: fid, i, num_t

! Read the theta1 input data
if (dyn_theta1) then
    ! Count number of entries and allocate
    num_t = count_lines(theta1_dat)
    allocate( theta1_time(num_t) )
    allocate( theta1_arr(nloc, num_t) )

    ! Read values from file
    open(newunit=fid, file=theta1_dat, status='unknown', form='formatted',     &
        position='rewind')
    do i = 1, num_t
        read(fid,*) theta1_time(i), theta1_arr(:,i)
    end do
end if

! Read the theta2 input data
if (dyn_theta2) then
    ! Count number of entries and allocate
    num_t = count_lines(theta2_dat)
    allocate( theta2_time(num_t) )
    allocate( theta2_arr(nloc, num_t) )

    ! Read values from file
    open(newunit=fid, file=theta2_dat, status='unknown', form='formatted',     &
        position='rewind')
    do i = 1, num_t
        read(fid,*) theta2_time(i), theta2_arr(:,i)
    end do
end if

! Read the Ct_prime input data
if (dyn_Ct_prime) then
    ! Count number of entries and allocate
    num_t = count_lines(Ct_prime_dat)
    allocate( Ct_prime_time(num_t) )
    allocate( Ct_prime_arr(nloc, num_t) )

    ! Read values from file
    open(newunit=fid, file=Ct_prime_dat, status='unknown', form='formatted',   &
        position='rewind')
    do i = 1, num_t
        read(fid,*) Ct_prime_time(i), Ct_prime_arr(:,i)
    end do
end if

end subroutine read_control_files

end module turbines
