program AmplitudeSourceLocation_DoubleDifference
  !!Relative amplitude source location (Double-differnce of amplitude)
  !!using depth-dependent 1D velocity structure, 3D heterogeneous attenuation structure
  !!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2020
  !!License  : MIT License (https://opensource.org/licenses/MIT)

  !!-DWITHOUT_ERROR: do not calculate estimation errors

  use nrtype,               only : fp, dp
  use constants,            only : rad2deg, deg2rad, pi, r_earth
  use rayshooting,          only : rayshooting3D
  use set_velocity_model,   only : set_velocity
  use linear_interpolation, only : linear_interpolation_1d, linear_interpolation_2d, linear_interpolation_3d, &
  &                                block_interpolation_3d
  use greatcircle,          only : greatcircle_dist, latctog, latgtoc
  use grdfile_io,           only : read_grdfile_2d
  !$ use omp_lib

#if defined (RAY_BENDING)
  use raybending,           only : pseudobending3D
#endif

#if defined (MKL)
  use lapack95
#else
  use f95_lapack
#endif

  use lsqr_kinds
  use lsqr_module,          only : lsqr_solver_ez

  implicit none

  integer,            parameter :: wavetype = 2          !!1 for P-wave, 2 for S-wave
  real(kind = fp),    parameter :: snratio_accept = 0.1_fp
  real(kind = fp),    parameter :: interevent_dist_max = 0.5_fp
  real(kind = fp),    parameter :: delta_residual_max = 1.0e-10_fp
  integer,            parameter :: nsta_use_min = 4
  integer,            parameter :: nconstraint = 4
  integer,            parameter :: loop_count_max = 10
  !!Range for velocity and attenuation structure
  !!whole Japan
  real(kind = fp),    parameter :: lon_str_w = 122.0_fp, lon_str_e = 150.0_fp
  real(kind = fp),    parameter :: lat_str_s = 23.0_fp, lat_str_n = 46.0_fp
  real(kind = fp),    parameter :: z_str_min = -3.0_fp, z_str_max = 50.0_fp
  real(kind = fp),    parameter :: dlon_str = 0.1_fp, dlat_str = 0.1_fp, dz_str = 0.1_fp
  integer,            parameter :: nlon_str = int((lon_str_e - lon_str_w) / dlon_str) + 2
  integer,            parameter :: nlat_str = int((lat_str_n - lat_str_s) / dlat_str) + 2
  integer,            parameter :: nz_str = int((z_str_max - z_str_min) / dz_str) + 2
  !!Ray shooting
  real(kind = fp),    parameter :: dvdlon = 0.0_fp, dvdlat = 0.0_fp         !!assume 1D structure
  integer,            parameter :: ninc_angle = 180                         !!grid search in incident angle
  integer,            parameter :: nrayshoot = 2                            !!number of grid search
  real(kind = fp),    parameter :: time_step = 0.005_fp
  real(kind = fp),    parameter :: rayshoot_dist_thr = 0.05_fp

  real(kind = fp),    allocatable :: stlon(:), stlat(:), stdp(:), ttime_cor(:, :)
  integer,            allocatable :: event_index(:), event_index_rev(:)
  character(len = 6), allocatable :: stname(:)
  logical,            allocatable :: stuse_flag(:), evflag(:, :), obsamp_flag(:, :)
  character(len = 10)             :: freq_t

  real(kind = fp),    parameter :: alt_to_depth = -1.0e-3_fp
  real(kind = dp),    parameter :: huge = 1.0e+6_dp

  real(kind = fp)               :: velocity(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                qinv(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                val_1d(1 : 2), val_2d(1 : 2, 1 : 2), val_3d(1 : 2, 1 : 2, 1 : 2), &
  &                                xgrid(1 : 2), ygrid(1 : 2), zgrid(1 : 2), inc_angle_ini_min(0 : nrayshoot), &
  &                                normal_vector_j(1 : 3), normal_vector_k(1 : 3)
  real(kind = dp),  allocatable :: topography(:, :), lon_topo(:), lat_topo(:)
  real(kind = fp),  allocatable :: obsamp(:, :), calamp(:, :), obsamp_noise(:), hypodist(:, :), ray_azinc(:, :, :), &
  &                                evlon(:, :), evlat(:, :), evdp(:, :), evamp(:, :), &
  &                                qinv_interpolate(:), velocity_interpolate(:), residual_tmp(:)
  real(kind = fp)               :: evlat_tmp, epdist, epdelta, az_ini, dinc_angle_org, dinc_angle, inc_angle_ini, &
  &                                lon_tmp, lat_tmp, depth_tmp, az_tmp, inc_angle_tmp, dist_tmp, ttime_tmp, width_tmp, &
  &                                dist_min, ttime_min, width_min, &
  &                                dvdz, lon_new, lat_new, depth_new, az_new, inc_angle_new, matrix_const_j, &
  &                                matrix_const_k, lon_min, lat_min, depth_min, delta_depth, delta_lon, delta_lat, &
  &                                residual_sum, residual_old, sigma_lon, sigma_lat, sigma_depth, sigma_amp, &
  &                                siteamp_tmp, interevent_dist
  real(kind = dp)               :: topography_interpolate, dlon_topo, dlat_topo, freq
  integer                       :: nlon_topo, nlat_topo, nsta, nevent, lon_index, lat_index, z_index, &
  &                                i, j, k, ii, jj, kk, loop_count, icount, nsta_use, ios, &
  &                                obsvector_count, irow_count, icol_count, nnonzero_elem, nonzero_elem_count, &
  &                                nevent_est, nobsamp_ratio, mainloop_count, event_count
  character(len = 129)          :: topo_grd, station_param, event_initloc_param, event_amp_param, resultfile
  character(len = 20)           :: cfmt, nsta_c

#if defined (RAY_BENDING)
  integer,                parameter :: ndiv_raypath = 10
  integer,                parameter :: nraypath_ini = 4
  real(kind = fp)                   :: raypath_lon((nraypath_ini - 1) * 2 ** ndiv_raypath + 1), &
  &                                    raypath_lat((nraypath_ini - 1) * 2 ** ndiv_raypath + 1), &
  &                                    raypath_dep((nraypath_ini - 1) * 2 ** ndiv_raypath + 1)
  integer                           :: nraypath
#endif

  !!LSQR variables
  type(lsqr_solver_ez)              :: solver
  real(kind = fp),      allocatable :: modelvector(:), obsvector(:), nonzero_elem(:)
  integer,              allocatable :: irow(:), icol(:)
  integer                           :: istop

  !!OpenMP variable
  !$ integer                    :: omp_thread

  icount = iargc()

  if(icount .ne. 6) then
    write(0, '(a)', advance="no") "usage: ./asl_dd "
    write(0, '(a)', advance="no") "(topography_grd) (station_param_file) (event_initial_location_file) (event_amplitude_file) "
    write(0, '(a)')               "(frequency) (result_file)"
    error stop
  endif
  call getarg(1, topo_grd)
  call getarg(2, station_param)
  call getarg(3, event_initloc_param)
  call getarg(4, event_amp_param)
  call getarg(5, freq_t); read(freq_t, *) freq
  call getarg(6, resultfile)

  !!read topography file (netcdf grd format)
  call read_grdfile_2d(topo_grd, lon_topo, lat_topo, topography)
  nlon_topo = ubound(lon_topo, 1)
  nlat_topo = ubound(lat_topo, 1)
  dlon_topo = lon_topo(2) - lon_topo(1)
  dlat_topo = lat_topo(2) - lat_topo(1)
  topography(1 : nlon_topo, 1 : nlat_topo) = topography(1 : nlon_topo, 1 : nlat_topo) * alt_to_depth

  !!read station parameter
  open(unit = 10, file = station_param)
  nsta = 0
  do
    read(10, *, iostat = ios)
    if(ios .ne. 0) exit
    nsta = nsta + 1
  enddo
  rewind(10)
  allocate(stlon(1 : nsta), stlat(1 : nsta), stdp(1 : nsta), stname(1 : nsta), ttime_cor(1 : nsta, 1 : 2), &
  &        stuse_flag(1 : nsta), obsamp_noise(1 : nsta))
  do i = 1, nsta
    read(10, *) stlon(i), stlat(i), stdp(i), stname(i), stuse_flag(i), &
    &           ttime_cor(i, 1), ttime_cor(i, 2), siteamp_tmp, obsamp_noise(i)
    write(0, '(a, i0, a, f9.4, a, f8.4, a, f6.3, 1x, a7, l2)') &
    &     "station(", i, ") lon(deg) = ", stlon(i), " lat(deg) = ", stlat(i), " depth(km) = ", stdp(i), &
    &     trim(stname(i)), stuse_flag(i)
  enddo
  close(10)
  if(.not. (nsta .ge. nsta_use_min .and. nsta_use .ge. nsta_use_min)) then
    close(10)
    write(0, '(a, i0)') "Number of stations for calculation should be larger than ", nsta_use_min
    error stop
  endif


  !!read event initial locations
  open(unit = 10, file = event_initloc_param)
  nevent = 0
  read(10, *)
  do
    read(10, *, iostat = ios)
    if(ios .ne. 0) exit
    nevent = nevent + 1
  enddo    
  rewind(10)
  write(0, '(a, i0)') "nevent = ", nevent
  allocate(evlon(1 : nevent, 0 : nevent + 1), evlat(1 : nevent, 0 : nevent + 1), evdp(1 : nevent, 0 : nevent + 1), &
  &        evamp(1 : nevent, 0 : nevent + 1), evflag(1 : nevent, 0 : nevent + 1))
  allocate(hypodist(1 : nsta, 1 : nevent), ray_azinc(1 : 2, 1 : nsta, 1 : nevent))
  allocate(obsamp_flag(1 : nsta, 1 : nevent))
  read(10, *)
  do i = 1, nevent
    read(10, *) evlon(i, 0), evlat(i, 0), evdp(i, 0), evamp(i, 0)
    evlon (i, 1 : nevent + 1) = evlon(i, 0)
    evlat (i, 1 : nevent + 1) = evlat(i, 0)
    evdp  (i, 1 : nevent + 1) = evdp (i, 0)
    evamp (i, 1 : nevent + 1) = evamp(i, 0)
    evflag(i, 0 : nevent + 1) = .true.
  enddo
  close(10)
    
  open(unit = 10, file = event_amp_param)
  allocate(obsamp(1 : nsta, 1 : nevent))
  read(10, *)
  do j = 1, nevent
    obsamp_flag(1 : nsta, j) = .true.
    read(10, *) (obsamp(i, j), i = 1, nsta)
    nsta_use = 0
    do i = 1, nsta
      if(.not. (stuse_flag(i) .eqv. .true. &
      &  .and.  obsamp(i, j)  .gt.  0.0_fp &
      &  .and.  obsamp(i, j) / obsamp_noise(i) .gt. snratio_accept)) then
        obsamp_flag(i, j) = .false.
        cycle
      endif
      nsta_use = nsta_use + 1
    enddo
    if(nsta_use .lt. nsta_use_min) then
      obsamp_flag(1 : nsta, j)  = .false.
      evflag(j, 0 : nevent + 1) = .false.
    endif
  enddo
  close(10)

  !!set velocity/attenuation structure
  call set_velocity(z_str_min, dz_str, velocity, qinv)

  !!Main loop: estimate \Delta_x, \Delta_y, \Delta_z, \Delta_amp until convergence
  main_loop: do mainloop_count = 1, nevent + 1
    if(mainloop_count .gt. 1) evflag(mainloop_count - 1, mainloop_count) = .false.
    residual_old = huge
    loop_count = 0
    iter_loop: do
      if(loop_count .gt. loop_count_max) exit iter_loop
      loop_count = loop_count + 1   
      !!calculate ray length, pulse width, unit vector of ray incident
      write(0, '(a)') "calculate ray length, pulse width, and ray incident vector for each event"

      nevent_est = 0
      !$omp parallel default(none), &
      !$omp&         shared(obsamp_flag, nsta, evlon, evlat, evdp, stlon, stlat, stdp, hypodist, ray_azinc, &
      !$omp&                lon_topo, lat_topo, dlon_topo, dlat_topo, topography, nlon_topo, nlat_topo, &
      !$omp&                calamp, velocity, qinv, evflag, velocity_interpolate, qinv_interpolate), &
      !$omp&         private(epdist, az_ini, epdelta, lon_index, lat_index, z_index, dinc_angle_org, dinc_angle, &
      !$omp&                 inc_angle_ini_min, inc_angle_ini, lon_tmp, lat_tmp, depth_tmp, az_tmp, dist_min, inc_angle_tmp, &
      !$omp&                 xgrid, ygrid, zgrid, val_1d, val_2d, val_3d, topography_interpolate, ttime_tmp, width_tmp, &
      !$omp&                 dist_tmp, lon_min, lat_min, depth_min, depth_max, depth_max_tmp, ttime_min, width_min, &
      !$omp&                 dvdz, &
#if defined (RAY_BENDING)
      !$omp&                 raypath_lon, raypath_lat, raypath_dep, nraypath, &
#endif
      !$omp&                 lon_new, lat_new, depth_new, az_new, inc_angle_new, i, ii, jj, kk, omp_thread), &
      !$omp&         reduction(+: nevent_est)
    
      !$ omp_thread = omp_get_thread_num()

      !!calculate traveltime from each event location to stations
      !$omp do
      event_loop: do j = 1, nevent
        !!check whether the depth of each event location is lower than the topo
        lon_index = int((evlon(j, mainloop_count) - lon_topo(1)) / dlon_topo) + 1
        lat_index = int((evlat(j, mainloop_count) - lat_topo(1)) / dlat_topo) + 1
        xgrid(1 : 2) = [lon_topo(lon_index), lon_topo(lon_index + 1)]
        ygrid(1 : 2) = [lat_topo(lat_index), lat_topo(lat_index + 1)]
        val_2d(1 : 2, 1 : 2) = topography(lon_index : lon_index + 1, lat_index : lat_index + 1)
        call linear_interpolation_2d(evlon(i, mainloop_count), evlat(i, mainloop_count), xgrid, ygrid, val_2d, &
        &                            topography_interpolate)

        if(evdp(j, mainloop_count) .lt. topography_interpolate) evflag(j, mainloop_count) = .false.
        if(evflag(j, mainloop_count) .eqv. .false.) cycle

        nevent_est = nevent_est + 1
  
        station_loop: do i = 1, nsta
          if(obsamp_flag(i, j) .eqv. .false.) cycle
          print *, stlon(i), stlat(i)
          !!calculate azimuth and hypocentral distance
          call greatcircle_dist(evlat(j, mainloop_count), evlon(j, mainloop_count), stlat(i), stlon(i), &
          &                     distance = epdist, azimuth = az_ini, delta_out = epdelta)
          lon_index = int((stlon(i) - lon_str_w) / dlon_str) + 1
          lat_index = int((stlat(i) - lat_str_s) / dlat_str) + 1
          z_index   = int((stdp(i) - z_str_min) / dz_str) + 1
          !print *, lon_sta(jj), lon_w + real(lon_index - 1) * dlon
          !print *, lat_sta(jj), lat_s + real(lat_index - 1) * dlat
          hypodist(i, j) = sqrt((r_earth - evdp(j, mainloop_count)) ** 2 + (r_earth - stdp(i)) ** 2 &
          &                    - 2.0_fp * (r_earth - evdp(j, mainloop_count)) * (r_earth - stdp(i)) * cos(epdelta))
          ray_azinc(1, i, j) = az_ini

#if defined (V_CONST)
          !!homogeneous structure: ray incident angle is calculated using cosine function (assuming cartesian coordinate)
          ray_azinc(2, i, j) = acos(epdist / hypodist(i, j)) + pi * 0.5_fp
          ttime_min = hypodist(i, j) / velocity(lon_index, lat_index, z_index, wavetype)
          width_min = ttime_min * qinv(lon_index, lat_index, z_index, wavetype)
#else

#if defined (RAY_BENDING)
          !!do ray tracing with pseudobending scheme
          raypath_lon(1) = evlon(j, mainloop_count)
          raypath_lat(1) = evlat(j, mainloop_count)
          raypath_dep(1) = evdp(j, mainloop_count)
          raypath_lon(nraypath_ini) = stlon(i)
          raypath_lat(nraypath_ini) = stlat(i)
          raypath_dep(nraypath_ini) = stdp(i)
          do jj = 2, nraypath_ini - 1
            raypath_lon(jj) = raypath_lon(jj - 1) + (raypath_lon(nraypath_ini) - raypath_lon(1)) / real(nraypath_ini, kind = fp)
            raypath_lat(jj) = raypath_lat(jj - 1) + (raypath_lat(nraypath_ini) - raypath_lat(1)) / real(nraypath_ini, kind = fp)
            raypath_dep(jj) = raypath_dep(jj - 1) + (raypath_dep(nraypath_ini) - raypath_dep(1)) / real(nraypath_ini, kind = fp)
          enddo
          nraypath = nraypath_ini
          call pseudobending3D(raypath_lon, raypath_lat, raypath_dep, nraypath, ndiv_raypath, velocity(:, :, :, wavetype), &
          &                    lon_str_w, lat_str_s, z_str_min, dlon_str, dlat_str, dz_str, ttime_min, &
          &                    qinv = qinv(:, :, :, wavetype), lon_w_qinv = lon_str_w, lat_s_qinv = lat_str_s, &
          &                    dep_min_qinv = z_str_min, dlon_qinv = dlon_str, dlat_qinv = dlat_str, ddep_qinv = dz_str, &
          &                    pulsewidth = width_min, ray_az = ray_azinc(1, i, j), ray_incangle = ray_azinc(2, i, j))

          write(0, '(a, 4(f8.4, 1x))') "event lon, lat, depth, az_ini = ", &
          &                            evlon(j, mainloop_count), evlat(j, mainloop_count), &
          &                            evdp(j, mainloop_count), az_ini * rad2deg
          write(0, '(a, 3(f8.4, 1x))') "station lon, lat, depth = ", stlon(i), stlat(i), stdp(i)
          write(0, '(a, 2(f8.4, 1x))') "ray azimuth and inc_angle (deg) = ", &
          &                             ray_azinc(1, i, j) * rad2deg, ray_azinc(2, i, j) * rad2deg
          write(0, '(a, f5.2)') "traveltime (s) = ", ttime_min

#else /* -DRAY_BENDING */

          !!do ray shooting
          dist_min = real(huge, kind = fp)
          incangle_loop2: do jj = 1, nrayshoot
            if(jj .eq. 1) then
              dinc_angle_org = pi * 0.5_fp
            else
              dinc_angle_org = dinc_angle
            endif
            dinc_angle = 2.0_fp * dinc_angle_org / real(ninc_angle, kind = fp)
            inc_angle_ini_min(0) = dinc_angle_org
            !print *, "dinc_angle = ", dinc_angle * rad2deg, inc_angle_ini_min(kk - 1) * rad2deg

            incangle_loop: do ii = 1, ninc_angle
              inc_angle_ini = (inc_angle_ini_min(jj - 1) - dinc_angle_org) + real(ii, kind = fp) * dinc_angle
              !print '(2(i0, 1x), a, e15.7)', ii, kk, "inc_angle_ini = ", inc_angle_ini * rad2deg
                   
              lon_tmp = evlon(j, mainloop_count)
              lat_tmp = evlat(j, mainloop_count)
              depth_tmp = evdp(j, mainloop_count)
              az_tmp = az_ini
              inc_angle_tmp = inc_angle_ini

              !!loop until ray arrives at surface/boundary
              ttime_tmp = 0.0_fp
              width_tmp = 0.0_fp
              shooting_loop: do

                !!exit if ray approaches to the surface
                lon_index = int((lon_tmp - lon_topo(1)) / dlon_topo) + 1
                lat_index = int((lat_tmp - lat_topo(1)) / dlat_topo) + 1
                !!exit if ray approaches to the boundary of the topography array
                if(lon_index .lt. 1 .or. lon_index .gt. nlon_topo - 1      &
                &  .or. lat_index .lt. 1 .or. lat_index .gt. nlat_topo - 1) exit shooting_loop

                xgrid(1 : 2) = [lon_topo(lon_index), lon_topo(lon_index + 1)]
                ygrid(1 : 2) = [lat_topo(lat_index), lat_topo(lat_index + 1)]
                val_2d(1 : 2, 1 : 2) = topography(lon_index : lon_index + 1, lat_index : lat_index + 1)
                call linear_interpolation_2d(lon_tmp, lat_tmp, xgrid, ygrid, val_2d, topography_interpolate)
                !print '(a, 3(f9.4, 1x))', "ray surface arrived, lon/lat = ", lon_tmp, lat_tmp, depth_tmp
                if(depth_tmp .lt. topography_interpolate) exit shooting_loop

                lon_index = int((lon_tmp - lon_str_w) / dlon_str) + 1
                lat_index = int((lat_tmp - lat_str_s) / dlat_str) + 1
                z_index   = int((depth_tmp - z_str_min) / dz_str) + 1

                !!exit if ray approaches to the boundary
                if(lon_index .lt. 1 .or. lon_index .gt. nlon_str - 1 .or. &
                &  lat_index .lt. 1 .or. lat_index .gt. nlat_str - 1 .or. &
                &  z_index   .lt. 1 .or. z_index   .gt. nz_str   - 1) exit shooting_loop

                xgrid(1) = lon_str_w + real(lon_index - 1, kind = fp) * dlon_str; xgrid(2) = xgrid(1) + dlon_str
                ygrid(1) = lat_str_s + real(lat_index - 1, kind = fp) * dlat_str; ygrid(2) = ygrid(1) + dlat_str
                zgrid(1) = z_str_min + real(z_index   - 1, kind = fp) * dz_str  ; zgrid(2) = zgrid(1) + dz_str

                !!calculate distance between ray and station
                call greatcircle_dist(lat_tmp, lon_tmp, stlat(i), stlon(i), delta_out = epdelta)
                dist_tmp = sqrt((r_earth - depth_tmp) ** 2 + (r_earth - stdp(i)) ** 2 &
                &              - 2.0_fp * (r_earth - depth_tmp) * (r_earth - stdp(i)) * cos(epdelta))
                !print *, real(ii, kind = fp) * dinc_angle, lat_tmp, lon_tmp, depth_tmp
                !print *, jj, lat_sta(jj), lon_sta(jj), z_sta(jj), dist_tmp

                if(dist_tmp .lt. dist_min) then
                  dist_min = dist_tmp
                  lon_min = lon_tmp
                  lat_min = lat_tmp
                  depth_min = depth_tmp
                  inc_angle_ini_min(jj) = inc_angle_ini
                  ray_azinc(2, i, j) = inc_angle_ini
                  ttime_min = ttime_tmp
                  width_min = width_tmp
                  !print '(a, 4(f8.4, 1x))', "rayshoot_tmp lon, lat, depth, dist_min = ", &
                  !&                                       lon_min, lat_min, depth_min, dist_min
                endif

                !!shooting the ray
                val_1d(1 : 2) = velocity(lon_index, lat_index, z_index : z_index + 1, wavetype)
                call linear_interpolation_1d(depth_tmp, zgrid, val_1d, velocity_interpolate(j))
                dvdz = (val_1d(2) - val_1d(1)) / dz_str
                call rayshooting3D(lon_tmp, lat_tmp, depth_tmp, az_tmp, inc_angle_tmp, time_step, velocity_interpolate(j), &
                &                  dvdlon, dvdlat, dvdz, lon_new, lat_new, depth_new, az_new, inc_angle_new)

                val_3d(1 : 2, 1 : 2, 1 : 2) &
                &  = qinv(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1, wavetype)
                call block_interpolation_3d(lon_tmp, lat_tmp, depth_tmp, xgrid, ygrid, zgrid, val_3d, qinv_interpolate(j))

                lon_tmp = lon_new
                lat_tmp = lat_new
                depth_tmp = depth_new
                az_tmp = az_new
                inc_angle_tmp = inc_angle_new
                ttime_tmp = ttime_tmp + time_step
                width_tmp = width_tmp + qinv_interpolate(j) * time_step
              enddo shooting_loop

            enddo incangle_loop
          enddo incangle_loop2
          write(0, '(a, 4(f8.4, 1x))') "masterevent lon, lat, depth, az_ini = ", &
          &                            evlon(j, mainloop_count), evlat(j, mainloop_count), evdp(j, mainloop_count), &
          &                            az_ini * rad2deg
          write(0, '(a, 3(f8.4, 1x))') "station lon, lat, depth = ", stlon(i), stlat(i), stdp(i)
          write(0, '(a, 4(f8.4, 1x))') "rayshoot lon, lat, depth, dist = ", lon_min, lat_min, depth_min, dist_min
          write(0, '(a, 2(f8.4, 1x))') "inc_angle (deg) = ", inc_angle_ini_min(nrayshoot) * rad2deg
          write(0, '(a, f5.2)') "traveltime (s) = ", ttime_min

#endif /* -DRAY_BENDING */
#endif /* -DV_CONST     */

          calamp(i, j) = evamp(j, mainloop_count) * exp(-pi * freq * width_min) / hypodist(i, j)

        enddo station_loop

        !!velocity and Qinv at each event location
        lon_index = int((evlon(j, mainloop_count) - lon_str_w) / dlon_str) + 1
        lat_index = int((evlat(j, mainloop_count) - lat_str_s) / dlat_str) + 1
        z_index   = int((evdp(j, mainloop_count)  - z_str_min) / dz_str  ) + 1
        xgrid(1) = lon_str_w + real(lon_index - 1, kind = fp) * dlon_str; xgrid(2) = xgrid(1) + dlon_str
        ygrid(1) = lat_str_s + real(lat_index - 1, kind = fp) * dlat_str; ygrid(2) = ygrid(1) + dlat_str
        zgrid(1) = z_str_min + real(z_index   - 1, kind = fp) * dz_str  ; zgrid(2) = zgrid(1) + dz_str
        val_3d(1 : 2, 1 : 2, 1 : 2) &
        &  = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1, wavetype)
        call linear_interpolation_3d(evlon(j, mainloop_count), evlat(j, mainloop_count), evdp(j, mainloop_count), &
        &                            xgrid, ygrid, zgrid, val_3d, velocity_interpolate(j))
        val_3d(1 : 2, 1 : 2, 1 : 2) &
        &  = qinv(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1, wavetype)
        call block_interpolation_3d(evlon(j, mainloop_count), evlat(j, mainloop_count), evdp(j, mainloop_count), &
        &                           xgrid, ygrid, zgrid, val_3d, qinv_interpolate(j))

      enddo event_loop
      !$omp end do
      !$omp end parallel

      event_count = 1
      do i = 1, nevent
        if(evflag(i, mainloop_count) .eqv. .true.) then
          event_index(i) = event_count
          event_index_rev(event_count) = i
          event_count = event_count + 1
        endif
      enddo

      !!Set up the observation vector and inversion matrix
      !!count number of amplitude ratio
      nobsamp_ratio = 0
      do j = 1, nevent - 1
        if(evflag(j, mainloop_count) .eqv. .false.) cycle
        do i = j + 1, nevent 
          if(evflag(i, mainloop_count) .eqv. .false.) cycle
          call greatcircle_dist(evlat(i, mainloop_count), evlon(i, mainloop_count), &
          &                     evlat(j, mainloop_count), evlon(j, mainloop_count), delta_out = epdelta)
          interevent_dist = sqrt((r_earth - evdp(i, mainloop_count)) ** 2 + (r_earth - evdp(j, mainloop_count)) ** 2 &
          &            - 2.0_fp * (r_earth - evdp(i, mainloop_count)) * (r_earth - evdp(j, mainloop_count)) * cos(epdelta))
          if(interevent_dist .le. interevent_dist_max) then
            do k = 1, nsta
              if(obsamp_flag(k, i) .eqv. .true. .and. obsamp_flag(k, j) .eqv. .true.) nobsamp_ratio = nobsamp_ratio + 1
            enddo
          endif
        enddo
      enddo

      !!set up matrix G
      nnonzero_elem = 8 * nobsamp_ratio + nconstraint * nevent_est
      allocate(irow(nnonzero_elem), icol(nnonzero_elem), nonzero_elem(nnonzero_elem), &
      &        obsvector(nobsamp_ratio + nconstraint), modelvector(nevent_est * 4))
      obsvector_count = 1
      do k = 1, nevent - 1
        if(evflag(k, mainloop_count) .eqv. .false.) cycle
        do j = k + 1, nevent
          if(evflag(j, mainloop_count) .eqv. .false.) cycle
          call greatcircle_dist(evlat(j, mainloop_count), evlon(j, mainloop_count), &
          &                     evlat(k, mainloop_count), evlon(k, mainloop_count), delta_out = epdelta)
          interevent_dist = sqrt((r_earth - evdp(j, mainloop_count)) ** 2 + (r_earth - evdp(k, mainloop_count)) ** 2 &
          &            - 2.0_fp * (r_earth - evdp(j, mainloop_count)) * (r_earth - evdp(k, mainloop_count)) * cos(epdelta))
          if(interevent_dist .le. interevent_dist_max) then
            do i = 1, nsta
              if(obsamp_flag(i, j) .eqv. .true. .and. obsamp_flag(i, k) .eqv. .true.) then
                obsvector(obsvector_count) = (log(obsamp(i, k)) - log(calamp(i, k))) &
                &                          - (log(obsamp(i, j)) - log(calamp(i, j))) &
                &                          +  log(evamp(k, mainloop_count)) - log(evamp(j, mainloop_count))

                normal_vector_j(1 : 3) = [sin(ray_azinc(2, i, j)) * cos(ray_azinc(1, i, j)), &
                &                         sin(ray_azinc(2, i, j)) * sin(ray_azinc(1, i, j)), &
                &                         cos(ray_azinc(2, i, j))]
                matrix_const_j = pi * freq * qinv_interpolate(j) / velocity_interpolate(j) + 1.0_fp / hypodist(i, j)
                normal_vector_k(1 : 3) = [sin(ray_azinc(2, i, k)) * cos(ray_azinc(1, i, k)), &
                &                         sin(ray_azinc(2, i, k)) * sin(ray_azinc(1, i, k)), &
                &                         cos(ray_azinc(2, i, k))]
                matrix_const_k = pi * freq * qinv_interpolate(k) / velocity_interpolate(k) + 1.0_fp / hypodist(i, k)

                irow(8 * (obsvector_count - 1) + 1 : 8 * (obsvector_count - 1) + 8) = obsvector_count
                icol(8 * (obsvector_count - 1) + 1)         = 4 * (event_index(k) - 1) + 1  !log(As_k)
                icol(8 * (obsvector_count - 1) + 2)         = 4 * (event_index(k) - 1) + 2  !Delta X_k
                icol(8 * (obsvector_count - 1) + 3)         = 4 * (event_index(k) - 1) + 3  !Delta Y_k
                icol(8 * (obsvector_count - 1) + 4)         = 4 * (event_index(k) - 1) + 4  !Delta Z_k
                icol(8 * (obsvector_count - 1) + 5)         = 4 * (event_index(j) - 1) + 1  !log(As_j)
                icol(8 * (obsvector_count - 1) + 6)         = 4 * (event_index(j) - 1) + 2  !Delta X_j
                icol(8 * (obsvector_count - 1) + 7)         = 4 * (event_index(j) - 1) + 3  !Delta Y_j
                icol(8 * (obsvector_count - 1) + 8)         = 4 * (event_index(j) - 1) + 4  !Delta Z_j
                nonzero_elem(8 * (obsvector_count - 1) + 1) =  1.0                                 !log(As_k(true))
                nonzero_elem(8 * (obsvector_count - 1) + 2) = -matrix_const_k * normal_vector_k(1) !Delta X_k
                nonzero_elem(8 * (obsvector_count - 1) + 3) = -matrix_const_k * normal_vector_k(2) !Delta Y_k
                nonzero_elem(8 * (obsvector_count - 1) + 4) = -matrix_const_k * normal_vector_k(3) !Delta Z_k
                nonzero_elem(8 * (obsvector_count - 1) + 5) = -1.0                                 !log(As_j(true))
                nonzero_elem(8 * (obsvector_count - 1) + 6) =  matrix_const_j * normal_vector_j(1) !Delta X_j
                nonzero_elem(8 * (obsvector_count - 1) + 7) =  matrix_const_j * normal_vector_j(2) !Delta Y_j
                nonzero_elem(8 * (obsvector_count - 1) + 8) =  matrix_const_j * normal_vector_j(3) !Delta Z_j
              
                obsvector_count = obsvector_count + 1
              endif
            enddo
          endif
        enddo
      enddo

      !!constraints
      do j = 1, nconstraint
        obsvector(obsvector_count + (j - 1)) = 0.0_fp
        do i = 1, nevent_est
          irow(8 * (obsvector_count - 1) + nevent_est * (j - 1) + i) = obsvector_count + (j - 1)
          icol(8 * (obsvector_count - 1) + nevent_est * (j - 1) + i) = 4 * (i - 1) + j
          nonzero_elem(8 * (obsvector_count - 1) + nevent_est * (j - 1) + i) = 1.0_fp
        enddo
      enddo

      !!call LSQR
      call solver%initialize(4 * nevent_est, nobsamp_ratio + nconstraint, nonzero_elem, irow, icol)
      call solver%solve(obsvector, zero, modelvector, istop)

      !!calculate residual
      residual_sum = 0.0_fp
      residual_tmp(1 : nevent_est * 4) = 0.0_fp
      do i = 1, 8 * nobsamp_ratio
        residual_tmp(irow(i)) = residual_tmp(irow(i)) + nonzero_elem(i) * modelvector(icol(i))
      enddo
      do i = 1, nobsamp_ratio
        residual_sum = residual_sum + (obsvector(i) - residual_tmp(i)) ** 2
      enddo
      residual_sum = residual_sum / real(nobsamp_ratio, kind = fp)
      if(residual_sum - residual_old .lt. 0.0_fp)             exit iter_loop
      if(residual_old - residual_sum .lt. delta_residual_max) exit iter_loop
      residual_old = residual_sum

      do i = 1, size(event_index_rev)
        delta_lat = modelvector(4 * (i - 1) + 2) / (r_earth - evdp(event_index_rev(i), mainloop_count)) * rad2deg
        delta_lon = modelvector(4 * (i - 1) + 3) / (r_earth - evdp(event_index_rev(i), mainloop_count) &
        &           * sin(0.5_fp * pi - evlat(event_index_rev(i), mainloop_count) * deg2rad))* rad2deg
        delta_depth = modelvector(4 * (i - 1) + 4)

        evamp(event_index_rev(i), mainloop_count) = exp(modelvector(4 * (i - 1) + 1))
        evlon(event_index_rev(i), mainloop_count) = evlon(event_index_rev(i), mainloop_count) - delta_lon
        call latgtoc(evlat(event_index_rev(i), mainloop_count), evlat_tmp)
        evlat_tmp = evlat_tmp - delta_lat
        call latctog(evlat_tmp, evlat(event_index_rev(i), mainloop_count)) 
        evdp(event_index_rev(i), mainloop_count) = evdp(event_index_rev(i), mainloop_count) - delta_depth
      enddo

    enddo iter_loop
  enddo main_loop


  !!output result
  sigma_amp = 0.0_fp
  sigma_lat = 0.0_fp
  sigma_lon = 0.0_fp
  sigma_depth = 0.0_fp
  open(unit = 10, file = trim(resultfile))
  write(10, '(a)') "# amp_ratio sigma_ampratio longitude sigma_lon latitude sigma_lat depth sigma_depth residual_sum nsta"
  do i = 1, nevent

    write(10, '(8(e15.8, 1x))') &
    &          evamp(i, 1), sigma_amp, evlon(i, 1), sigma_lon, evlat(i, 1), sigma_lat, evdp(i, 1), sigma_depth

    write(0, '(a, i0)')           "subevent index = ", i
    write(0, '(a, 2(e14.7, 1x))') "amp_ratio and sigma_amp = ", evamp(i, 1), sigma_amp
    write(0, '(a, 2(e14.7, 1x))') "longitude and sigma_lon = ", evlon(i, 1), sigma_lon
    write(0, '(a, 2(e14.7, 1x))') "latitude and sigma_lat = ",  evlat(i, 1), sigma_lat
    write(0, '(a, 2(e14.7, 1x))') "depth and sigma_depth = ",   evdp(i, 1), sigma_depth
  enddo
  close(10)

  stop
end program AmplitudeSourceLocation_DoubleDifference
