program AmplitudeSourceLocation_DoubleDifference
  !!Relative amplitude source location (Double-differnce of amplitude)
  !!using depth-dependent 1D velocity structure, 3D heterogeneous attenuation structure
  !!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2022
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
  use raybending,           only : pseudobending3D

  use lsqr_kinds
  use lsqr_module,          only : lsqr_solver_ez

#if defined (MPI)
  use mpi_module
#endif


  implicit none

  integer,            parameter :: wavetype = 2          !!1 for P-wave, 2 for S-wave
  real(kind = fp),    parameter :: snratio_accept = 0.1_fp
  real(kind = fp),    parameter :: delta_residual_max = 1.0e-6_fp
  real(kind = fp),    parameter :: constraint_weight_ini = 1.0_fp
  real(kind = fp),    parameter :: neventpair_ratio_max = 0.5_fp
  integer,            parameter :: nsta_use_min = 5
  integer,            parameter :: nconstraint = 4
  !integer,            parameter :: iter_loop_count_max = 500
  integer,            parameter :: iter_loop_count_max = 3
  !!Range for velocity and attenuation structure
  !!whole Japan
  real(kind = fp),    parameter :: lon_str_w = 122.0_fp, lon_str_e = 150.0_fp
  real(kind = fp),    parameter :: lat_str_s = 23.0_fp, lat_str_n = 46.0_fp
  real(kind = fp),    parameter :: z_str_min = -3.0_fp, z_str_max = 50.0_fp
  real(kind = fp),    parameter :: dlon_str = 0.1_fp, dlat_str = 0.1_fp, dz_str = 0.1_fp
  !!Meakandake volcano
  !real(kind = fp),    parameter :: lon_str_w = 143.5_fp, lon_str_e = 144.1_fp
  !real(kind = fp),    parameter :: lat_str_s = 43.0_fp, lat_str_n = 43.5_fp
  !real(kind = fp),    parameter :: z_str_min = -1.5_fp, z_str_max = 10.0_fp
  !real(kind = fp),    parameter :: dlon_str = 0.001_fp, dlat_str = 0.001_fp, dz_str = 0.1_fp

  integer,            parameter :: nlon_str = int((lon_str_e - lon_str_w) / dlon_str) + 2
  integer,            parameter :: nlat_str = int((lat_str_n - lat_str_s) / dlat_str) + 2
  integer,            parameter :: nz_str = int((z_str_max - z_str_min) / dz_str) + 2
  !!Ray shooting
  real(kind = fp),    parameter :: dvdlon = 0.0_fp, dvdlat = 0.0_fp         !!assume 1D structure
  integer,            parameter :: ninc_angle = 180                         !!grid search in incident angle
  integer,            parameter :: nrayshoot = 2                            !!number of grid search
  real(kind = fp),    parameter :: time_step = 0.005_fp
  real(kind = fp),    parameter :: rayshoot_dist_thr = 0.05_fp

  real(kind = fp),    parameter :: alt_to_depth = -1.0e-3_fp
  real(kind = fp),    parameter :: huge = 1.0e+6_fp
  real(kind = fp),    parameter :: epsilon = 1.0e-6_fp

  real(kind = fp),     allocatable :: stlon(:), stlat(:), stdp(:), ttime_cor(:, :)
  real(kind = dp),     allocatable :: topography(:, :), lon_topo(:), lat_topo(:)
  real(kind = fp),     allocatable :: obsamp(:, :), calamp(:, :), obsamp_noise(:), hypodist(:, :), ray_azinc(:, :, :), &
  &                                   evlon(:, :), evlat(:, :), evdp(:, :), evamp(:, :), interevent_dist(:, :), &
  &                                   interevent_dist_max(:, :), neventpair(:), siteamp(:), ray_baz(:, :), &
  &                                   qinv_interpolate(:), velocity_interpolate(:), residual_tmp(:), residual_min(:)
  integer,             allocatable :: event_index(:), event_index_rev(:)
  character(len = 6),  allocatable :: stname(:)
  character(len = 30), allocatable :: evid(:)
  logical,             allocatable :: stuse_flag(:), evflag(:, :), obsamp_flag(:, :), obsamp_flag_ini(:, :)
  character(len = 10)              :: freq_t

  real(kind = fp)               :: velocity(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                qinv(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                val_1d(1 : 2), val_2d(1 : 2, 1 : 2), val_3d(1 : 2, 1 : 2, 1 : 2), &
  &                                xgrid(1 : 2), ygrid(1 : 2), zgrid(1 : 2), inc_angle_ini_min(0 : nrayshoot), &
  &                                normal_vector_j(1 : 3), normal_vector_k(1 : 3)
  real(kind = fp)               :: damp, evlat_tmp, epdist, epdelta, az_ini, dinc_angle_org, dinc_angle, inc_angle_ini, &
  &                                lon_tmp, lat_tmp, depth_tmp, az_tmp, inc_angle_tmp, dist_tmp, ttime_tmp, width_tmp, &
  &                                dist_min, ttime_min, width_min, dvdz, lon_new, lat_new, depth_new, az_new, inc_angle_new, &
  &                                matrix_const_j, matrix_const_k, lon_min, lat_min, depth_min, delta_depth, delta_lon, &
  &                                delta_lat, residual_sum, residual_old, sigma_lon, sigma_lat, sigma_depth, sigma_amp, &
  &                                evlon_mean, evlat_mean, evdp_mean, evamp_mean, interevent_dist_max1, baz_diff, baz_diff_max, &
  &                                interevent_dist_max2, neventpair_maxval, var_siteamp_tmp, dist_weight, constraint_weight
  real(kind = dp)               :: topography_interpolate, dlon_topo, dlat_topo, freq
  integer                       :: nlon_topo, nlat_topo, nsta, nevent, lon_index, lat_index, z_index, &
  &                                i, j, k, ii, jj, kk, iter_loop_count, icount, nsta_use, ios, nampratio_each, &
  &                                obsvector_count, irow_count, icol_count, nnonzero_elem, nonzero_elem_count, &
  &                                nevent_est, nobsamp_ratio, mainloop_count, event_count, mainloop_count_max
  character(len = 129)          :: topo_grd, station_param, event_initloc_param, event_amp_param, resultfile
  character(len = 20)           :: cfmt, nsta_c, interevent_dist_max1_t, interevent_dist_max2_t, damp_t, baz_diff_max_t
  logical                       :: eventpair_flag

  !!Psuedo ray bending
  integer,                parameter :: ndiv_raypath = 10
  integer,                parameter :: nraypath_ini = 4
  real(kind = fp)                   :: raypath_lon((nraypath_ini - 1) * 2 ** ndiv_raypath + 1), &
  &                                    raypath_lat((nraypath_ini - 1) * 2 ** ndiv_raypath + 1), &
  &                                    raypath_dep((nraypath_ini - 1) * 2 ** ndiv_raypath + 1)
  integer                           :: nraypath

  !!LSQR variables
  type(lsqr_solver_ez)              :: solver
  real(kind = fp),      allocatable :: modelvector(:), obsvector(:), nonzero_elem(:)
  integer,              allocatable :: irow(:), icol(:)
  integer                           :: istop

#if defined (MPI)
  integer,              allocatable :: mainloop_count_begin(:), mainloop_count_end(:)
#endif

#if defined (MPI)
  call mpi_ini_rank_size
  call mpi_get_hostname
#endif

#if defined (MPI)
  if(mpi_rank .eq. 0) then
#endif
    icount = command_argument_count()

    if(icount .ne. 10) then
      write(0, '(a)', advance="no") "usage: ./asl_dd "
      write(0, '(a)', advance="no") "(topography_grd) (station_param_file) (event_initial_location_file) (event_amplitude_file) "
      write(0, '(a)', advance="no") "(frequency) (maximum interevent distance1) (maximum_interevent_distance2) "
      write(0, '(a)')               "(baz_diff_max) (damping factor) (result_file)"
#if defined (MPI)
      call mpiabort
#endif
      error stop
    endif
    call getarg(1, topo_grd)
    call getarg(2, station_param)
    call getarg(3, event_initloc_param)
    call getarg(4, event_amp_param)
    call getarg(5, freq_t); read(freq_t, *) freq
    call getarg(6, interevent_dist_max1_t); read(interevent_dist_max1_t, *) interevent_dist_max1
    call getarg(7, interevent_dist_max2_t); read(interevent_dist_max2_t, *) interevent_dist_max2
    call getarg(8, baz_diff_max_t); read(baz_diff_max_t, *) baz_diff_max
    call getarg(9, damp_t); read(damp_t, *) damp
    call getarg(10, resultfile)

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
    &        stuse_flag(1 : nsta), obsamp_noise(1 : nsta), siteamp(1 : nsta))
    do i = 1, nsta
      read(10, *) stlon(i), stlat(i), stdp(i), stname(i), stuse_flag(i), &
      &           ttime_cor(i, 1), ttime_cor(i, 2), siteamp(i), var_siteamp_tmp, obsamp_noise(i)
      write(0, '(a, i0, a, f9.4, a, f8.4, a, f6.3, 1x, a7, l2)') &
      &     "station(", i, ") lon(deg) = ", stlon(i), " lat(deg) = ", stlat(i), " depth(km) = ", stdp(i), &
      &     trim(stname(i)), stuse_flag(i)
#if defined (TESTDATA)
      obsamp_noise(i) = 1.0e-6_fp
#endif
    enddo
    close(10)
    if(nsta .lt. nsta_use_min) then
      close(10)
      write(0, '(a, i0)') "Number of stations for calculation should be larger than ", nsta_use_min
#if defined (MPI)
      call mpiabort
#endif
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
    allocate(evlon(1 : nevent, 0 : nsta + 1), evlat(1 : nevent, 0 : nsta + 1), evdp(1 : nevent, 0 : nsta + 1), &
    &        evamp(1 : nevent, 0 : nsta + 1), evflag(1 : nevent, 0 : nsta + 1), evid(1 : nevent))
    allocate(hypodist(1 : nsta, 1 : nevent), ray_azinc(1 : 2, 1 : nsta, 1 : nevent), obsamp_flag(1 : nsta, 1 : nevent), &
    &        obsamp_flag_ini(1 : nsta, 1 : nevent), obsamp(1 : nsta, 1 : nevent), calamp(1 : nsta, 1 : nevent), &
    &        ray_baz(1 : nsta, 1 : nevent))
    allocate(velocity_interpolate(1 : nevent), qinv_interpolate(1 : nevent), event_index(1 : nevent), &
    &        event_index_rev(1 : nevent), neventpair(1 : nevent))
    read(10, *)
    do i = 1, nevent
      !!output from AmplitudeSourceLocation_masterevent.F90
      !read(10, *) evamp(i, 0), sigma_amp, evlon(i, 0), sigma_lon, evlat(i, 0), sigma_lat, evdp(i, 0), sigma_depth, evid(i)
      !!output from AmplitudeSourceLocation_PulseWidth.F90
      read(10, *) evid(i), evlon(i, 0), evlat(i, 0), evdp(i, 0), evamp(i, 0)

      do j = 1, nsta + 1
        evlon (i, j) = evlon(i, 0)
        evlat (i, j) = evlat(i, 0)
        evdp  (i, j) = evdp (i, 0)
        evamp (i, j) = evamp(i, 0)
        evflag(i, j) = .true.
      enddo
    enddo
    close(10)
    
    open(unit = 10, file = event_amp_param)
    read(10, *)
    do j = 1, nevent
      obsamp_flag_ini(1 : nsta, j) = .false.
      read(10, *) (obsamp(i, j), i = 1, nsta)
      nsta_use = 0
      do i = 1, nsta
        if(stuse_flag(i) .eqv. .true.) then
          if(obsamp(i, j) / obsamp_noise(i) .ge. snratio_accept) then
            nsta_use = nsta_use + 1
            obsamp_flag_ini(i, j) = .true.
          endif
        endif
      enddo
      if(nsta_use .lt. nsta_use_min) then
        obsamp_flag_ini(1 : nsta, j)  = .false.
        evflag(j, 0 : nsta + 1) = .false.
      endif
      !print *, j, evflag(j, 0), nsta_use
    enddo
    close(10)

    !!set velocity/attenuation structure
    call set_velocity(z_str_min, dz_str, velocity, qinv)

    write(0, '(2a)') "Damping factor = ", trim(damp_t)

#if defined (MPI)
  endif
#endif
#if defined (MPI)
  !!Topography
  call mpi_bcast_int(nlon_topo, 0)
  call mpi_bcast_int(nlat_topo, 0)
  call mpi_bcast_dp(dlon_topo, 0)
  call mpi_bcast_dp(dlat_topo, 0)
  if(mpi_rank .ne. 0) allocate(topography(1 : nlon_topo, 1 : nlat_topo), lon_topo(1 : nlon_topo), lat_topo(1 : nlat_topo))
  call mpi_bcast_dp_1d(lon_topo, 0)
  call mpi_bcast_dp_1d(lat_topo, 0)
  call mpi_bcast_dp_2d(topography, 0)
  !!station
  call mpi_bcast_int(nsta, 0)
  if(mpi_rank .ne. 0) allocate(stlon(1 : nsta), stlat(1 : nsta), stdp(1 : nsta), siteamp(1 : nsta), stuse_flag(1 : nsta), &
  &                            obsamp_noise(1 : nsta))
  call mpi_bcast_fp_1d(stlon, 0)
  call mpi_bcast_fp_1d(stlat, 0)
  call mpi_bcast_fp_1d(stdp, 0)
  call mpi_bcast_logical_1d(stuse_flag, 0)
  call mpi_bcast_fp_1d(siteamp, 0)
  call mpi_bcast_fp_1d(obsamp_noise, 0)
  !!event
  call mpi_bcast_int(nevent, 0)
  if(mpi_rank .ne. 0) then
    allocate(evlon(1 : nevent, 0 : nsta + 1), evlat(1 : nevent, 0 : nsta + 1), evdp(1 : nevent, 0 : nsta + 1), &
    &        evamp(1 : nevent, 0 : nsta + 1), evflag(1 : nevent, 0 : nsta + 1), ray_baz(1 : nsta, 1 : nevent), &
    &        hypodist(1 : nsta, 1 : nevent), ray_azinc(1 : 2, 1 : nsta, 1 : nevent), obsamp_flag_ini(1 : nsta, 1 : nevent), &
    &        obsamp_flag(1 : nsta, 1 : nevent), obsamp(1 : nsta, 1 : nevent), calamp(1 : nsta, 1 : nevent), &
    &        velocity_interpolate(1 : nevent), qinv_interpolate(1 : nevent), event_index(1 : nevent), &
    &        event_index_rev(1 : nevent), neventpair(1 : nevent))
  endif
  call mpi_bcast_fp_2d(evlon, 0)
  call mpi_bcast_fp_2d(evlat, 0)
  call mpi_bcast_fp_2d(evdp, 0)
  call mpi_bcast_fp_2d(evamp, 0)
  call mpi_bcast_logical_2d(evflag, 0)
  call mpi_bcast_fp_2d(obsamp, 0)
  call mpi_bcast_logical_2d(obsamp_flag_ini, 0)
  !!velocity and attenuation structure
  call mpi_bcast_fp_4d(velocity, 0)
  call mpi_bcast_fp_4d(qinv, 0)
  !!calculation parameter
  call mpi_bcast_fp(freq, 0)
  call mpi_bcast_fp(damp, 0)
  call mpi_bcast_fp(baz_diff_max, 0)
  call mpi_bcast_fp(interevent_dist_max1, 0)
  call mpi_bcast_fp(interevent_dist_max2, 0)
#endif

  

#if defined (WITHOUT_ERROR)
  mainloop_count_max = 1
#else
  mainloop_count_max = nsta + 1
#endif

  allocate(residual_min(1 : mainloop_count_max))

#if defined (MPI)
  allocate(mainloop_count_begin(0 : mpi_size - 1), mainloop_count_end(0 : mpi_size - 1))

#if defined (WITHOUT_ERROR)
  mainloop_count_begin(0 : mpi_size - 1) = 1
  mainloop_count_end(0) = 1
  mainloop_count_end(1 : mpi_size - 1) = 0
#else
  do i = 0, mpi_size - 1
    mainloop_count_begin(i) = int(mainloop_count_max / mpi_size) * i + 1
  enddo
  do i = 0, mpi_size - 1
    if(i .ne. mpi_size - 1) then
      mainloop_count_end(i) = mainloop_count_begin(i + 1) - 1
    else
      mainloop_count_end(i) = mainloop_count_max
    endif
  enddo

  write(0, '(a, i0)', advance = "no") "mpi_rank = ", mpi_rank
  write(0, '(a, 2(i0, 1x))') " Mainloop_count_range = ", mainloop_count_begin(mpi_rank), mainloop_count_end(mpi_rank)
#endif

#endif

  !!Main loop: estimate \Delta_x, \Delta_y, \Delta_z, \Delta_amp until convergence
#if defined (MPI)
  main_loop: do mainloop_count = mainloop_count_begin(mpi_rank), mainloop_count_end(mpi_rank)
#else
  !main_loop: do mainloop_count = 1, mainloop_count_max
  main_loop: do mainloop_count = 1, mainloop_count_max
#endif

    obsamp_flag(1 : nsta, 1 : nevent) = obsamp_flag_ini(1 : nsta, 1 : nevent)
    if(mainloop_count .gt. 1) then
      do i = 1, nevent
        obsamp_flag(mainloop_count - 1, i) = .false.
      enddo
    endif
   

    residual_old = huge
    residual_min(mainloop_count) = huge
    iter_loop_count = 0
    iter_loop: do
      iter_loop_count = iter_loop_count + 1   
      if(iter_loop_count .gt. iter_loop_count_max) exit iter_loop
#if defined (MPI)
      write(0, '(a, i0, a)', advance = "no") "mpi_rank = ", mpi_rank, " "
#endif
      write(0, '(2(a, i0))') "Mainloop count = ", mainloop_count, " Iteration loop count = ", iter_loop_count
      !!calculate ray length, pulse width, unit vector of ray incident
      write(0, '(a)') "calculate ray length, pulse width, and ray incident vector for each event"

      nevent_est = 0

      !!calculate traveltime from each event location to stations
      event_loop: do j = 1, nevent
        !!check whether the depth of each event location is lower than the topo
        lon_index = int((evlon(j, mainloop_count) - lon_topo(1)) / dlon_topo) + 1
        lat_index = int((evlat(j, mainloop_count) - lat_topo(1)) / dlat_topo) + 1
        xgrid(1 : 2) = [lon_topo(lon_index), lon_topo(lon_index + 1)]
        ygrid(1 : 2) = [lat_topo(lat_index), lat_topo(lat_index + 1)]
        val_2d(1 : 2, 1 : 2) = topography(lon_index : lon_index + 1, lat_index : lat_index + 1)
        call linear_interpolation_2d(evlon(j, mainloop_count), evlat(j, mainloop_count), xgrid, ygrid, val_2d, &
        &                            topography_interpolate)

        if(evdp(j, mainloop_count) .lt. topography_interpolate) evflag(j, mainloop_count) = .false.
        if(evflag(j, mainloop_count) .eqv. .false.) cycle

        nevent_est = nevent_est + 1
  
        station_loop: do i = 1, nsta
          if(obsamp_flag(i, j) .eqv. .false.) cycle
          !!calculate azimuth and hypocentral distance
          call greatcircle_dist(evlat(j, mainloop_count), evlon(j, mainloop_count), stlat(i), stlon(i), &
          &                     distance = epdist, azimuth = az_ini, backazimuth = ray_baz(i, j), delta_out = epdelta)
          lon_index = int((stlon(i) - lon_str_w) / dlon_str) + 1
          lat_index = int((stlat(i) - lat_str_s) / dlat_str) + 1
          z_index   = int((stdp(i) - z_str_min) / dz_str) + 1
          hypodist(i, j) = sqrt((r_earth - evdp(j, mainloop_count)) ** 2 + (r_earth - stdp(i)) ** 2 &
          &                    - 2.0_fp * (r_earth - evdp(j, mainloop_count)) * (r_earth - stdp(i)) * cos(epdelta))
          ray_azinc(1, i, j) = az_ini

#if defined (V_CONST)
          !!homogeneous structure: ray incident angle is calculated using cosine function (assuming cartesian coordinate)
          ray_azinc(2, i, j) = acos(epdist / hypodist(i, j)) + pi * 0.5_fp
          ttime_min = hypodist(i, j) / velocity(lon_index, lat_index, z_index, wavetype)
          width_min = ttime_min * qinv(lon_index, lat_index, z_index, wavetype)
#else
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

#endif /* -DV_CONST     */

          calamp(i, j) = evamp(j, mainloop_count) * exp(-pi * freq * width_min) / hypodist(i, j) * siteamp(i)

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

#if defined (MPI)
      write(0, '(a, i0, a)', advance = "no") "mpi_rank = ", mpi_rank, " "
#endif
      write(0, '(a, i0)') "nevent_est = ", nevent_est

      event_count = 1
      do i = 1, nevent
        if(evflag(i, mainloop_count) .eqv. .true.) then
          event_index(i) = event_count
          event_index_rev(event_count) = i
          event_count = event_count + 1
        endif
      enddo

      !!Set up the observation vector and inversion matrix
      !!count number of amplitude ratio, event pairs
      nobsamp_ratio = 0
      neventpair(1 : nevent) = 0.0_fp
      allocate(interevent_dist(1 : nevent, 1 : nevent), interevent_dist_max(1 : nevent, 1 : nevent))
      interevent_dist(1 : nevent, 1 : nevent) = huge
      interevent_dist_max(1 : nevent, 1 : nevent) = interevent_dist_max2
      do j = 1, nevent - 1
        if(evflag(j, mainloop_count) .eqv. .false.) cycle
        do i = j + 1, nevent 
          if(evflag(i, mainloop_count) .eqv. .false.) cycle
          call greatcircle_dist(evlat(i, mainloop_count), evlon(i, mainloop_count), &
          &                     evlat(j, mainloop_count), evlon(j, mainloop_count), delta_out = epdelta)
          interevent_dist(i, j) = sqrt((r_earth - evdp(i, mainloop_count)) ** 2 + (r_earth - evdp(j, mainloop_count)) ** 2 &
          &           - 2.0_fp * (r_earth - evdp(i, mainloop_count)) * (r_earth - evdp(j, mainloop_count)) * cos(epdelta))

          if(interevent_dist(i, j) .le. interevent_dist_max(i, j)) then
            eventpair_flag= .false.
            do k = 1, nsta
              if(obsamp_flag(k, i) .eqv. .false.) cycle
              if(obsamp_flag(k, j) .eqv. .false.) cycle
              baz_diff = abs(ray_baz(k, j) - ray_baz(k, i))
              if(baz_diff .gt. pi) baz_diff = 2.0_fp * pi - baz_diff
              if(baz_diff .gt. baz_diff_max * deg2rad) cycle
              nobsamp_ratio = nobsamp_ratio + 1
              eventpair_flag = .true.
            enddo
            if(eventpair_flag .eqv. .true.) then
              neventpair(i) = neventpair(i) + 1.0_fp
              neventpair(j) = neventpair(j) + 1.0_fp
            endif
          endif
        enddo
      enddo
#if defined (MPI)
      write(0, '(a, i0, a)', advance = "no") "mpi_rank = ", mpi_rank, " "
#endif
      write(0, '(a, i0)') "Number of obsamp_ratio (within interevent_dist_max2) = ", nobsamp_ratio
      neventpair_maxval = maxval(neventpair)

      !!recount amplitude ratio, event pairs
      nobsamp_ratio = 0
      do j = 1, nevent - 1
        if(evflag(j, mainloop_count) .eqv. .false.) cycle
        do i = j + 1, nevent 
          if(evflag(i, mainloop_count) .eqv. .false.) cycle

          if(.not. (neventpair(i) / neventpair_maxval .lt. neventpair_ratio_max .and. &
          &         neventpair(j) / neventpair_maxval .lt. neventpair_ratio_max)) &
          &  interevent_dist_max(i, j) = interevent_dist_max1

          if(interevent_dist(i, j) .le. interevent_dist_max(i, j)) then
            do k = 1, nsta
              if(obsamp_flag(k, i) .eqv. .false.) cycle
              if(obsamp_flag(k, j) .eqv. .false.) cycle
              baz_diff = abs(ray_baz(k, j) - ray_baz(k, i))
              if(baz_diff .gt. pi) baz_diff = 2.0_fp * pi - baz_diff
              if(baz_diff .gt. baz_diff_max * deg2rad) cycle
              nobsamp_ratio = nobsamp_ratio + 1
            enddo
          endif
        enddo
      enddo

#if defined (MPI)
      write(0, '(a, i0, a)', advance = "no") "mpi_rank = ", mpi_rank, " "
#endif
      write(0, '(a, i0)') "number of obsamp_ratio = ", nobsamp_ratio

      !!set up matrix G
      nnonzero_elem = 8 * nobsamp_ratio + nconstraint * nevent_est
      allocate(irow(nnonzero_elem), icol(nnonzero_elem), nonzero_elem(nnonzero_elem), &
      &        obsvector(nobsamp_ratio + nconstraint), modelvector(nevent_est * 4))
      obsvector_count = 0
      obsvector(1 : nobsamp_ratio + nconstraint) = 0.0_fp
      modelvector(1 : nevent_est * 4) = 0
      do k = 1, nevent - 1
        do j = k + 1, nevent

          if(interevent_dist(j, k) .le. interevent_dist_max(j, k)) then
            dist_weight = exp(-(interevent_dist(j, k) ** 2) / (interevent_dist_max2 ** 2))

            do i = 1, nsta
              if(obsamp_flag(i, j) .eqv. .false.) cycle
              if(obsamp_flag(i, k) .eqv. .false.) cycle
              baz_diff = abs(ray_baz(i, k) - ray_baz(i, j))
              if(baz_diff .gt. pi) baz_diff = 2.0_fp * pi - baz_diff
              if(baz_diff .gt. baz_diff_max * deg2rad) cycle

              obsvector_count = obsvector_count + 1

              obsvector(obsvector_count) = (log(obsamp(i, k) / calamp(i, k)) - log(obsamp(i, j) / calamp(i, j))) * dist_weight

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

              nonzero_elem(8 * (obsvector_count - 1) + 1 : 8 * (obsvector_count - 1) + 8) &
              &  = nonzero_elem(8 * (obsvector_count - 1) + 1 : 8 * (obsvector_count - 1) + 8) * dist_weight

            enddo
          endif
        enddo
      enddo
      deallocate(interevent_dist, interevent_dist_max)

      !!constraints
      !constraint_weight = constraint_weight_ini / real(iter_loop_count, kind = fp)
      constraint_weight = &
      &  constraint_weight_ini * (1.0_fp - real(iter_loop_count, kind = fp) / real(iter_loop_count_max, kind = fp))
      do j = 1, nconstraint
        obsvector(obsvector_count + j) = 0.0_fp
        do i = 1, nevent_est
          irow(8 * obsvector_count + nevent_est * (j - 1) + i) = obsvector_count + j
          icol(8 * obsvector_count + nevent_est * (j - 1) + i) = 4 * (i - 1) + j
          nonzero_elem(8 * obsvector_count + nevent_est * (j - 1) + i) = 1.0_fp * constraint_weight
        enddo
      enddo

      if(mpi_rank .eq. 7) then
        open(unit = 10, file = "log")
        do i = 1, size(obsvector)
          write(10, *) obsvector(i)
        enddo
        close(10)
      endif

      !!call LSQR
      call solver%initialize(nobsamp_ratio + nconstraint, 4 * nevent_est, nonzero_elem, irow, icol)
      call solver%solve(obsvector, damp, modelvector, istop)
      write(0, '(a, i0)') "LSQR istop = ", istop
      if(istop .eq. 5) then
        write(0, '(a)') "LSQR did not converge, damping factor should be changed"
#ifdef MPI
        call mpiabort
#endif
        error stop
      endif

      !!calculate residual
      allocate(residual_tmp(1 : nobsamp_ratio))
      residual_sum = 0.0_fp
      residual_tmp(1 : nobsamp_ratio) = 0.0_fp
      do i = 1, 8 * nobsamp_ratio
        residual_tmp(irow(i)) = residual_tmp(irow(i)) + nonzero_elem(i) * modelvector(icol(i))
      enddo
      do i = 1, nobsamp_ratio
        residual_sum = residual_sum + (obsvector(i) - residual_tmp(i)) ** 2
      enddo
      residual_sum = residual_sum / real(nobsamp_ratio, kind = fp)
      deallocate(residual_tmp)
#if defined (MPI)
      write(0, '(a, i0, a)', advance = "no") "mpi_rank = ", mpi_rank, " "
#endif
      write(0, '(a, 2(e15.7, 1x))') "residual_old and residual_sum = ", residual_old, residual_sum
      residual_min(mainloop_count) = min(residual_min(mainloop_count), residual_sum)
      if(residual_old - residual_sum .lt. delta_residual_max) then
        deallocate(irow, icol, nonzero_elem, obsvector, modelvector)
        exit iter_loop
      endif
      residual_old = residual_sum

      do i = 1, nevent_est
        delta_lat = modelvector(4 * (i - 1) + 2) / (r_earth - evdp(event_index_rev(i), mainloop_count)) * rad2deg
        delta_lon = modelvector(4 * (i - 1) + 3) / (r_earth - evdp(event_index_rev(i), mainloop_count) &
        &           * sin(0.5_fp * pi - evlat(event_index_rev(i), mainloop_count) * deg2rad))* rad2deg
        delta_depth = modelvector(4 * (i - 1) + 4)

        evamp(event_index_rev(i), mainloop_count) = exp(modelvector(4 * (i - 1) + 1)) * evamp(event_index_rev(i), mainloop_count)
        evlon(event_index_rev(i), mainloop_count) = evlon(event_index_rev(i), mainloop_count) - delta_lon
        call latgtoc(evlat(event_index_rev(i), mainloop_count), evlat_tmp)
        evlat_tmp = (evlat_tmp * rad2deg - delta_lat) * deg2rad
        call latctog(evlat_tmp, evlat(event_index_rev(i), mainloop_count)) 
        evdp(event_index_rev(i), mainloop_count) = evdp(event_index_rev(i), mainloop_count) - delta_depth
      enddo
      deallocate(irow, icol, nonzero_elem, obsvector, modelvector)


    enddo iter_loop
  enddo main_loop

#if defined (WITHOUT_ERROR)
#else
#if defined (MPI)
  !!send results to other processes
  do i = 0, mpi_size - 1
    call mpi_bcast_fp_2d(evlon(:, mainloop_count_begin(i) : mainloop_count_end(i)), i)
    call mpi_bcast_fp_2d(evlat(:, mainloop_count_begin(i) : mainloop_count_end(i)), i)
    call mpi_bcast_fp_2d(evdp(:, mainloop_count_begin(i) : mainloop_count_end(i)), i)
    call mpi_bcast_fp_2d(evamp(:, mainloop_count_begin(i) : mainloop_count_end(i)), i)
    call mpi_bcast_logical_2d(evflag(:, mainloop_count_begin(i) : mainloop_count_end(i)), i)
    call mpi_bcast_fp_1d(residual_min(mainloop_count_begin(i) : mainloop_count_end(i)), i)
  enddo
#endif 
#endif

  !!output result
#if defined (MPI)
  if(mpi_rank .eq. 0) then
#endif

    open(unit = 10, file = trim(resultfile))
    write(10, '(4a)', advance = "no") "# intereven_dist_max1 = ", trim(interevent_dist_max1_t)
    write(10, '(4a)', advance = "no") " intereven_dist_max2 = ", trim(interevent_dist_max2_t), " damping factor = ", trim(damp_t)
    write(10, '(a, f5.2)', advance = "no") " freq (Hz) = ", freq
    write(10, '(a, e15.7)') " Residual_min = ", residual_min(1)
    write(10, '(a)') "# amp sigma_amp(in log scale) longitude sigma_lon latitude sigma_lat depth sigma_depth evflag evid"

    do j = 1, nevent
      sigma_amp = 0.0_fp
      sigma_lat = 0.0_fp
      sigma_lon = 0.0_fp
      sigma_depth = 0.0_fp

#if defined (WITHOUT_ERROR)
#else
      !!calculate estimation errors
      icount = 0
      evlon_mean = 0.0_fp
      evlat_mean = 0.0_fp
      evdp_mean = 0.0_fp
      evamp_mean = 0.0_fp
      do i = 2, nsta + 1
        if(evflag(j, i) .eqv. .false.) cycle
        if(.not. (abs(evlon(j, 1) - evlon(j, i)) .le. epsilon .and. &
        &         abs(evlat(j, 1) - evlat(j, i)) .le. epsilon .and. &
        &         abs(evdp(j, 1)  - evdp(j, i)) .le. epsilon .and. &
        &         abs(log(evamp(j, 1)) - log(evamp(j, i))) .le. epsilon)) then
          icount = icount + 1
          evlon_mean = evlon_mean + evlon(j, i)
          evlat_mean = evlat_mean + evlat(j, i)
          evdp_mean  = evdp_mean  + evdp(j, i)
          evamp_mean = evamp_mean + log(evamp(j, i))
        endif
      enddo
      if(icount .ge. 1) then
        evlon_mean = evlon_mean / real(icount, kind = fp)
        evlat_mean = evlat_mean / real(icount, kind = fp)
        evdp_mean  = evdp_mean  / real(icount, kind = fp)
        evamp_mean = evamp_mean / real(icount, kind = fp)
        do i = 2, nsta + 1
          if(evflag(j, i) .eqv. .false.) cycle
          if(.not. (abs(evlon(j, 1) - evlon(j, i)) .le. epsilon .and. &
          &         abs(evlat(j, 1) - evlat(j, i)) .le. epsilon .and. &
          &         abs(evdp(j, 1)  - evdp(j, i)) .le. epsilon .and. &
          &         abs(log(evamp(j, 1)) - log(evamp(j, i))) .le. epsilon)) then
            sigma_lon   = sigma_lon   + (evlon_mean - evlon(j, i)) ** 2
            sigma_lat   = sigma_lat   + (evlat_mean - evlat(j, i)) ** 2
            sigma_depth = sigma_depth + (evdp_mean  - evdp(j, i)) ** 2
            sigma_amp   = sigma_amp   + (evamp_mean - log(evamp(j, i))) ** 2
          endif
        enddo
        sigma_lon   = sqrt(real(icount - 1, kind = fp) / real(icount, kind = fp) * sigma_lon)
        sigma_lat   = sqrt(real(icount - 1, kind = fp) / real(icount, kind = fp) * sigma_lat)
        sigma_depth = sqrt(real(icount - 1, kind = fp) / real(icount, kind = fp) * sigma_depth)
        sigma_amp   = sqrt(real(icount - 1, kind = fp) / real(icount, kind = fp) * sigma_amp)
      endif

#endif

      write(10, '(8(e15.8, 1x), l1, 1x, a)') &
      &          evamp(j, 1), sigma_amp, evlon(j, 1), sigma_lon, evlat(j, 1), sigma_lat, evdp(j, 1), sigma_depth, &
      &          evflag(j, 1), trim(evid(j))

      !write(0, '(a, i0, 1x, a, 1x, l1)') "subevent index = ", j, trim(evid(j)), evflag(j, 1)
      !write(0, '(a, 2(e14.7, 1x))')      "amp_ratio and sigma_amp = ", evamp(j, 1), sigma_amp
      !write(0, '(a, 2(e14.7, 1x))')      "longitude and sigma_lon = ", evlon(j, 1), sigma_lon
      !write(0, '(a, 2(e14.7, 1x))')      "latitude and sigma_lat = ",  evlat(j, 1), sigma_lat
      !write(0, '(a, 2(e14.7, 1x))')      "depth and sigma_depth = ",   evdp(j, 1), sigma_depth
    enddo
    close(10)

#if defined (MPI)
  endif
  call mpifinalize
#endif

  stop
end program AmplitudeSourceLocation_DoubleDifference

