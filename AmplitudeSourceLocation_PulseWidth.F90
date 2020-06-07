program AmplitudeSourceLocation_PulseWidth
  !!Amplitude Source Location using depth-dependent 1D velocity structure, 3D heterogeneous attenuation structure

  use nrtype,               only : fp, sp, dp
  use constants,            only : rad2deg, deg2rad, pi, r_earth
  use rayshooting,          only : rayshooting3D
  use read_sacfile,         only : read_sachdr, read_sacdata
  use set_velocity_model,   only : set_velocity
  use linear_interpolation, only : linear_interpolation_1d, linear_interpolation_2d, block_interpolation_3d
  use greatcircle,          only : greatcircle_dist
  use itoa,                 only : int_to_char
  use GMT,                  only : grd_create
#ifdef WIN
  use m_win
  use m_winch
#endif
  !$ use omp_lib

  implicit none

  real(kind = fp),    parameter :: dlon = 0.001_fp, dlat = 0.001_fp, dz = 0.1_fp
  real(kind = fp),    parameter :: lon_w = 143.90_fp, lon_e = 144.05_fp
  real(kind = fp),    parameter :: lat_s = 43.35_fp, lat_n = 43.410_fp
  real(kind = fp),    parameter :: z_min = -1.5_fp, z_max = 3.2_fp
  real(kind = fp),    parameter :: dvdlon = 0.0_fp, dvdlat = 0.0_fp
  integer,            parameter :: ninc_angle = 200                         !!ray shooting
  integer,            parameter :: nsta = 4
  character(len = 6), parameter :: stname(1 : nsta) = ["V.MEAB", "V.MEAA", "V.PMNS", "V.NSYM"]
#ifdef WIN
  character(len = 4), parameter :: st_winch(1 : nsta) = ["2724", "13F1", "274D", "2720"]
#endif
  real(kind = dp),    parameter :: siteamp(1 : nsta) = [1.0_dp, 0.738_dp, 2.213_dp, 1.487_dp]

  real(kind = fp),    parameter :: time_step = 0.01_fp
  integer,            parameter :: maxlen = 4
  real(kind = dp),    parameter :: fl = 5.0_dp, fh = 10.0_dp, fs = 12.0_dp  !!bandpass filter parameters
  real(kind = dp),    parameter :: ap = 0.5_dp, as = 5.0_dp                 !!bandpass filter parameters
  real(kind = dp),    parameter :: order = 1.0e-3_dp                        !!nano m/s to micro m/s
  real(kind = fp),    parameter :: alt_to_depth = -1.0e-3_fp
  real(kind = dp),    parameter :: huge = 1.0e+5_dp

  real(kind = fp),    parameter :: dinc_angle = pi / real(ninc_angle, kind = fp)
  integer,            parameter :: nlon = int((lon_e - lon_w) / dlon) + 1
  integer,            parameter :: nlat = int((lat_n - lat_s) / dlat) + 2
  integer,            parameter :: nz   = int((z_max - z_min) / dz) + 2
  real(kind = dp),    parameter :: freq = (fl + fh) * 0.5_dp

  real(kind = fp)               :: velocity(1 : nlon, 1 : nlat, 1 : nz), qinv(1 : nlon, 1 : nlat, 1 : nz), &
  &                                topography(1 : nlon, 1 : nlat), sampling(1 : nsta), begin(1 : nsta), &
  &                                val_1d(1 : 2), val_3d(1 : 2, 1 : 2, 1 : 2), xgrid(1 : 2), ygrid(1 : 2), zgrid(1 : 2), &
  &                                lon_sta(1 : nsta), lat_sta(1 : nsta), z_sta(1 : nsta), &
  &                                width_min(1 : nsta, 1 : nlon, 1 : nlat, 1 : nz), &
  &                                ttime_min(1 : nsta, 1 : nlon, 1 : nlat, 1 : nz), &
  &                                hypodist(1 : nsta, 1 : nlon, 1 : nlat, 1 : nz)
  real(kind = dp)               :: rms_amp_obs(1 : nsta), residual(1 : nlon, 1 : nlat, 1 : nz), &
  &                                xrange(1 : 2), yrange(1 : 2), spacing(1 : 2), source_amp(1 : nlon, 1 : nlat, 1 : nz)
  integer                       :: npts(1 : nsta), residual_minloc(3)
  real(kind = sp), allocatable  :: residual_grd(:, :)
  real(kind = dp), allocatable  :: waveform_obs(:, :)
  
  real(kind = fp)               :: ttime_tmp, width_tmp, velocity_interpolate, qinv_interpolate, &
  &                                az_tmp, inc_angle_tmp, az_new, inc_angle_new, az_ini, &
  &                                lon_tmp, lat_tmp, depth_tmp, lon_new, lat_new, depth_new, origintime, &
  &                                lon_grid, lat_grid, depth_grid, dist_min, dist_tmp, dvdz, epdelta, &
  &                                ot_begin, ot_end, ot_shift, rms_tw
  real(kind = sp)               :: lon_r, lat_r, topo_r
  real(kind = dp)               :: residual_normalize, amp_avg
  
  integer                       :: i, j, k, ii, jj, icount, wave_index, time_count, lon_index, lat_index, z_index, &
  &                                grd_status, npts_max
  character(len = 129)          :: dem_file, sacfile, sacfile_index, ot_begin_t, ot_end_t, rms_tw_t, &
  &                                grdfile, resultfile, resultdir
  character(len = maxlen)       :: time_count_char

  !!filter variables
  real(kind = dp), allocatable  :: h(:), waveform_tmp(:)
  real(kind = dp)               :: gn, c
  integer                       :: m, n

#ifdef WIN
  !!in case of win file input
  character(len = 129)          :: win_filename, win_chfilename
  real(dp), parameter           :: order_um = 1.0e+6_dp                     !! m/s to micro-m/s 
  integer                       :: sampling_int(1 : nsta)
  integer                       :: nsec, tim, nch_chtbl
  integer,          allocatable :: waveform_obs_int(:, :), npts_win(:, :)
  type(winch__hdr), allocatable :: chtbl(:)
#endif

  !!OpenMP variable
  !$ integer                    :: omp_thread

  icount = iargc()
#ifdef WIN
  if(icount .ne. 8) then
    write(0, '(a)') "usage: ./a.out winfile win_chfile dem_file_name ot_begin ot_end rms_time_window resultdir result_file_name"
    error stop
  endif
  
  call getarg(1, win_filename)
  call getarg(2, win_chfilename)
  call getarg(3, dem_file)
  call getarg(4, ot_begin_t); read(ot_begin_t, *) ot_begin
  call getarg(5, ot_end_t)  ; read(ot_end_t, *) ot_end
  call getarg(6, rms_tw_t)  ; read(rms_tw_t, *) rms_tw
  call getarg(7, resultdir)
  call getarg(8, resultfile)
#else
  if(icount .ne. 7) then
    write(0, '(a)') "usage: ./a.out sacfile_index dem_file_name ot_begin ot_end rms_time_window resultdir result_file_name"
    error stop
  endif
  
  call getarg(1, sacfile_index)
  call getarg(2, dem_file)
  call getarg(3, ot_begin_t); read(ot_begin_t, *) ot_begin
  call getarg(4, ot_end_t)  ; read(ot_end_t, *) ot_end
  call getarg(5, rms_tw_t)  ; read(rms_tw_t, *) rms_tw
  call getarg(6, resultdir)
  call getarg(7, resultfile)
#endif

  ot_shift = rms_tw / 2.0_fp

  write(0, '(a, 3(1x, f8.3))') "lon_w, lat_s, z_min =", lon_w, lat_s, z_min
  write(0, '(a, 3(1x, i0))') "nlon, nlat, nz =", nlon, nlat, nz

  !!read topography file
  !!pay attention to the order
  open(unit = 10, file = dem_file, form = "unformatted", access = "direct", recl = 12)
  icount = 1
  do j = nlat, 1, -1
    do i = 1, nlon
      read(10, rec = icount) lon_r, lat_r, topo_r
      topography(i, j) = real(topo_r, kind = fp) * alt_to_depth
      icount = icount + 1
    enddo
  enddo
  close(10)

  !!set velocity/attenuation structure
  call set_velocity(z_min, dz, velocity, qinv)

#ifdef WIN
  !!read waveform_obs_int from winfile
  call win__read_file(trim(win_filename), st_winch, sampling_int, nsec, tim, waveform_obs_int, npts_win)
  !!read channel table
  call winch__read_tbl(trim(win_chfilename), chtbl)
  npts_max = ubound(waveform_obs_int, 1)
  nch_chtbl = ubound(chtbl, 1)
  allocate(waveform_obs(1 : npts_max, 1 : nsta))
  waveform_obs(1 : npts_max, 1 : nsta) = 0.0_dp
  do j = 1, nsta
    chtbl_loop: do i = 1, nch_chtbl
      if(chtbl(i)%achid .eq. st_winch(j)) then
        write(0, '(3a)') "chid ", st_winch(j), " found"
        lon_sta(j) = real(chtbl(i)%lon, kind = fp)
        lat_sta(j) = real(chtbl(i)%lat, kind = fp)
        z_sta(j)   = real(chtbl(i)%elev, kind = fp) * alt_to_depth
        npts(j) = nsec * sampling_int(j)
        sampling(j) = 1.0_fp / real(sampling_int(j), kind = fp)
        do ii = 1, npts(j)
          waveform_obs(ii, j) = waveform_obs_int(ii, j) * chtbl(i)%conv * order_um
        enddo
        write(0, '(a, i0, 3a, f8.4, 1x, f7.4, 1x, f6.2)') &
        &     "station(", j, ") name = ", trim(chtbl(i)%stnm), " lon/lat = ", lon_sta(j), lat_sta(j), z_sta(j)
        exit chtbl_loop
      endif
    enddo chtbl_loop
  enddo
  
#else

  !!read sac file
  npts_max = 0
  do j = 1, nsta
    sacfile = trim(sacfile_index) // trim(stname(j)) // "__U__.sac"
    call read_sachdr(sacfile, delta=sampling(j), stlat = lat_sta(j), stlon = lon_sta(j), stdp = z_sta(j), &
    &                npts = npts(j), begin = begin(j))
    if(npts(j) .gt. npts_max) npts_max = npts(j)
    z_sta(j) = z_sta(j) * 0.001_fp
    if(j .ne. 1) then
      do i = 1, j - 1
        if(begin(j) .ne. begin(i)) then
          write(0, '(3a)') "beginning time is different: ", trim(stname(j)), trim(stname(i))
          error stop
        endif
      enddo
    endif
    write(0, '(a, i0, 3a, f8.4, 1x, f7.4, 1x, f6.2)') &
    &     "station(", j, ") name = ", trim(stname(j)), " lon/lat/dep = ", lon_sta(j), lat_sta(j), z_sta(j)
  enddo
  allocate(waveform_obs(npts_max, nsta))
  do i = 1, nsta
    call read_sacdata(sacfile, npts_max, waveform_obs(:, i))
  enddo
  do i = 1, nsta
    waveform_obs(1 : npts(i), i) = waveform_obs(1 : npts(i), i) * order
  enddo
#endif

  !!remove offset
  do j = 1, nsta
    amp_avg = 0.0_dp
    icount = 0
    do i = 1, int(rms_tw / sampling(j))
      amp_avg = amp_avg + waveform_obs(i, j)
      icount = icount + 1
    enddo
    amp_avg = amp_avg / real(icount, kind = dp)
    do i = 1, npts(j)
      waveform_obs(i, j) = (waveform_obs(i, j) - amp_avg)
    enddo
  enddo

  !!bandpass filter
  do j = 1, nsta
    !!design
    call calc_bpf_order(fl, fh, fs, ap, as, sampling(j), m, n, c)
    allocate(h(4 * m), waveform_tmp(npts(j)))
    call calc_bpf_coef(fl, fh, sampling(j), m, n, h, c, gn)

    !!apply
    waveform_tmp(1 : npts(j)) = waveform_obs(1 : npts(j), j)
    call tandem1(waveform_tmp, waveform_tmp, npts(j), h, m, 1)
    waveform_obs(1 : npts(j), j) = waveform_tmp(1 : npts(j)) * gn
    deallocate(h, waveform_tmp)
  enddo

 
  !!make traveltime/pulse width table for each grid point
  write(0, '(a)') "making traveltime / pulse width table..."
  !$omp parallel default(none), &
  !$omp&         shared(topography, lat_sta, lon_sta, z_sta, velocity, qinv, ttime_min, width_min, hypodist), &
  !$omp&         private(i, j, k, ii, jj, lon_grid, lat_grid, depth_grid, az_ini, epdelta, lon_index, lat_index, z_index, &
  !$omp&                dist_min, inc_angle_tmp, lon_tmp, lat_tmp, depth_tmp, az_tmp, &
  !$omp&                ttime_tmp, width_tmp, xgrid, ygrid, zgrid, dist_tmp, val_1d, velocity_interpolate, val_3d, &
  !$omp&                qinv_interpolate, dvdz, lon_new, lat_new, depth_new, az_new, inc_angle_new, omp_thread)

  !$ omp_thread = omp_get_thread_num()

  !$omp do schedule(guided)
  z_loop: do k = 1, nz - 1
    depth_grid = z_min + dz * real(k - 1, kind = fp)
    !$ write(0, '(2(a, i0))') "omp_thread_num = ", omp_thread, " k = ", k
    lat_loop: do j = 1, nlat
      lat_grid = lat_s + dlat * real(j - 1, kind = fp)
      lon_loop: do i = 1, nlon
        lon_grid = lon_w + dlon * real(i - 1, kind = fp)

        !!check the grid is lower than the topo
        if(depth_grid .lt. topography(i, j)) cycle
        !write(0, '(a, (f7.3, 1x, f6.3, 1x, f6.2), a)') "(lon, lat, z) = (", lon_grid, lat_grid, depth_grid, ")"

        !!calculate traveltime to station        
        station_loop: do jj = 1, nsta
          !!calculate azimuth and hypocentral distance
          call greatcircle_dist(lat_grid, lon_grid, lat_sta(jj), lon_sta(jj), azimuth = az_ini, delta_out = epdelta)
          lon_index = int((lon_sta(jj) - lon_w) / dlon) + 1
          lat_index = int((lat_sta(jj) - lat_s) / dlat) + 1
          !print *, lon_sta(jj), lon_w + real(lon_index - 1) * dlon
          !print *, lat_sta(jj), lat_s + real(lat_index - 1) * dlat
          hypodist(jj, i, j, k) = sqrt((r_earth - depth_grid) ** 2 + (r_earth - z_sta(jj)) ** 2 &
          &                     - 2.0_fp * (r_earth - depth_grid) * (r_earth - z_sta(jj)) * cos(epdelta))

#ifdef VEL_CONST
          !!homogeneous structure: using velocity/qinv at the grid
          ttime_min(jj, i, j, k) = hypodist(jj, i, j, k) / velocity(i, j, k)
          width_min(jj, i, j, k) = ttime_min(jj, i, j, k) * qinv(i, j, k)

#else
          
          !!do ray shooting
          dist_min = huge
          ttime_min(jj, i, j, k) = real(huge, kind = fp)
          width_min(jj, i, j, k) = 0.0_fp
          incangle_loop: do ii = 1, ninc_angle - 1
            inc_angle_tmp = real(ii, kind = fp) * dinc_angle
   
            lon_tmp = lon_grid
            lat_tmp = lat_grid
            depth_tmp = depth_grid
            az_tmp = az_ini

            ttime_tmp = 0.0_fp
            width_tmp = 0.0_fp
            !!loop until ray arrives at surface/boundary
            shooting_loop: do
              lon_index = int((lon_tmp - lon_w) / dlon) + 1
              lat_index = int((lat_tmp - lat_s) / dlat) + 1
              z_index   = int((depth_tmp - z_min) / dz) + 1

              !!exit if ray approaches to the boundary
              if(lon_index .lt. 1 .or. lon_index .gt. nlon - 1      &
              &  .or. lat_index .lt. 1 .or. lat_index .gt. nlat - 1 &
              &  .or. z_index .lt. 1 .or. z_index .gt. nz - 1) then
                exit shooting_loop
              endif

              xgrid(1) = lon_w + real(lon_index - 1, kind = fp) * dlon; xgrid(2) = xgrid(1) + dlon
              ygrid(1) = lat_s + real(lat_index - 1, kind = fp) * dlat; ygrid(2) = ygrid(1) + dlat
              zgrid(1) = z_min + real(z_index - 1, kind = fp) * dz;     zgrid(2) = zgrid(1) + dz

              !!calculate distance between ray and station
              call greatcircle_dist(lat_tmp, lon_tmp, lat_sta(jj), lon_sta(jj), delta_out = epdelta)
              dist_tmp = sqrt((r_earth - depth_tmp) ** 2 + (r_earth - z_sta(jj)) ** 2 &
              &        - 2.0_fp * (r_earth - depth_tmp) * (r_earth - z_sta(jj)) * cos(epdelta))
              !print *, real(ii, kind = fp) * dinc_angle, lat_tmp, lon_tmp, depth_tmp
              !print *, jj, lat_sta(jj), lon_sta(jj), z_sta(jj), dist_tmp

              if(dist_tmp .lt. dist_min) then
                dist_min = dist_tmp
                ttime_min(jj, i, j, k) = ttime_tmp
                width_min(jj, i, j, k) = width_tmp
              endif
 
              !!shooting the ray
              val_1d(1 : 2) = velocity(lon_index, lat_index, z_index : z_index + 1)
              call linear_interpolation_1d(depth_tmp, zgrid, val_1d, velocity_interpolate)
              val_3d(1 : 2, 1 : 2, 1 : 2) = qinv(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
              call block_interpolation_3d(lon_tmp, lat_tmp, depth_tmp, xgrid, ygrid, zgrid, val_3d, qinv_interpolate)

              dvdz = (val_1d(2) - val_1d(1)) / dz
              call rayshooting3D(lon_tmp, lat_tmp, depth_tmp, az_tmp, inc_angle_tmp, time_step, velocity_interpolate, &
              &                  dvdlon, dvdlat, dvdz, lon_new, lat_new, depth_new, az_new, inc_angle_new)
              ttime_tmp = ttime_tmp + time_step
              width_tmp = width_tmp + qinv_interpolate * time_step

              lon_tmp = lon_new
              lat_tmp = lat_new
              depth_tmp = depth_new
              az_tmp = az_new
              inc_angle_tmp = inc_angle_new
            enddo shooting_loop

          enddo incangle_loop
          !print *, "station lon, lat, dist = ", lon_sta(jj), lat_sta(jj), dist_min, hypodist(jj), ttime_min, width_min(jj)
#endif
        enddo station_loop
      enddo lon_loop
    enddo lat_loop
  enddo z_loop
  !$omp end do
  !$omp end parallel

  !!find minimum residual grid for seismic source 
  resultfile = trim(resultdir) // "/" // trim(resultfile)
  open(unit = 10, file = resultfile)
  write(10, '(a)') "# OT min_lon min_lat min_dep source_amp residual"
  residual(1 : nlon, 1 : nlat, 1 : nz) = huge

  time_count = 0
  time_loop: do
    origintime = ot_begin + ot_shift * real(time_count, kind = fp)
    if(origintime .gt. ot_end) exit time_loop
    !write(0, '(a, f6.1)') "origin time = ", origintime

    source_amp(1 : nlon, 1 : nlat, 1 : nz) = 0.0_dp
    residual(1 : nlon, 1 : nlat, 1 : nz) = huge
    z_loop2: do k = 1, nz - 1
      depth_grid = z_min + dz * real(k - 1, kind = fp)
      lat_loop2: do j = 1, nlat
        lon_loop2: do i = 1, nlon

          if(depth_grid .lt. topography(i, j)) cycle

          !!caluculate site-corrected amplitude
          do jj = 1, nsta
            if(ttime_min(jj, i, j, k) .eq. real(huge, kind = fp)) then
              cycle lon_loop2
            endif

            wave_index = int((origintime + ttime_min(jj, i, j, k)) / sampling(jj) + 0.5_fp) + 1
            if(wave_index .gt. npts(jj)) then
              write(0, '(a)') "wave_index is larger than npts"
              error stop
            endif
            rms_amp_obs(jj) = 0.0_fp
            icount = 0
            do ii = wave_index, wave_index + int(rms_tw / sampling(jj) + 0.5_fp) - 1
              if(ii .lt. npts(jj)) then
                rms_amp_obs(jj) = rms_amp_obs(jj) + waveform_obs(ii, jj) * waveform_obs(ii, jj)
                icount = icount + 1
              endif
            enddo
            rms_amp_obs(jj) = sqrt(rms_amp_obs(jj) / real(icount, kind = dp))
  
            !!calculate source amplitude and residual
            source_amp(i, j, k) = source_amp(i, j, k) &
            &            + rms_amp_obs(jj) / siteamp(jj) &
            &              * real(hypodist(jj, i, j, k) * exp(width_min(jj, i, j, k) * (pi * freq)), kind = dp)
          enddo
          source_amp(i, j, k) = source_amp(i, j, k) / real(nsta, kind = dp)
          residual(i, j, k) = 0.0_dp
          residual_normalize = 0.0_dp
          do ii = 1, nsta
            residual(i, j, k) = residual(i, j, k) &
            &                 + (rms_amp_obs(ii) / siteamp(ii) &
            &                 - source_amp(i, j, k) / real(hypodist(ii, i, j, k), kind = dp) &
            &                   * real(exp(-pi * freq * width_min(ii, i, j, k)), kind = dp)) ** 2
            residual_normalize = residual_normalize + (rms_amp_obs(ii) / siteamp(ii)) ** 2
          enddo
          residual(i, j, k) = residual(i, j, k) / residual_normalize
  
        enddo lon_loop2
      enddo lat_loop2
    enddo z_loop2

    !!search minimum residual
    residual_minloc = minloc(residual)
    lon_grid = lon_w + real(residual_minloc(1) - 1, kind = fp) * dlon
    lat_grid = lat_s + real(residual_minloc(2) - 1, kind = fp) * dlat
    depth_grid = z_min + real(residual_minloc(3) - 1, kind = fp) * dz
    write(0, '(a, f0.1, a, f0.3, 1x, f0.3, 1x, f0.2, a, 2(a, e15.7))') &
    &                 "OT = ", origintime, " residual_minimum (lon, lat, dep) = (", lon_grid, lat_grid, depth_grid, ")", &
    &                 " source_amp = ", source_amp(residual_minloc(1), residual_minloc(2), residual_minloc(3)), &
    &                 " residual = ", residual(residual_minloc(1), residual_minloc(2), residual_minloc(3))
    write(10, '(f0.1, 1x, f0.3, 1x, f0.3, 1x, f0.2, 2(1x, e15.7))') &
    &                 origintime, lon_grid, lat_grid, depth_grid, &
    &                 source_amp(residual_minloc(1), residual_minloc(2), residual_minloc(3)), &
    &                 residual(residual_minloc(1), residual_minloc(2), residual_minloc(3))

    !!output grd files
    call int_to_char(time_count, maxlen, time_count_char)
    !!output horizontal slice
    grdfile = trim(resultdir) // "/min_err_lon-lat_" // trim(time_count_char) // ".grd"
    allocate(residual_grd(nlon, nlat))
    xrange(1) = real(lon_w, kind = dp)
    xrange(2) = real(lon_w, kind = dp) + real(nlon - 1, kind = dp) * dlon
    yrange(1) = real(lat_s, kind = dp)
    yrange(2) = real(lat_s, kind = dp) + real(nlat - 1, kind = dp) * dlat
    spacing(1 : 2) = [real(dlon, kind = dp), real(dlat, kind = dp)]
    residual_grd(1 : nlon, 1 : nlat) = real(residual(1 : nlon, 1 : nlat, residual_minloc(3)), kind = sp)
    grd_status = grd_create(grdfile, residual_grd, xrange, yrange, spacing, jscan = 1, &
    &                       NaN = real(huge, kind = sp), overwrite = .true.)
    deallocate(residual_grd)

    !!output depth slice -- lon-dep
    grdfile = trim(resultdir) // "/min_err_lon-dep_" // trim(time_count_char) // ".grd"
    allocate(residual_grd(nlon, nz))
    xrange(1) = real(lon_w, kind = dp)
    xrange(2) = real(lon_w, kind = dp) + real(nlon - 1, kind = dp) * dlon
    yrange(1) = real(z_min, kind = dp)
    yrange(2) = real(z_min, kind = dp) + real(nz - 1, kind = dp) * dz
    spacing(1 : 2) = [real(dlon, kind = dp), real(dz, kind = dp)]
    residual_grd(1 : nlon, 1 : nz) = real(residual(1 : nlon, residual_minloc(2), 1 : nz), kind = sp)
    grd_status = grd_create(grdfile, residual_grd, xrange, yrange, spacing, jscan = 1, &
    &                       NaN = real(huge, kind = sp), overwrite = .true.)
    deallocate(residual_grd)

    !!output depth slice -- dep-lat
    grdfile = trim(resultdir) // "/min_err_dep-lat_" // trim(time_count_char) // ".grd"
    allocate(residual_grd(nz, nlat))
    xrange(1) = real(z_min, kind = dp)
    xrange(2) = real(z_min, kind = dp) + real(nz - 1, kind = dp) * dz
    yrange(1) = real(lat_s, kind = dp)
    yrange(2) = real(lat_s, kind = dp) + real(nlat - 1, kind = dp) * dlat
    spacing(1 : 2) = [real(dz, kind = dp), real(dlat, kind = dp)]
    do j = 1, nz
      do i = 1, nlat
        residual_grd(j, i) = real(residual(residual_minloc(1), i, j), kind = sp)
      enddo
    enddo
    grd_status = grd_create(grdfile, residual_grd, xrange, yrange, spacing, jscan = 1, &
    &                       NaN = real(huge, kind = sp), overwrite = .true.)
    deallocate(residual_grd)

    time_count = time_count + 1
  enddo time_loop
  close(10)

  stop
end program AmplitudeSourceLocation_PulseWidth

