program AmplitudeSourceLocation_PulseWidth
  !!Amplitude Source Location using depth-dependent 1D velocity structure, 3D heterogeneous attenuation structure
  !!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2020
  !!License  : MIT License (https://opensource.org/licenses/MIT)

#if ! defined (WIN) && ! defined (SAC) && ! defined (AMP_TXT)  /* define default waveform data format */
#define SAC
#endif

  use nrtype,               only : fp, sp, dp
  use constants,            only : rad2deg, deg2rad, pi, r_earth
  use rayshooting,          only : rayshooting3D
  use set_velocity_model,   only : set_velocity
  use linear_interpolation, only : linear_interpolation_1d, linear_interpolation_2d, block_interpolation_3d
  use greatcircle,          only : greatcircle_dist
  use itoa,                 only : int_to_char
  use grdfile_io,           only : read_grdfile_2d, write_grdfile_2d
#if defined (WIN)
  use m_win
  use m_winch
#endif
#if defined (SAC)
  use read_sacfile,         only : read_sachdr, read_sacdata
#endif
#if defined (RAYBENDING)
  use raybending,           only : pseudobending3D
#endif
  !$ use omp_lib

  implicit none
  integer,                parameter :: wavetype = 2           !!1 for P-wave, 2 for S-wave
  integer,                parameter :: nsta_use_minimum = 8
  real(kind = dp),        parameter :: snratio_accept = 0.0_dp
  real(kind = fp),        parameter :: do_rayshooting_threshold = 300.0_fp
  real(kind = fp),        parameter :: conv_rayshooting_threshold = 0.1_fp
  !!Use station
  integer                           :: nsta
  integer,              allocatable :: npts(:)
  character(len = 6),   allocatable :: stname(:)
  real(kind = dp),      allocatable :: siteamp(:)
  real(kind = fp),      allocatable :: ttime_cor(:, :)
  logical,              allocatable :: use_flag(:), use_flag_tmp(:)
#if defined (WIN)    /* use win-format waveform file for input waveforms */
  character(len = 129)              :: win_filename, win_chfilename
  character(len = 10)               :: cmpnm
  real(dp),               parameter :: order = 1.0e+6_dp                       !! m/s to micro-m/s 
  integer                           :: nsec, tim, nch_chtbl
  integer,              allocatable :: waveform_obs_int(:, :), npts_win(:, :), sampling_int(:), ikey(:)
  character(len = 4),   allocatable :: st_winch(:)
  type(winch__hdr),     allocatable :: chtbl(:)
#elif defined (SAC)         /* use sac binary files as input waveforms */
  real(kind = dp),        parameter :: order = 1.0e+6_dp
  character(len = 3),     parameter :: sacfile_extension = "sac"
  character(len = 10)               :: cmpnm
  character(len = 129)              :: sacfile, sacfile_index
#elif defined (AMP_TXT)
  integer                           :: namp
  character(len = 129)              :: amp_filename, freq_t
  character(len = 129), allocatable :: eventindex(:)
  real(kind = dp),      allocatable :: amp_txt(:, :)
#endif

  !!Search range
  real(kind = fp),        parameter :: lon_w = 135.5_fp, lon_e = 137.5_fp
  real(kind = fp),        parameter :: lat_s = 32.7_fp, lat_n = 33.7_fp
  real(kind = fp),        parameter :: z_min = 0.0_fp, z_max = 20.0_fp
  real(kind = fp),        parameter :: dlon = 0.02_fp, dlat = 0.02_fp, dz = 2.0_fp
  !!structure range
  real(kind = fp),        parameter :: lon_str_w = 135.0_fp, lon_str_e = 138.0_fp
  real(kind = fp),        parameter :: lat_str_s = 32.0_fp,  lat_str_n = 34.5_fp
  real(kind = fp),        parameter :: z_str_min = 0.0_fp, z_str_max = 30.0_fp
  real(kind = fp),        parameter :: dlon_str = 0.005_fp, dlat_str = 0.005_fp, dz_str = 0.1_fp
  !!Ray shooting
  real(kind = fp),        parameter :: dvdlon = 0.0_fp, dvdlat = 0.0_fp         !!assume 1D structure
  integer,                parameter :: ninc_angle = 180                         !!grid search in incident angle
  integer,                parameter :: nrayshoot = 2                            !!number of grid search
  real(kind = fp),        parameter :: time_step = 0.01_fp

  integer,                parameter :: maxlen = 4
  real(kind = fp),        parameter :: alt_to_depth = -1.0e-3_fp
  real(kind = dp),        parameter :: huge = 1.0e+5_dp

  integer,                parameter :: nlon = int((lon_e - lon_w) / dlon) + 2
  integer,                parameter :: nlat = int((lat_n - lat_s) / dlat) + 2
  integer,                parameter :: nz   = int((z_max - z_min) / dz) + 2
  integer,                parameter :: nlon_str = int((lon_str_e - lon_str_w) / dlon_str) + 2
  integer,                parameter :: nlat_str = int((lat_str_n - lat_str_s) / dlat_str) + 2
  integer,                parameter :: nz_str   = int((z_str_max - z_str_min) / dz_str) + 2

  real(kind = fp)                   :: velocity(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                    qinv(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                    val_1d(1 : 2), val_2d(1 : 2, 1 : 2), val_3d(1 : 2, 1 : 2, 1 : 2), &
  &                                    xgrid(1 : 2), ygrid(1 : 2), zgrid(1 : 2), inc_angle_ini_min(0 : nrayshoot)
  real(kind = dp)                   :: residual(1 : nlon, 1 : nlat, 1 : nz), source_amp(1 : nlon, 1 : nlat, 1 : nz)
  integer                           :: residual_minloc(3), nsta_use_grid(1 : nlon, 1 : nlat, 1 : nz)
  real(kind = sp),      allocatable :: residual_grd(:, :)
  real(kind = fp),      allocatable :: lon_sta(:), lat_sta(:), z_sta(:), sampling(:), begin(:), &
  &                                    width_min(:, :, :, :), ttime_min(:, :, :, :), hypodist(:, :, :, :)
  real(kind = dp),      allocatable :: waveform_obs(:, :), topography(:, :), lon_topo(:), lat_topo(:), rms_amp_obs(:), &
  &                                    rms_amp_obs_noise(:)
  
  real(kind = fp)                   :: ttime_tmp, width_tmp, velocity_interpolate, qinv_interpolate, &
  &                                    az_tmp, inc_angle_tmp, inc_angle_min, az_new, inc_angle_new, az_ini, &
  &                                    lon_tmp, lat_tmp, depth_tmp, lon_new, lat_new, depth_new, origintime, &
  &                                    lon_grid, lat_grid, depth_grid, dist_min, dist_tmp, dvdz, epdelta, &
  &                                    ot_begin, ot_end, ot_shift, rms_tw, lon_min, lat_min, depth_min, &
  &                                    dinc_angle, dinc_angle_org, inc_angle_ini
  real(kind = sp)                   :: lon_r, lat_r, topo_r
  real(kind = dp)                   :: residual_normalize, amp_avg, topography_interpolate, &
  &                                    dlon_topo, dlat_topo, freq, rms_amp_cal
  
  integer                           :: i, j, k, ii, jj, kk, icount, wave_index, time_count, lon_index, lat_index, z_index, &
  &                                    npts_max, nlon_topo, nlat_topo, ios
  character(len = 129)              :: station_param, dem_file, ot_begin_t, ot_end_t, &
  &                                    rms_tw_t, ot_shift_t, grdfile, resultfile, resultdir, ampfile
  character(len = maxlen)           :: time_count_char

#if defined (RAYBENDING)
  integer,                parameter :: ndiv_raypath = 10
  integer                 parameter :: nraypath_ini = 4
  real(kind = fp)                   :: raypath_lon(nraypath_ini + 2 * (2 ** ndiv_raypath - 1), &
  &                                    raypath_lat(nraypath_ini + 2 * (2 ** ndiv_raypath - 1), &
  &                                    raypath_dep(nraypath_ini + 2 * (2 ** ndiv_raypath - 1)
  integer                           :: nraypath

#if defined (AMP_RATIO)
  real(kind = dp)                   :: amp_ratio_obs, amp_ratio_cal
#endif
#if defined

  !!filter variables
  real(kind = dp),        parameter :: ap = 0.5_dp, as = 5.0_dp                 !!bandpass filter parameters
  real(kind = dp),      allocatable :: h(:), waveform_tmp(:)
  real(kind = dp)                   :: gn, c, fl, fh, fs
  integer                           :: m, n
  character(len = 6)                :: fl_t, fh_t, fs_t

  !!OpenMP variable
  !$ integer                        :: omp_thread

  icount = iargc()

#if defined (AMP_TXT)
  if(icount .ne. 6) then
    write(0, '(a)', advance="no") "usage: ./asl_pw (topography_grd) (station_param_file) (txtfile_amplitude) (frequency)"
    write(0, '(a)')               " (result_dir) (result_file_name)"
    error stop
  endif
  
  call getarg(1, dem_file)
  call getarg(2, station_param)
  call getarg(3, amp_filename)
  call getarg(4, freq_t); read(freq_t, *) freq
  call getarg(5, resultdir)
  call getarg(6, resultfile)
#elif defined (WIN)
  if(icount .ne. 14) then
    write(0, '(a)', advance="no") "usage: ./asl_pw (topography_grd) (station_param_file) (winfile) (win_chfile) (component_name)"
    write(0, '(a)', advance="no") " (fl) (fh) (fs) (ot_begin) (ot_end) (ot_shift) (rms_time_window_length) (result_dir)"
    write(0, '(a)')               " (result_file_name)"
    error stop
  endif
  
  call getarg(1, dem_file)
  call getarg(2, station_param)
  call getarg(3, win_filename)
  call getarg(4, win_chfilename)
  call getarg(5, cmpnm)
  call getarg(6, fl_t)       ; read(fl_t, *) fl
  call getarg(7, fh_t)       ; read(fh_t, *) fh
  call getarg(8, fs_t)       ; read(fs_t, *) fs
  call getarg(9, ot_begin_t) ; read(ot_begin_t, *) ot_begin
  call getarg(10, ot_end_t)  ; read(ot_end_t, *) ot_end
  call getarg(11, ot_shift_t); read(ot_shift_t, *) ot_shift
  call getarg(12, rms_tw_t)  ; read(rms_tw_t, *) rms_tw
  call getarg(13, resultdir)
  call getarg(14, resultfile)
#elif defined (SAC)
  if(icount .ne. 13) then
    write(0, '(a)', advance="no") "usage: ./asl_pw (topography_grd) (station_param_file) (sacfile_prefix) (compnent_name)"
    write(0, '(a)', advance="no") " (fl) (fh) (fs) (ot_begin) (ot_end) (ot_shift) (rms_time_window_length) (result_dir)"
    write(0, '(a)')               " (result_file_name)"
    error stop
  endif
  
  call getarg(1, dem_file)
  call getarg(2, station_param)
  call getarg(3, sacfile_index)
  call getarg(4, cmpnm)
  call getarg(5, fl_t)       ; read(fl_t, *) fl
  call getarg(6, fh_t)       ; read(fh_t, *) fh
  call getarg(7, fs_t)       ; read(fs_t, *) fs
  call getarg(8, ot_begin_t) ; read(ot_begin_t, *) ot_begin
  call getarg(9, ot_end_t)   ; read(ot_end_t, *) ot_end
  call getarg(10, ot_shift_t); read(ot_shift_t, *) ot_shift
  call getarg(11, rms_tw_t)  ; read(rms_tw_t, *) rms_tw
  call getarg(12, resultdir)
  call getarg(13, resultfile)
#endif

  write(0, '(a, 3(1x, f8.3))') "lon_w, lat_s, z_min =", lon_w, lat_s, z_min
  write(0, '(a, 3(1x, i0))') "nlon, nlat, nz =", nlon, nlat, nz

#if defined (WIN) || defined (SAC)
  write(0, '(a, 3(f5.2, 1x))') "Bandpass filter parameter fl, fh, fs (Hz) = ", fl, fh, fs
  freq = (fl + fh) * 0.5_dp
#endif

  !!read topography file (netcdf grd format)
  call read_grdfile_2d(dem_file, lon_topo, lat_topo, topography)
  nlon_topo = ubound(lon_topo, 1)
  nlat_topo = ubound(lat_topo, 1)
  dlon_topo = lon_topo(2) - lon_topo(1)
  dlat_topo = lat_topo(2) - lat_topo(1)
  topography(1 : nlon_topo, 1 : nlat_topo) = topography(1 : nlon_topo, 1 : nlat_topo) * alt_to_depth

  !!set velocity/attenuation structure
  call set_velocity(z_str_min, dz_str, velocity, qinv)

  write(0, '(2a)') "read station paramters from ", trim(station_param)
  open(unit = 40, file = station_param)
  nsta = 0
  do
    read(40, *, iostat = ios)
    if(ios .ne. 0) exit
    nsta = nsta + 1
  enddo
  rewind(40)
  allocate(lon_sta(1 : nsta), lat_sta(1 : nsta), z_sta(1 : nsta), stname(1 : nsta), rms_amp_obs(1 : nsta), &
  &        ttime_cor(1 : nsta, 1 : 2), siteamp(1 : nsta), use_flag(1 : nsta), use_flag_tmp(1 : nsta), &
  &        hypodist(1 : nsta, 1 : nlon, 1 : nlat, 1 : nz), ttime_min(1 : nsta, 1 : nlon, 1 : nlat, 1 : nz), &
  &        width_min(1 : nsta, 1 : nlon, 1 : nlat, 1 : nz))
  do i = 1, nsta
    read(40, *) lon_sta(i), lat_sta(i), z_sta(i), stname(i), use_flag(i), ttime_cor(i, 1), ttime_cor(i, 2), siteamp(i)

#if defined (TESTDATA)
    ttime_cor(i, 1) = 0.0_fp
    ttime_cor(i, 2) = 0.0_fp
    siteamp(i)   = 1.0_dp
#endif

  enddo
  close(40)

#if defined (AMP_TXT)
  open(unit = 40, file = amp_filename)
  read(40, *)
  namp = 0
  do
    read(40, *, iostat = ios)
    if(ios .ne. 0) exit
    namp = namp + 1
  enddo
  rewind(40)
  read(40, *)
  allocate(amp_txt(1 : nsta, 1 : namp), eventindex(1 : namp))
  do j = 1, namp
    read(40, *) (amp_txt(i, j), i = 1, nsta), eventindex(j)
  enddo
  close(40)
  write(0, '(3a, i0)') "read amplitude data from ", trim(amp_filename), " namp = ", namp
#else 
#if defined (WIN)
  allocate(st_winch(1 : nsta), sampling_int(1 : nsta), ikey(1 : nsta), &
  &        sampling(1 : nsta), begin(1 : nsta), npts(1 : nsta))
  !!read channel_table
  call winch__read_tbl(trim(win_chfilename), chtbl)
  !!find chid from given stname
  do i = 1, nsta
    call winch__st2chid(chtbl, stname(i), trim(cmpnm), st_winch(i), ikey(i))
    if(ikey(i) .eq. 0) then
      write(0, '(5(a, 1x))') "station/comp =", trim(stname(i)), trim(cmpnm), "does not exist in", trim(win_chfilename)
      error stop
    endif
    write(0, '(6a)') "station name = ", trim(stname(i)), " comp = ", trim(cmpnm), " chid = ", st_winch(i)
  enddo

  !!read waveform_obs_int from winfile
  call win__read_file(trim(win_filename), st_winch, sampling_int, nsec, tim, waveform_obs_int, npts_win)
  npts_max = ubound(waveform_obs_int, 1)
  nch_chtbl = ubound(chtbl, 1)
  begin(1 : nsta) = 0.0_fp
  allocate(waveform_obs(1 : npts_max, 1 : nsta))
  waveform_obs(1 : npts_max, 1 : nsta) = 0.0_dp
  do j = 1, nsta
    npts(j) = nsec * sampling_int(j)
    sampling(j) = 1.0_fp / real(sampling_int(j), kind = fp)
    if(npts(j) .eq. 0) then
      use_flag(j) = .false.
    else
      do i = 1, npts(j)
        waveform_obs(i, j) = waveform_obs_int(i, j) * chtbl(ikey(j))%conv * order
      enddo
    endif
  enddo
  deallocate(sampling_int, chtbl, waveform_obs_int, ikey, st_winch)
#elif defined (SAC)
  !!read sac file
  allocate(sampling(1 : nsta), begin(1 : nsta), npts(1 : nsta))
  npts_max = 0
  do i = 1, nsta
    sacfile = trim(sacfile_index) // "." // trim(stname(i)) // "." // trim(cmpnm) // "." // sacfile_extension
    call read_sachdr(sacfile, delta = sampling(i), npts = npts(i), begin = begin(i))
    if(npts(i) .gt. npts_max) npts_max = npts(i)
    !write(0, '(a, i0, 3a, f8.4, 1x, f7.4, 1x, f6.2)') &
    !&     "station(", i, ") name = ", trim(stname(i)), " lon/lat/dep = ", lon_sta(i), lat_sta(i), z_sta(i)
  enddo
  allocate(waveform_obs(npts_max, nsta))
  do i = 1, nsta
    sacfile = trim(sacfile_index) // "." // trim(stname(i)) // "." // trim(cmpnm) // "." // sacfile_extension
    call read_sacdata(sacfile, npts_max, waveform_obs(:, i))
    waveform_obs(1 : npts(i), i) = waveform_obs(1 : npts(i), i) * order
  enddo
#endif   /* -DWIN or -DSAC */

  !!remove offset
  write(0, '(a)') "Removing offset"
  do j = 1, nsta
    if(use_flag(j) .eqv. .false.) cycle
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

#if defined (TESTDATA)
#else
  !!bandpass filter
  write(0, '(a, 3(f5.2, a))') "Applying bandpass filter fl = ", fl, " (Hz), fh = ", fh, " (Hz), fs = ", fs, " (Hz)"
  do j = 1, nsta
    if(use_flag(j) .eqv. .false.) cycle
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
#endif   /* -DTESTDATA */  

#endif   /* -DAMP_TXT */

  !!calculate noise level
  allocate(rms_amp_obs_noise(1 : nsta))
#if defined (AMP_TXT)
  rms_amp_obs_noise(1 : nsta) = 1.0_dp / real(huge, kind = dp)
#else
  do j = 1, nsta
    icount = 0
    rms_amp_obs_noise(j) = 0.0_dp
    if(use_flag(j) .eqv. .false.) cycle
    do i = int(rms_tw / sampling(j)) + 1, int(rms_tw / sampling(j)) * 2
      rms_amp_obs_noise(j) = rms_amp_obs_noise(j) + waveform_obs(i, j) ** 2
      icount = icount + 1
    enddo
    rms_amp_obs_noise(j) = sqrt(rms_amp_obs_noise(j) / real(icount, kind = dp))
  enddo
#endif


  !!station list
  do i = 1, nsta
    write(0, '(a, i0, a, f9.4, a, f8.4, a, f6.3, 1x, g0)') &
    &     "station(", i, ") lon(deg) = ", lon_sta(i), " lat(deg) = ", lat_sta(i), " depth(km) = ", z_sta(i), use_flag(i)
  enddo

#if defined (READ_TTIMETABLE)
  open(unit = 10, file = "ttime_pulsewidth_table.dat", form = "unformatted", access = "direct", recl = fp, iostat = ios)
  if(ios .ne. 0) then
    close(10)
    write(0, '(a)') "error: cannot read ttime_pulsewidth_table.dat"
    error stop
  else
    icount = 1
    do k = 1, nz
      do j = 1, nlat
        do i = 1, nlon
          do ii = 1, nsta
            read(10, rec = icount) ttime_min(ii, i, j, k)
            icount = icount + 1
            read(10, rec = icount) width_min(ii, i, j, k)
            icount = icount + 1
            read(10, rec = icount) hypodist(ii, i, j, k)
            icount = icount + 1
          enddo
        enddo
      enddo
    enddo
    close(10)
  endif

#else /* -DREAD_TTIMETABLE */

  !!make traveltime/pulse width table for each grid point
  write(0, '(a)') "making traveltime / pulse width table..."
  !$omp parallel default(none), &
  !$omp&         shared(topography, nsta, stname, lat_sta, lon_sta, z_sta, velocity, qinv, ttime_min, width_min, hypodist, &
  !$omp&                lon_topo, lat_topo, dlon_topo, dlat_topo, nlon_topo, nlat_topo, use_flag), &
  !$omp&         private(i, j, k, ii, jj, kk, lon_grid, lat_grid, depth_grid, az_ini, epdelta, lon_index, lat_index, z_index, &
  !$omp&                dist_min, inc_angle_tmp, inc_angle_min, lon_tmp, lat_tmp, depth_tmp, az_tmp, val_2d, &
  !$omp&                topography_interpolate, ttime_tmp, width_tmp, xgrid, ygrid, zgrid, dist_tmp, val_1d, &
  !$omp&                velocity_interpolate, val_3d, qinv_interpolate, dvdz, lon_new, lat_new, depth_new, &
  !$omp&                az_new, inc_angle_new, omp_thread, lon_min, lat_min, depth_min, dinc_angle, dinc_angle_org, &
#if defined (RAYBENDING)
  !$omp&                raypath_lon, raypath_lat, raypath_dep, nraypath, &
#endif
  !$omp&                inc_angle_ini, inc_angle_ini_min)

  !$ omp_thread = omp_get_thread_num()

  !$omp do schedule(guided)
  z_loop: do k = 1, nz
  !z_loop: do k = nz / 2, nz / 2
    depth_grid = z_min + dz * real(k - 1, kind = fp)
    !$ write(0, '(2(a, i0))') "omp_thread_num = ", omp_thread, " depth index k = ", k
    lat_loop: do j = 1, nlat
    !lat_loop: do j = nlat / 2, nlat / 2
      lat_grid = lat_s + dlat * real(j - 1, kind = fp)
      lon_loop: do i = 1, nlon
      !lon_loop: do i = nlon / 2, nlon / 2
        lon_grid = lon_w + dlon * real(i - 1, kind = fp)

        ttime_min(1 : nsta, i, j, k) = real(huge, kind = fp)
        width_min(1 : nsta, i, j, k) = 0.0_fp

        !!check the grid is lower than the topo
        lon_index = int((lon_grid - lon_topo(1)) / dlon_topo) + 1
        lat_index = int((lat_grid - lat_topo(1)) / dlat_topo) + 1
        xgrid(1 : 2) = [lon_topo(lon_index), lon_topo(lon_index + 1)]
        ygrid(1 : 2) = [lat_topo(lat_index), lat_topo(lat_index + 1)]
        val_2d(1 : 2, 1 : 2) = topography(lon_index : lon_index + 1, lat_index : lat_index + 1)
        call linear_interpolation_2d(lon_grid, lat_grid, xgrid, ygrid, val_2d, topography_interpolate)
        if(depth_grid .lt. topography_interpolate) cycle lon_loop
        !write(0, '(a, (f7.3, 1x, f6.3, 1x, f6.2), a)') "(lon, lat, z) = (", lon_grid, lat_grid, depth_grid, ")"

        !!calculate traveltime to station        
        station_loop: do jj = 1, nsta
          !!calculate azimuth and hypocentral distance
          call greatcircle_dist(lat_grid, lon_grid, lat_sta(jj), lon_sta(jj), azimuth = az_ini, delta_out = epdelta)
          hypodist(jj, i, j, k) = sqrt((r_earth - depth_grid) ** 2 + (r_earth - z_sta(jj)) ** 2 &
          &                     - 2.0_fp * (r_earth - depth_grid) * (r_earth - z_sta(jj)) * cos(epdelta))

          if(use_flag(jj) .eqv. .false.) cycle

#if defined (V_CONST)
          !!homogeneous structure: using velocity/qinv at the grid
          lon_index = int((lon_grid - lon_str_w) / dlon_str) + 1
          lat_index = int((lat_grid - lat_str_s) / dlat_str) + 1
          z_index   = int((depth_grid - z_str_min) / dz_str) + 1
          ttime_min(jj, i, j, k) = hypodist(jj, i, j, k) / velocity(lon_index, lat_index, z_index, wavetype)
          width_min(jj, i, j, k) = ttime_min(jj, i, j, k) * qinv(lon_index, lat_index, z_index, wavetype)

#else

#if defined (RAYBENDING)
          !!do ray tracing with pseudobending scheme
          raypath_lon(1) = lon_grid
          raypath_lat(1) = lat_grid
          raypath_dep(1) = depth_grid
          raypath_lon(nraypath_ini) = lon_sta(jj)
          raypath_lat(nraypath_ini) = lat_sta(jj)
          raypath_dep(nraypath_ini) = z_sta(jj)
          do ii = 2, nraypath_ini - 1
            raypath_lon(ii) = raypath_lon(1) + (raypath_lon(nraypath_ini) - raypath_lon(1)) * dble(ii - 1)
            raypath_lat(ii) = raypath_lat(1) + (raypath_lat(nraypath_ini) - raypath_lat(1)) * dble(ii - 1)
            raypath_dep(ii) = raypath_dep(1) + (raypath_dep(nraypath_ini) - raypath_dep(1)) * dble(ii - 1)
          enddo
          nraypath = nraypath_ini
          call pseudobending3D(raypath_lon, raypath_lat, raypath_dep, nraypath, ndiv_raypath, &
          &                    ttime_min(jj, i, j, k), &
          &                    velocity(:, :, :, wavetype), lon_str_w, lat_str_s, z_str_min, dlon_str, dlat_str, dz_str, &
          &                    qinv = qinv(:, :, :, wavetype), lon_w_qinv = lon_str_w, lat_s_qinv = lat_str_s, &
          &                    dep_min_qinv = z_str_min, dlon_qinv = dlon_str, dlat_qinv = dlat_str, ddep_qinv = dz_str, &
          &                    pulsewidth = width_min(jj, i, j, k)

#else

          !!do not rayshooting if hypodist is longer than threshold
          if(hypodist(jj, i, j, k) .gt. do_rayshooting_threshold) cycle station_loop
          !!do ray shooting
          dist_min = huge
          incangle_loop2: do kk = 1, nrayshoot
            if(kk .eq. 1) then
              dinc_angle_org = pi / 2.0_fp
            else
              dinc_angle_org = dinc_angle
            endif
            dinc_angle = 2.0_fp * dinc_angle_org / real(ninc_angle, kind = fp)
            inc_angle_ini_min(0) = dinc_angle_org
            !print *, "dinc_angle = ", dinc_angle * rad2deg, inc_angle_ini_min(kk - 1) * rad2deg

            incangle_loop: do ii = 1, ninc_angle
              inc_angle_ini = (inc_angle_ini_min(kk - 1) - dinc_angle_org) + real(ii, kind = fp) * dinc_angle
              !print '(2(i0, 1x), a, e15.7)', ii, kk, "inc_angle_ini = ", inc_angle_ini * rad2deg
                
              lon_tmp = lon_grid
              lat_tmp = lat_grid
              depth_tmp = depth_grid
              az_tmp = az_ini
              inc_angle_tmp = inc_angle_ini

              ttime_tmp = 0.0_fp
              width_tmp = 0.0_fp
              !!loop until ray arrives at surface/boundary
              shooting_loop: do

                !!exit if ray approaches to the surface
                lon_index = int((lon_tmp - lon_topo(1)) / dlon_topo) + 1
                lat_index = int((lat_tmp - lat_topo(1)) / dlat_topo) + 1
                !!exit if ray approaches to the boundary of the topography array
                if(lon_index .lt. 1 .or. lon_index .gt. nlon_topo - 1      &
                &  .or. lat_index .lt. 1 .or. lat_index .gt. nlat_topo - 1) then
                  exit shooting_loop
                endif
                xgrid(1 : 2) = [lon_topo(lon_index), lon_topo(lon_index + 1)]
                ygrid(1 : 2) = [lat_topo(lat_index), lat_topo(lat_index + 1)]
                val_2d(1 : 2, 1 : 2) = topography(lon_index : lon_index + 1, lat_index : lat_index + 1)
                call linear_interpolation_2d(lon_tmp, lat_tmp, xgrid, ygrid, val_2d, topography_interpolate)
                if(depth_tmp .lt. topography_interpolate) then
                  !print '(a, 3(f9.4, 1x))', "ray surface arrived, lon/lat = ", lon_tmp, lat_tmp, depth_tmp
                  exit shooting_loop
                endif

                lon_index = int((lon_tmp - lon_str_w) / dlon_str) + 1
                lat_index = int((lat_tmp - lat_str_s) / dlat_str) + 1
                z_index   = int((depth_tmp - z_str_min) / dz_str) + 1

                !!exit if ray approaches to the boundary of the structure array
                if(lon_index .lt. 1 .or. lon_index .gt. nlon_str - 1      &
                &  .or. lat_index .lt. 1 .or. lat_index .gt. nlat_str - 1 &
                &  .or. z_index .lt. 1 .or. z_index .gt. nz_str - 1) then
                  exit shooting_loop
                endif

                xgrid(1) = lon_str_w + real(lon_index - 1, kind = fp) * dlon_str; xgrid(2) = xgrid(1) + dlon_str
                ygrid(1) = lat_str_s + real(lat_index - 1, kind = fp) * dlat_str; ygrid(2) = ygrid(1) + dlat_str
                zgrid(1) = z_str_min + real(z_index - 1, kind = fp) * dz_str;     zgrid(2) = zgrid(1) + dz_str

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
                  lon_min = lon_tmp
                  lat_min = lat_tmp
                  depth_min = depth_tmp
                  inc_angle_ini_min(kk) = inc_angle_ini
                  !print '(a, 4(f8.4, 1x))', "rayshoot_tmp lon, lat, depth, dist_min = ", lon_min, lat_min, depth_min, &
                  !&                                                                      dist_min
                endif
 
                !!shooting the ray
                val_1d(1 : 2) = velocity(lon_index, lat_index, z_index : z_index + 1, wavetype)
                call linear_interpolation_1d(depth_tmp, zgrid, val_1d, velocity_interpolate)
                val_3d(1 : 2, 1 : 2, 1 : 2) &
                &  = qinv(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1, wavetype)
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
          enddo incangle_loop2
          !print '(a, 4(f8.4, 1x))', "grid lon, lat, depth, az_ini = ", lon_grid, lat_grid, depth_grid, az_ini * rad2deg
          !print '(a, i0, 3a, 3(f8.4, 1x))', "jj = ", jj, " station ", &
          !&                                  trim(stname(jj)), " lon, lat, depth = ", lon_sta(jj), lat_sta(jj), z_sta(jj)
          !print '(a, 4(f8.4, 1x))', "rayshoot lon, lat, depth, inc_angle = ", lon_min, lat_min, depth_min, &
          !&                                                                   inc_angle_ini_min(nrayshoot) * rad2deg
          !print '(a, 3(f8.4, 1x))', "dist_min, ttime, width = ", dist_min, ttime_min(jj, i, j, k), width_min(jj, i, j, k)
          !print *, ""
          if(dist_min .gt. conv_rayshooting_threshold) then
            ttime_min(jj, i, j, k) = real(huge, kind = fp)
            width_min(jj, i, j, k) = 0.0_fp
          endif

#endif   /* -DRAYBENDING or not */
#endif   /* -DREAD_TTIME or not */

        enddo station_loop
      enddo lon_loop
      !$ write(0, '(3(a, i0), a)') "omp_thread_num = ", omp_thread, " lat index j = ", j, " depth index k = ", k, " end"
    enddo lat_loop
  enddo z_loop
  !$omp end do
  !$omp end parallel

#if defined (OUT_TTIMETABLE)
  open(unit = 10, file = "ttime_pulsewidth_table.dat", access = "direct", form = "unformatted", &
  &    status = "new", recl = fp, iostat = ios)
  if(ios .ne. 0) then
    close(10)
    write(0, '(a)') "Travel time table (ttime_pulsewidth_table.dat) has already existed."
    error stop
  else
    icount = 1
    do k = 1, nz
      do j = 1, nlat
        do i = 1, nlon
          do ii = 1, nsta
            write(10, rec = icount) ttime_min(ii, i, j, k)
            icount = icount + 1
            write(10, rec = icount) width_min(ii, i, j, k)
            icount = icount + 1
            write(10, rec = icount) hypodist(ii, i, j, k)
            icount = icount + 1
          enddo
        enddo
      enddo
    enddo
    close(10)
  endif
#endif

#endif /* -DREAD_TTIMETABLE */

  !!find minimum residual grid for seismic source 
  resultfile = trim(resultdir) // "/" // trim(resultfile)
  open(unit = 10, file = resultfile)
  write(10, '(a)') "# OT min_lon min_lat min_dep source_amp residual nsta_use"
  residual(1 : nlon, 1 : nlat, 1 : nz) = huge

#if defined (OUT_AMPLITUDE)
  ampfile = trim(resultdir) // "/" // "subevent_amplitude.txt"
  open(unit = 20, file = ampfile)
  write(20, '(a)', advance = "no") "# "
  do i = 1, nsta
    write(20, '(a, 1x)', advance = "no") stname(i)
  enddo
  write(20, *)
#endif

!!Do grid search
  time_count = 0
  time_loop: do
#if defined (AMP_TXT)
    if(time_count + 1 .gt. namp) exit time_loop
#else
    origintime = ot_begin + ot_shift * real(time_count, kind = fp)
    if(origintime .gt. ot_end) exit time_loop
    !write(0, '(a, f6.1)') "origin time = ", origintime
#endif

    !$omp parallel default(none), &
    !$omp&         shared(nsta, freq, ttime_min, origintime, sampling, npts, waveform_obs, source_amp, begin, nsta_use_grid, &
#if defined (AMP_TXT)
    !$omp&                amp_txt, time_count, &
#endif
    !$omp&                rms_amp_obs_noise, rms_tw, hypodist, ttime_cor, use_flag, siteamp, width_min, residual, &
    !$omp&                lon_topo, lat_topo, dlon_topo, dlat_topo, topography), &
    !$omp&         private(omp_thread, i, j, ii, jj, depth_grid, wave_index, rms_amp_obs, icount, residual_normalize, &
#if defined (AMP_RATIO)
    !$omp&                 amp_ratio_obs, amp_ratio_cal, &
#endif
    !$omp&                 lon_grid, lat_grid, lon_index, lat_index, xgrid, ygrid, val_2d, topography_interpolate, &
    !$omp&                 use_flag_tmp, rms_amp_cal)

    !$ omp_thread = omp_get_thread_num()


    !$omp do schedule(guided)
    z_loop2: do k = 1, nz
      depth_grid = z_min + dz * real(k - 1, kind = fp)
      !!$ write(0, '(2(a, i0))') "omp_thread_num = ", omp_thread, " depth index k = ", k
      lat_loop2: do j = 1, nlat
        lat_grid = lat_s + dlat * real(j - 1, kind = fp)
        lon_loop2: do i = 1, nlon
          lon_grid = lon_w + dlon * real(i - 1, kind = fp)

          source_amp(i, j, k) = 0.0_dp
          residual(i, j, k) = huge
          nsta_use_grid(i, j, k) = 0

          !!check the grid is lower than the topo
          lon_index = int((lon_grid - lon_topo(1)) / dlon_topo) + 1
          lat_index = int((lat_grid - lat_topo(1)) / dlat_topo) + 1
          xgrid(1 : 2) = [lon_topo(lon_index), lon_topo(lon_index + 1)]
          ygrid(1 : 2) = [lat_topo(lat_index), lat_topo(lat_index + 1)]
          val_2d(1 : 2, 1 : 2) = topography(lon_index : lon_index + 1, lat_index : lat_index + 1)
          call linear_interpolation_2d(lon_grid, lat_grid, xgrid, ygrid, val_2d, topography_interpolate)
          if(depth_grid .lt. topography_interpolate) cycle lon_loop2

          !!caluculate site-corrected amplitude
          do jj = 1, nsta
            rms_amp_obs(jj) = 0.0_fp

            if(ttime_min(jj, i, j, k) .eq. real(huge, kind = fp)) cycle
            if(use_flag(jj) .eqv. .false.) cycle

#if defined (AMP_TXT)
            rms_amp_obs(jj) = amp_txt(jj, time_count + 1)
#else
#if defined (WITHOUT_TTIME)
            wave_index = int((origintime - begin(jj)) / sampling(jj) + 0.5_fp) + 1
#else
            wave_index = int((origintime - begin(jj) + ttime_min(jj, i, j, k) + ttime_cor(jj, wavetype)) &
            &                / sampling(jj) + 0.5_fp) + 1
#endif
            icount = 0
            do ii = wave_index, wave_index + int(rms_tw / sampling(jj) + 0.5_fp) - 1
              if(ii .lt. npts(jj)) then
                rms_amp_obs(jj) = rms_amp_obs(jj) + waveform_obs(ii, jj) * waveform_obs(ii, jj)
                icount = icount + 1
              endif
            enddo
            if(icount .ne. 0) rms_amp_obs(jj) = sqrt(rms_amp_obs(jj) / real(icount, kind = dp))
#endif /* -DAMP_TXT */

            !!calculate source amplitude temporally from high-s/n ratio staitons
            if(rms_amp_obs(jj) .gt. snratio_accept * rms_amp_obs_noise(jj)) then
              nsta_use_grid(i, j, k) = nsta_use_grid(i, j, k) + 1
              source_amp(i, j, k) = source_amp(i, j, k) &
              &            + rms_amp_obs(jj) / siteamp(jj) &
              &            * real(hypodist(jj, i, j, k) * exp(width_min(jj, i, j, k) * (pi * freq)), kind = dp)
              !print *, source_amp(i, j, k)
            endif
          enddo
          if(nsta_use_grid(i, j, k) .ne. 0) source_amp(i, j, k) = source_amp(i, j, k) / real(nsta_use_grid(i, j, k), kind = dp)

          !!check s/n ratio
          nsta_use_grid(i, j, k) = 0
          do jj = 1, nsta
            if(ttime_min(jj, i, j, k) .eq. real(huge, kind = fp)) cycle
            !!check whether expected amplitude is large or not (Doi et al., 2020, SSJ meeting)
            !rms_amp_cal = source_amp(i, j, k) * siteamp(jj) &
            !&             / real(hypodist(jj, i, j, k), kind = dp) * real(exp(-pi * freq * width_min(jj, i, j, k)), kind = dp)
            !if(use_flag(jj) .eqv. .false.) then
            !  use_flag_tmp(jj) = .false.
            !else
            !  if(rms_amp_cal .gt. snratio_accept * rms_amp_obs_noise(jj) .or. &
            !  &  rms_amp_obs(jj) .gt. snratio_accept * rms_amp_obs_noise(jj)) then
            !    use_flag_tmp(jj) = .true.
            !    nsta_use_grid(i, j, k) = nsta_use_grid(i, j, k) + 1
            !  else
            !    use_flag_tmp(jj) = .false.
            !  endif
            !endif
            if(use_flag(jj) .eqv. .false.) then
              use_flag_tmp(jj) = .false.
            else
              !!just compare amplitude  of signal and noise
              if(rms_amp_obs(jj) .gt. snratio_accept * rms_amp_obs_noise(jj)) then
                use_flag_tmp(jj) = .true.
                nsta_use_grid(i, j, k) = nsta_use_grid(i, j, k) + 1
              else
                use_flag_tmp(jj) = .false.
              endif
            endif
          enddo
          if(nsta_use_grid(i, j, k) .lt. nsta_use_minimum) cycle lon_loop2

          !!calculate source amplitude again
          source_amp(i, j, k) = 0.0_dp
          nsta_use_grid(i, j, k) = 0
          do jj = 1, nsta
            if(use_flag_tmp(jj) .eqv. .false.) cycle
            if(ttime_min(jj, i, j, k) .eq. real(huge, kind = fp)) cycle

            nsta_use_grid(i, j, k) = nsta_use_grid(i, j, k) + 1
            source_amp(i, j, k) = source_amp(i, j, k) &
            &            + rms_amp_obs(jj) / siteamp(jj) &
            &            * real(hypodist(jj, i, j, k) * exp(width_min(jj, i, j, k) * (pi * freq)), kind = dp)
          enddo
          if(nsta_use_grid(i, j, k) .ne. 0) source_amp(i, j, k) = source_amp(i, j, k) / real(nsta_use_grid(i, j, k), kind = dp)
           

#if defined (AMP_RATIO)
          !!calculate amplitude ration and residual
          residual(i, j, k) = 0.0_dp
          nsta_use_grid(i, j, k) = 0
          do jj = 1, nsta - 1
            if(use_flag_tmp(jj) .eqv. .false.) cycle
            if(ttime_min(jj, i, j, k) .eq. real(huge, kind = fp)) cycle
            do ii = jj + 1, nsta
              if(use_flag_tmp(ii) .eqv. .false.) cycle
              if(ttime_min(ii, i, j, k) .eq. real(huge, kind = fp)) cycle
              nsta_use_grid(i, j, k) = nsta_use_grid(i, j, k) + 1
              amp_ratio_obs = rms_amp_obs(ii) / rms_amp_obs(jj)
              amp_ratio_cal = (siteamp(ii) / siteamp(jj)) &
              &  * (real(exp(-pi * freq * width_min(ii, i, j, k)), kind = dp) / real(hypodist(ii, i, j, k), kind = dp)) &
              &  / (real(exp(-pi * freq * width_min(jj, i, j, k)), kind = dp) / real(hypodist(jj, i, j, k), kind = dp))
              residual(i, j, k) = residual(i, j, k) + ((amp_ratio_obs - amp_ratio_cal) / amp_ratio_cal) ** 2
            enddo
          enddo

          if(nsta_use_grid(i, j, k) .ne. 0) then
            residual(i, j, k) = sqrt(2.0_dp / real(nsta_use_grid(i, j, k), kind = dp) * residual(i, j, k))
          else
            residual(i, j, k) = huge
          endif
#else
          residual(i, j, k) = 0.0_dp
          residual_normalize = 0.0_dp
          nsta_use_grid(i, j, k) = 0
          do ii = 1, nsta
            if(use_flag_tmp(ii) .eqv. .false.) cycle
            if(ttime_min(ii, i, j, k) .eq. real(huge, kind = fp)) cycle
            nsta_use_grid(i, j, k) = nsta_use_grid(i, j, k) + 1
            residual(i, j, k) = residual(i, j, k) &
            &                 + (rms_amp_obs(ii) / siteamp(ii) &
            &                 - source_amp(i, j, k) / real(hypodist(ii, i, j, k), kind = dp) &
            &                   * real(exp(-pi * freq * width_min(ii, i, j, k)), kind = dp)) ** 2
            residual_normalize = residual_normalize + (rms_amp_obs(ii) / siteamp(ii)) ** 2
            !residual_normalize = source_amp(i, j, k) / real(hypodist(ii, i, j, k), kind = dp) &
            !&                    * real(exp(-pi * freq * width_min(ii, i, j, k)), kind = dp)
            !residual(i, j, k) = residual(i, j, k) &
            !&                 + ((rms_amp_obs(ii) / siteamp(ii) - residual_normalize) / residual_normalize) ** 2
          enddo
          residual(i, j, k) = residual(i, j, k) / residual_normalize
          !residual(i, j, k) = residual(i, j, k) / real(nsta_use_grid(i, j, k), kind = dp)
#endif

        enddo lon_loop2
      enddo lat_loop2
    enddo z_loop2
    !$omp end do
    !$omp end parallel

    !!search minimum residual
    residual_minloc = minloc(residual)
    lon_grid = lon_w + real(residual_minloc(1) - 1, kind = fp) * dlon
    lat_grid = lat_s + real(residual_minloc(2) - 1, kind = fp) * dlat
    depth_grid = z_min + real(residual_minloc(3) - 1, kind = fp) * dz

#if defined (AMP_TXT)
    write(0, '(2a, a, f0.4, 1x, f0.4, 1x, f0.2, a, 2(a, e15.7), a, i0)') &
    &                 "Index = ", trim(eventindex(time_count + 1)), " residual_minimum (lon, lat, dep) = (", &
    &                             lon_grid, lat_grid, depth_grid, ")", &
    &                 " source_amp = ", source_amp(residual_minloc(1), residual_minloc(2), residual_minloc(3)), &
    &                 " residual = ", residual(residual_minloc(1), residual_minloc(2), residual_minloc(3)), &
    &                 " nsta_use = ", nsta_use_grid(residual_minloc(1), residual_minloc(2), residual_minloc(3))
    
    write(10, '(a, 1x, f0.4, 1x, f0.4, 1x, f0.2, 2(1x, e15.7), 1x, i0)') &
    &                 trim(eventindex(time_count + 1)), lon_grid, lat_grid, depth_grid, &
    &                 source_amp(residual_minloc(1), residual_minloc(2), residual_minloc(3)), &
    &                 residual(residual_minloc(1), residual_minloc(2), residual_minloc(3)), &
    &                 nsta_use_grid(residual_minloc(1), residual_minloc(2), residual_minloc(3))
#else
    write(0, '(a, f0.1, a, f0.4, 1x, f0.4, 1x, f0.2, a, 2(a, e15.7), a, i0)') &
    &                 "OT = ", origintime, " residual_minimum (lon, lat, dep) = (", lon_grid, lat_grid, depth_grid, ")", &
    &                 " source_amp = ", source_amp(residual_minloc(1), residual_minloc(2), residual_minloc(3)), &
    &                 " residual = ", residual(residual_minloc(1), residual_minloc(2), residual_minloc(3)), &
    &                 " nsta_use = ", nsta_use_grid(residual_minloc(1), residual_minloc(2), residual_minloc(3))
    write(10, '(f0.1, 1x, f0.4, 1x, f0.4, 1x, f0.2, 2(1x, e15.7), 1x, i0)') &
    &                 origintime, lon_grid, lat_grid, depth_grid, &
    &                 source_amp(residual_minloc(1), residual_minloc(2), residual_minloc(3)), &
    &                 residual(residual_minloc(1), residual_minloc(2), residual_minloc(3)), &
    &                 nsta_use_grid(residual_minloc(1), residual_minloc(2), residual_minloc(3))
#endif

    !!output grd files
    call int_to_char(time_count, maxlen, time_count_char)
    !!output horizontal slice
    grdfile = trim(resultdir) // "/min_err_lon-lat_" // trim(time_count_char) // ".grd"
    allocate(residual_grd(nlon, nlat))
    residual_grd(1 : nlon, 1 : nlat) = real(residual(1 : nlon, 1 : nlat, residual_minloc(3)), kind = sp)
    call write_grdfile_2d(lon_w, lat_s, dlon, dlat, nlon, nlat, residual_grd, grdfile, nanval = real(huge, kind = sp))
    deallocate(residual_grd)

    !!output depth slice -- lon-dep
    grdfile = trim(resultdir) // "/min_err_lon-dep_" // trim(time_count_char) // ".grd"
    allocate(residual_grd(nlon, nz))
    residual_grd(1 : nlon, 1 : nz) = real(residual(1 : nlon, residual_minloc(2), 1 : nz), kind = sp)
    call write_grdfile_2d(lon_w, z_min, dlon, dz, nlon, nz, residual_grd, grdfile, nanval = real(huge, kind = sp))
    deallocate(residual_grd)

    !!output depth slice -- dep-lat
    grdfile = trim(resultdir) // "/min_err_dep-lat_" // trim(time_count_char) // ".grd"
    allocate(residual_grd(nz, nlat))
    do j = 1, nz
      do i = 1, nlat
        residual_grd(j, i) = real(residual(residual_minloc(1), i, j), kind = sp)
      enddo
    enddo
    call write_grdfile_2d(z_min, lat_s, dz, dlat, nz, nlat, residual_grd, grdfile, nanval = real(huge, kind = sp))
    deallocate(residual_grd)

#if defined (OUT_AMPLITUDE)
    do i = 1, nsta
#if defined (WITHOUT_TTIME)
      wave_index = int((origintime - begin(i)) / sampling(i) + 0.5_fp) + 1
#else
      wave_index = int((origintime - begin(i) &
      &          + ttime_min(i, residual_minloc(1), residual_minloc(2), residual_minloc(3)) + ttime_cor(i, wavetype)) &
      &          / sampling(i) + 0.5_fp) + 1
#endif

#if defined (WIN) || defined (SAC)
      rms_amp_obs(i) = 0.0_fp
      icount = 0
      do ii = wave_index, wave_index + int(rms_tw / sampling(i) + 0.5_fp) - 1
         if(ii .lt. npts(i)) then
           rms_amp_obs(i) = rms_amp_obs(i) + waveform_obs(ii, i) * waveform_obs(ii, i)
           icount = icount + 1
         endif
      enddo
      rms_amp_obs(i) = sqrt(rms_amp_obs(i) / real(icount, kind = dp))
      if(rms_amp_obs(i) .le. snratio_accept * rms_amp_obs_noise(i)) rms_amp_obs(i) = 0.0_fp
#elif defined (AMP_TXT)
      rms_amp_obs(i) = amp_txt(i, time_count + 1)
#endif

      write(20, '(e15.7, 1x)', advance = "no") rms_amp_obs(i)
    enddo

#if defined (WIN) || defined (SAC)
    write(20, '(f0.1)') origintime
#else
    write(20, '(a)') trim(eventindex(time_count + 1))
#endif

#endif /* -DOUT_AMPLITUDE */
    

    time_count = time_count + 1
  enddo time_loop
  close(10)

#if defined (OUT_AMPLITUDE)
  close(20)
#endif

  stop
end program AmplitudeSourceLocation_PulseWidth

