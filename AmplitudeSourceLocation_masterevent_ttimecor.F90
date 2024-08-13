program AmplitudeSourceLocation_masterevent
#define RAY_BENDING
  !!Relative amplitude source location (master event)
  !!using depth-dependent 1D velocity structure, 3D heterogeneous attenuation structure
  !!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2020
  !!License  : MIT License (https://opensource.org/licenses/MIT)

  !!-DWIN: use win-formatted waveform file for input
  !!-DSAC: use sac-formatted(binary) waveform file  for input
  !!default: txt format
  !!-DWITHOUT_ERROR: do not calculate estimation errors

  !!calculate traveltime correction table for determining events which
  !interevent distance are long

  use nrtype,               only : fp, dp
  use constants,            only : rad2deg, deg2rad, pi, r_earth
  use rayshooting,          only : rayshooting3D
  use set_velocity_model,   only : set_velocity
  use linear_interpolation, only : linear_interpolation_1d, linear_interpolation_2d, linear_interpolation_3d, &
  &                                block_interpolation_3d
  use greatcircle,          only : greatcircle_dist
  use grdfile_io,           only : read_grdfile_2d
  use raybending,           only : pseudobending3D

#if defined (WIN)
  use m_win
  use m_winch
#endif

#if defined (SAC)
  use read_sacfile,         only : read_sachdr, read_sacdata
#endif

#if defined (MKL)
  use lapack95
#else
  use f95_lapack
#endif

  implicit none

  integer,            parameter :: wavetype = 2          !!1 for P-wave, 2 for S-wave
  integer,            parameter :: nsta_use_minimum = 5
  real(kind = fp),    parameter :: snratio_accept = 0.0_fp
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

  !!traveltime correction table
  real(kind = fp),    parameter :: lon_cor_w = 136.5_fp, lon_cor_e = 137.5_fp
  real(kind = fp),    parameter :: lat_cor_s = 37.2_fp,  lat_cor_n = 37.4_fp
  real(kind = fp),    parameter :: z_cor_min = 0.0_fp,   z_cor_max = 0.0_fp
  real(kind = fp),    parameter :: dlon_cor  = 0.02_fp,   dlat_cor = 0.02_fp, dz_cor = 1.0_fp
  integer,            parameter :: nlon_cor = int((lon_cor_e - lon_cor_w) / dlon_cor) + 1
  integer,            parameter :: nlat_cor = int((lat_cor_n - lat_cor_s) / dlat_cor) + 1
  integer,            parameter :: nz_cor   = int((z_cor_max - z_cor_min) / dz_cor) + 1
  real(kind = dp)               :: residual_loc_cor(1 : nlon_cor, 1 : nlat_cor, 1 : nz_cor)

#if defined (WIN) /* use win-format waveform file for input waveforms */
  real(kind = dp),    parameter   :: order = 1.0e+6_dp
  real(kind = fp),    allocatable :: sampling(:), begin(:), ttime(:), ttime_loc_cor(:, :, :, :)
  real(kind = dp),    allocatable :: waveform_obs(:, :)
  integer,            allocatable :: ikey(:), sampling_int(:), waveform_obs_int(:, :), npts_win(:, :), npts(:)
  character(len = 4), allocatable :: st_winch(:)
  type(winch__hdr),   allocatable :: chtbl(:)
  integer                         :: nsec, tim, time_index
  character(len = 129)            :: win_filename, win_chfilename
  character(len = 10)             :: cmpnm
#endif

#if defined (SAC) /* use sac-format wavefome files (NVHDR=6) */
  !real(kind = dp),    parameter   :: order = 1.0e+6_dp
  real(kind = dp),    parameter   :: order = 1.0_dp
  character(len = 3), parameter   :: sacfile_extension = "sac"
  real(kind = dp),    allocatable :: waveform_obs(:, :)
  real(kind = fp),    allocatable :: begin(:), sampling(:), stime(:), ttime(:), ttime_loc_cor(:, :, :, :)
  integer,            allocatable :: npts(:)
  integer                         :: time_index
  character(len = 129)            :: sacfile, sacfile_index
  character(len = 10)             :: cmpnm
#endif

  real(kind = fp),    allocatable :: stlon(:), stlat(:), stdp(:), ttime_cor(:, :)
  character(len = 6), allocatable :: stname(:)
  logical,            allocatable :: use_flag(:), use_flag_tmp(:)
  character(len = 10)             :: freq_t

#if defined (WIN) || defined (SAC)
  real(kind = fp)                 :: ot_begin, ot_end, ot_shift, ot_tmp, rms_tw, amp_avg
  character(len = 10)             :: ot_begin_t, ot_end_t, ot_shift_t, rms_tw_t
#endif

  !!bandpass filter
  real(kind = dp),    parameter   :: ap = 0.5_dp, as = 5.0_dp
  real(kind = dp),    allocatable :: h(:)
  character(len = 6)              :: fl_t, fh_t, fs_t
  real(kind = dp)                 :: fl, fh, fs, gn, c
  integer                         :: m, n

#if defined (DAMPED) /* constraints */
  real(kind = fp),    parameter :: damp(4) = [0.0_fp, 0.0_fp, 0.0_fp, 1.0_fp]  !![amp, dx, dy, dz]
#endif

  real(kind = fp),    parameter :: alt_to_depth = -1.0e-3_fp
  real(kind = dp),    parameter :: huge = 1.0e+6_dp

  real(kind = fp)               :: velocity(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                qinv(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                val_1d(1 : 2), val_2d(1 : 2, 1 : 2), val_3d(1 : 2, 1 : 2, 1 : 2), &
  &                                xgrid(1 : 2), ygrid(1 : 2), zgrid(1 : 2), inc_angle_ini_min(0 : nrayshoot), &
  &                                normal_vector(1 : 3) 
  integer                       :: min_residual(3)
  real(kind = dp),  allocatable :: topography(:, :), lon_topo(:), lat_topo(:)
  real(kind = fp),  allocatable :: obsamp_master(:), obsamp_sub(:, :), obsamp_noise(:), hypodist(:), ray_azinc(:, :), &
  &                                dist_min(:), obsvector(:, :), obsvector_copy(:), data_residual(:), &
  &                                inversion_matrix(:, :), inversion_matrix_copy(:, :), inversion_matrix_sub(:, :), &
  &                                sigma_inv_data(:, :), error_matrix(:, :), error_matrix_sub(:, :), obsvector_tmp(:, :, :, :)
  integer,          allocatable :: ipiv(:), nsta_use_tmp(:)
  real(kind = fp)               :: evlon_master, evlat_master, evdp_master, &
  &                                epdist, epdelta, az_ini, dinc_angle_org, dinc_angle, inc_angle_ini, &
  &                                lon_tmp, lat_tmp, depth_tmp, az_tmp, inc_angle_tmp, dist_tmp, velocity_interpolate, &
  &                                dvdz, lon_new, lat_new, depth_new, az_new, inc_angle_new, ttime_tmp, matrix_const, &
  &                                lon_min, lat_min, depth_min, delta_depth, delta_lon, delta_lat, &
  &                                data_variance, sigma_lon, sigma_lat, sigma_depth, sigma_amp, &
  &                                depth_max_tmp, depth_max, siteamp_tmp
  real(kind = dp)               :: topography_interpolate, dlon_topo, dlat_topo, qinv_interpolate, freq
  integer                       :: nlon_topo, nlat_topo, nsta, nsubevent, lon_index, lat_index, z_index, &
  &                                i, j, k, ii, jj, kk, icount, ncount, nsta_use, ios
  character(len = 129)          :: topo_grd, station_param, masterevent_param, subevent_param, resultfile
  character(len = 20)           :: cfmt, nsta_c

#if defined (RAY_BENDING)
  integer,                parameter :: ndiv_raypath = 10
  integer,                parameter :: nraypath_ini = 4
  real(kind = fp)                   :: raypath_lon((nraypath_ini - 1) * 2 ** ndiv_raypath + 1), &
  &                                    raypath_lat((nraypath_ini - 1) * 2 ** ndiv_raypath + 1), &
  &                                    raypath_dep((nraypath_ini - 1) * 2 ** ndiv_raypath + 1)
  integer                           :: nraypath
#endif

  !!OpenMP variable
  !$ integer                    :: omp_thread

  icount = iargc()

#if defined (WIN)
  if(icount .ne. 14) then
    write(0, '(a)', advance="no") "usage: ./asl_masterevent "
    write(0, '(a)', advance="no") "(topography_grd) (station_param_file) (masterevent_param_file) (win_waveform) (win_chfile) "
    write(0, '(a)', advance="no") "(component name) (fl) (fh) (fs) (ot_begin) (ot_end) (ot_shift) (rms_time_window_length) "
    write(0, '(a)')               "(result_file)"
    error stop
  endif
  call getarg(1, topo_grd)
  call getarg(2, station_param)
  call getarg(3, masterevent_param)
  call getarg(4, win_filename)
  call getarg(5, win_chfilename)
  call getarg(6, cmpnm)
  call getarg(7, fl_t)       ; read(fl_t, *) fl
  call getarg(8, fh_t)       ; read(fh_t, *) fh
  call getarg(9, fs_t)       ; read(fs_t, *) fs
  call getarg(10, ot_begin_t); read(ot_begin_t, *) ot_begin
  call getarg(11, ot_end_t)  ; read(ot_end_t, *) ot_end
  call getarg(12, ot_shift_t); read(ot_shift_t, *) ot_shift
  call getarg(13, rms_tw_t)  ; read(rms_tw_t, *) rms_tw
  call getarg(14, resultfile)
#elif defined (SAC)
  if(icount .ne. 13) then
    write(0, '(a)', advance="no") "usage: ./asl_masterevent "
    write(0, '(a)', advance="no") "(topography_grd) (station_param_file) (masterevent_param_file) (sacfile_index) "
    write(0, '(a)', advance="no") "(component_name) (fl) (fh) (fs) (ot_begin) (ot_end) (ot_shift) (rms_time_window_length) "
    write(0, '(a)')               "(result_file)"
    error stop
  endif
  call getarg(1, topo_grd)
  call getarg(2, station_param)
  call getarg(3, masterevent_param)
  call getarg(4, sacfile_index)
  call getarg(5, cmpnm)
  call getarg(6, fl_t)       ; read(fl_t, *) fl
  call getarg(7, fh_t)       ; read(fh_t, *) fh
  call getarg(8, fs_t)       ; read(fs_t, *) fs
  call getarg(9, ot_begin_t) ; read(ot_begin_t, *) ot_begin
  call getarg(10, ot_end_t)  ; read(ot_end_t, *) ot_end
  call getarg(11, ot_shift_t); read(ot_shift_t, *) ot_shift
  call getarg(12, rms_tw_t)  ; read(rms_tw_t, *) rms_tw
  call getarg(13, resultfile)
#else
  if(icount .ne. 6) then
    write(0, '(a)', advance="no") "usage: ./asl_masterevent "
    write(0, '(a)', advance="no") "(topography_grd) (station_param_file) (masterevent_param_file) (subevent_param_file) "
    write(0, '(a)')               "(frequency) (result_file)"
    error stop
  endif
  call getarg(1, topo_grd)
  call getarg(2, station_param)
  call getarg(3, masterevent_param)
  call getarg(4, subevent_param)
  call getarg(5, freq_t); read(freq_t, *) freq
  call getarg(6, resultfile)
#endif

#if defined (DAMPED)
  write(0, '(a)') "-DDAMPED set"
  write(0, '(a, 4(f5.2, 1x), a)') "Damped matrix = [ ", (damp(i), i = 1, 4), "]"
#endif

#if defined (WIN) || defined (SAC)
  write(0, '(a, 3(f5.2, 1x))') "Bandpass filter parameter fl, fh, fs (Hz) = ", fl, fh, fs
  freq = (fl + fh) * 0.5_dp
#endif

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
  if(nsta .lt. 4) then
    close(10)
    write(0, '(a)') "Number of stations, nsta, should be larger than 4"
    error stop
  endif
  rewind(10)
  allocate(stlon(1 : nsta), stlat(1 : nsta), stdp(1 : nsta), stname(1 : nsta), ttime_cor(1 : nsta, 1 : 2), &
  &        use_flag(1 : nsta), use_flag_tmp(1 : nsta), obsamp_noise(1 : nsta))
  nsta_use = 0
  do i = 1, nsta
    read(10, *) stlon(i), stlat(i), stdp(i), stname(i), use_flag(i), ttime_cor(i, 1), ttime_cor(i, 2), siteamp_tmp, &
    &           obsamp_noise(i)
    write(0, '(a, i0, a, f9.4, a, f8.4, a, f6.3, 1x, a7, l2)') &
    &     "station(", i, ") lon(deg) = ", stlon(i), " lat(deg) = ", stlat(i), " depth(km) = ", stdp(i), trim(stname(i)), &
    &     use_flag(i)
    if(use_flag(i) .eqv. .true.) nsta_use = nsta_use + 1
  enddo
  close(10)
  if(nsta_use .lt. 4) then
    close(10)
    write(0, '(a)') "Number of stations for calculation, nsta_use, should be larger than 4"
    error stop
  endif


  !!read masterevent parameter
  open(unit = 10, file = masterevent_param)
  nev_master = 0
  read(10. *)
  do
    read(10, *, iostat = ios)
    if(ios .ne. 0) exit
    nev_master = nev_master + 1
  enddo
  nev_master = nev_master / 2
  allocate(evlon_master(1 : nev_master), evlat_master(1 : nev_master), evdp_master(1 : nev_master))
  allocate(obsamp_master(1 : nsta, 1 : nev_master))
  rewind(10)
  read(10, *)
  do j = 1, nev_master
    read(10, *) evlon_master(j), evlat_master(j), evdp_master(j)
    read(10, *) (obsamp_master(i, j), i = 1, nsta)
  enddo
  close(10)
  nnearest_evmaster = min(2, nev_master)
  allocate(nearest_evmaster_index(1 : nnearest_evmaster, 1 : nlon_cor, 1 : nlat_cor, 1 : nz_cor))
  allocate(evmaster_dist(1 : nnearest_evmaster))

#if defined (WIN) || defined (SAC)
  allocate(sampling(1 : nsta), npts(1 : nsta), begin(1 : nsta), ttime(1 : nsta))
  allocate(ttime_loc_cor(1 : nsta, 1 : nlon_cor, 1 : nlat_cor, 1 : nz_cor))

#if defined (WIN)
  allocate(st_winch(nsta), sampling_int(1 : nsta))
  !!read channel table
  call winch__read_tbl(trim(win_chfilename), chtbl)
  !!find chid from given stname(nsta) and cmpnm
  allocate(ikey(nsta))
  do i = 1, nsta
    call winch__st2chid(chtbl, stname(i), trim(cmpnm), st_winch(i), ikey(i))
    if(ikey(i) .eq. 0) then
      write(0, '(5(a, 1x))') "station/comp =", trim(stname(i)), trim(cmpnm), "does not exist in", trim(win_chfilename)
      error stop
    endif
    write(0, '(6a)') "station name = ", trim(stname(i)), " comp = ", trim(cmpnm), " chid = ", st_winch(i)
  enddo
  !!read all waveform data from win-formatted file
  call win__read_file(trim(win_filename), st_winch, sampling_int, nsec, tim, waveform_obs_int, npts_win)
  sampling(1 : nsta) = 1.0_fp / real(sampling_int(1 : nsta), kind = fp)
  npts(1 : nsta) = sampling_int(1 : nsta) * nsec
  begin(1 : nsta) = 0.0_fp
  !!convert waveform data from digitized value to physical value, remove offset
  allocate(waveform_obs(1 : maxval(npts), 1 : nsta))
  waveform_obs(1 : maxval(npts), 1 : nsta) = 0.0_dp
  do j = 1, nsta
    do i = 1, npts(j)
      waveform_obs(i, j) = real(waveform_obs_int(i, j), kind = dp) * chtbl(ikey(j))%conv * order
    enddo
  enddo
#elif defined (SAC)
  !!read waveform data from sac
  allocate(stime(1 : nsta))
  do i = 1, nsta
    !sacfile = trim(sacfile_index) // "." // trim(stname(i)) // "." // trim(cmpnm) // "." // sacfile_extension
    sacfile = trim(sacfile_index) // trim(stname(i)) // "_env_cal." // sacfile_extension
    call read_sachdr(sacfile, delta=sampling(i), npts=npts(i), begin=begin(i), t0=stime(i))
  enddo
  allocate(waveform_obs(maxval(npts), nsta))
  do i = 1, nsta
    !sacfile = trim(sacfile_index) // "." // trim(stname(i)) // "." // trim(cmpnm) // "." // sacfile_extension
    sacfile = trim(sacfile_index) // trim(stname(i)) // "_env_cal." // sacfile_extension
    call read_sacdata(sacfile, npts(i), waveform_obs(:, i))
    waveform_obs(1 : npts(i), i) = waveform_obs(1 : npts(i), i) * order
  enddo
#endif

#if defined (TESTDATA)
#else
  !!remove offset
  do j = 1, nsta
    amp_avg = 0.0_dp
    icount = 0
    do i = 1, int(rms_tw / sampling(j))
      amp_avg = amp_avg + waveform_obs(i, j)
      icount = icount + 1
    enddo
    amp_avg = amp_avg / real(icount, kind = dp)
    waveform_obs(1 : npts(j), j) = waveform_obs(1 : npts(j), j) - amp_avg
  enddo
  !!bandpass filter
  write(0, '(a)') "Applying bandpass filter"
  do j = 1, nsta
    call calc_bpf_order(fl, fh, fs, ap, as, real(sampling(j), kind = dp), m, n, c)
    allocate(h(4 * m))
    call calc_bpf_coef(fl, fh, real(sampling(j), kind = dp), m, n, h, c, gn)
    call tandem1(waveform_obs(:, j), waveform_obs(:, j), npts(j), h, m, 1)
    waveform_obs(1 : npts(j), j) = waveform_obs(1 : npts(j), j) * gn
    deallocate(h)
  enddo
  allocate(obsamp_noise(1 : nsta))
  do j = 1, nsta
    obsamp_noise(j) = 0.0_fp
    icount = 0
    do i = int(rms_tw / sampling(j)) + 1, int(rms_tw / sampling(j)) * 2 
      obsamp_noise(j) = obsamp_noise(j) + waveform_obs(i, j) ** 2
    enddo
    obsamp_noise(j) = sqrt(obsamp_noise(j) / real(icount, kind = fp))
  enddo
#endif

#else /* !(-DWIN) && !(-DSAC) */
  !!read amplitudes of subevents
  open(unit = 10, file = subevent_param)
  read(10, *)
  nsubevent = 0
  do 
    read(10, *, iostat = ios)
    if(ios .ne. 0) exit
    nsubevent = nsubevent + 1
  enddo
  write(0, '(a, i0)') "nsubevent = ", nsubevent
  rewind(10)
  read(10, *)
  allocate(obsamp_sub(1 : nsta, 1 : nsubevent), nsta_use_tmp(1 : nsubevent))
  do j = 1, nsubevent
    read(10, *) (obsamp_sub(i, j), i = 1, nsta)
  enddo
  close(10)
  allocate(obsamp_noise(1 : nsta))
  obsamp_noise(1 : nsta) = 0.0_fp
#endif /* -DWIN || -DSAC */

  !!set velocity/attenuation structure
  call set_velocity(z_str_min, dz_str, velocity, qinv)

  !!calculate ray length, pulse width, unit vector of ray incident
  write(0, '(a)') "calculate ray length, pulse width, and ray incident vector for master event"

  !!check whether the depth of master event location is lower than the topo
  do i = 1, nev_master
    lon_index = int((evlon_master(i) - lon_topo(1)) / dlon_topo) + 1
    lat_index = int((evlat_master(i) - lat_topo(1)) / dlat_topo) + 1
    xgrid(1 : 2) = [lon_topo(lon_index), lon_topo(lon_index + 1)]
    ygrid(1 : 2) = [lat_topo(lat_index), lat_topo(lat_index + 1)]
    val_2d(1 : 2, 1 : 2) = topography(lon_index : lon_index + 1, lat_index : lat_index + 1)
    call linear_interpolation_2d(evlon_master(i), evlat_master(i), xgrid, ygrid, val_2d, topography_interpolate)
    if(evdp_master(i) .lt. topography_interpolate) then
      write(0, '(a, f5.2, a)') "Depth of master event ", evdp_master(i), " is higher than the altitude there."
      error stop
    endif
  enddo

  allocate(hypodist(1 : nsta), ray_azinc(1 : 2, 1 : nsta), dist_min(1 : nsta))
  
  !!calculate traveltime to station        
  masterevent_loop: do kk = 1, nev_master
    station_loop: do jj = 1, nsta
      !!calculate azimuth and hypocentral distance
      call greatcircle_dist(evlat_master(kk), evlon_master(kk), stlat(jj), stlon(jj), &
      &                     distance = epdist, azimuth = az_ini, delta_out = epdelta)
      lon_index = int((stlon(jj) - lon_str_w) / dlon_str) + 1
      lat_index = int((stlat(jj) - lat_str_s) / dlat_str) + 1
      z_index   = int((stdp(jj) - z_str_min) / dz_str) + 1
      !print *, lon_sta(jj), lon_w + real(lon_index - 1) * dlon
      !print *, lat_sta(jj), lat_s + real(lat_index - 1) * dlat
      hypodist(jj, kk) = sqrt((r_earth - evdp_master(kk)) ** 2 + (r_earth - stdp(jj)) ** 2 &
      &                       - 2.0_fp * (r_earth - evdp_master(kk)) * (r_earth - stdp(jj)) * cos(epdelta))
      ray_azinc(1, jj, kk) = az_ini

#if defined (V_CONST)
      !!homogeneous structure: ray incident angle is calculated using cosine function (assuming cartesian coordinate)
      ray_azinc(2, jj, kk) = acos(epdist / hypodist(jj, kk)) + pi / 2.0_fp
#else
      !!do ray tracing with pseudobending scheme
      raypath_lon(1) = evlon_master(kk)
      raypath_lat(1) = evlat_master(kk)
      raypath_dep(1) = evdp_master(kk)
      raypath_lon(nraypath_ini) = stlon(jj)
      raypath_lat(nraypath_ini) = stlat(jj)
      raypath_dep(nraypath_ini) = stdp(jj)
      do i = 2, nraypath_ini - 1
        raypath_lon(i) = raypath_lon(i - 1) + (raypath_lon(nraypath_ini) - raypath_lon(1)) / real(nraypath_ini, kind = fp)
        raypath_lat(i) = raypath_lat(i - 1) + (raypath_lat(nraypath_ini) - raypath_lat(1)) / real(nraypath_ini, kind = fp)
        raypath_dep(i) = raypath_dep(i - 1) + (raypath_dep(nraypath_ini) - raypath_dep(1)) / real(nraypath_ini, kind = fp)
      enddo
      nraypath = nraypath_ini
      call pseudobending3D(raypath_lon, raypath_lat, raypath_dep, nraypath, ndiv_raypath, &
      &                    velocity(:, :, :, wavetype), lon_str_w, lat_str_s, z_str_min, dlon_str, dlat_str, dz_str, &
      &                    ttime_tmp, ray_az = ray_azinc(1, jj, kk), ray_incangle = ray_azinc(2, jj, kk))

#if defined (WIN) || defined (SAC)
      ttime(jj, kk) = ttime_tmp + ttime_cor(jj, wavetype)
#endif

      write(0, '(a, 4(f8.4, 1x))') "masterevent lon, lat, depth, az_ini = ", &
      &                            evlon_master(kk), evlat_master(kk), evdp_master(kk), az_ini * rad2deg
      write(0, '(a, 3(f8.4, 1x))') "station lon, lat, depth = ", stlon(jj), stlat(jj), stdp(jj)
      write(0, '(a, 2(f8.4, 1x))') "ray azimuth and inc_angle (deg) = ", &
      &                             ray_azinc(1, jj, kk) * rad2deg, ray_azinc(2, jj, kk) * rad2deg
#if defined (WIN) || defined (SAC)
      write(0, '(a, f5.2)') "traveltime (s) = ", ttime(jj, kk)
#endif

#endif /* -DV_CONST or not */

  enddo station_loop
!  !$omp end do


  !!make traveltime correction table
  do k = 1, nz_cor
    depth_tmp = z_cor_min + dz_cor * real(k, kind = fp)
    do j = 1, nlat_cor
      lat_tmp = lat_cor_s + dlat_cor * real(j, kind = fp)
      do i = 1, nlon_cor
        lon_tmp = lon_cor_w + dlon_cor * real(i, kind = fp)
        evmaster_dist(1 : nnearest_evmaster) = huge
        nearest_evmaster_index(1 : nnearest_evmaster, i, j, k) = 1
        do jj = 1, nev_master
          call greatcircle_dist(evlat_master(jj), evlon_master(jj), lat_tmp, lon_tmp, delta_out = epdelta)
          dist_tmp = sqrt((r_earth - depth_tmp) ** 2 + (r_earth - evdp_master(jj)) ** 2 &
          &    - 2.0_fp * (r_earth - depth_tmp) *      (r_earth - evdp_master(jj)) * cos(epdelta))
          if(nnearest_evmaster .gt. 1) then
            if(evmaster_dist(1) .gt. dist_tmp) then
              nearest_evmaster_index(2, i, j, k) = nearest_evmaster_index(1, i, j, k)
              evmaster_dist(2) = evmaster_dist(1)
              nearest_evmaster_index(1, i, j, k) = jj
              evmaster_dist(1) = dist_tmp
              cycle
            endif
            if(evmaster_dist(2) .gt. dist_tmp) then
              nearest_evmaster_index(2, i, j, k) = jj
              evmaster_dist(2) = dist_tmp
              cycle
            endif
          else
            if(evmaster_dist(1) .gt. dist_tmp) then
              nearest_evmaster_index(1, i, j, k) = jj
              evmaster_dist(1) = dist_tmp
              cycle
            endif
          endif
        enddo
          
        do jj = 1, nsta
          !!do ray tracing with pseudobending scheme
          raypath_lon(1) = lon_tmp
          raypath_lat(1) = lat_tmp
          raypath_dep(1) = depth_tmp
          raypath_lon(nraypath_ini) = stlon(jj)
          raypath_lat(nraypath_ini) = stlat(jj)
          raypath_dep(nraypath_ini) = stdp(jj)
          do ii = 2, nraypath_ini - 1
            raypath_lon(ii) = raypath_lon(ii - 1) &
            &              + (raypath_lon(nraypath_ini) - raypath_lon(1)) / real(nraypath_ini, kind = fp)
            raypath_lat(ii) = raypath_lat(ii - 1) &
            &              + (raypath_lat(nraypath_ini) - raypath_lat(1)) / real(nraypath_ini, kind = fp)
            raypath_dep(ii) = raypath_dep(ii - 1) &
            &              + (raypath_dep(nraypath_ini) - raypath_dep(1)) / real(nraypath_ini, kind = fp)
          enddo
          nraypath = nraypath_ini
          call pseudobending3D(raypath_lon, raypath_lat, raypath_dep, nraypath, ndiv_raypath, &
          &                    velocity(:, :, :, wavetype), lon_str_w, lat_str_s, z_str_min, dlon_str, dlat_str, dz_str, &
          &                    ttime_loc_cor(jj, i, j, k))
          ttime_loc_cor(jj, i, j, k) = ttime_loc_cor(jj, i, j, k) - ttime(jj, nearest_evmaster_index(1, i, j, k))
        enddo  !!station loop
      enddo    !!longitude loop
    enddo      !!latitude loop
  enddo        !!depth loop

!  open(unit = 10, file = "ttime_cor_table.txt")
!  do kk = -nz_cor, nz_cor
!    depth_tmp = evdp_master + dz_cor * real(kk, kind = fp)
!    do jj = -nlat_cor, nlat_cor
!      lat_tmp = evlat_master + dlat_cor * real(jj, kind = fp)
!      do ii = -nlon_cor, nlon_cor
!        lon_tmp = evlon_master + dlon_cor * real(ii, kind = fp)
!        write(10, '(4(f8.4, 1x))') lon_tmp, lat_tmp, depth_tmp, ttime_loc_cor(1, ii, jj, kk)
!      enddo
!      write(10, *) ""
!    enddo
!  enddo
!  close(10)
  !stop


#if defined (WITHOUT_TTIME)
  ttime(1 : nsta, 1 : nev_master) = 0.0_fp
  ttime_loc_cor(1 : nsta, 1 : nlon_cor, 1 : nlat_cor, 1 : nz_cor) = 0.0_fp
#endif

#endif /* (-DWIN || -DSAC) or not */

  !!make amplitude data from win- or sac-formatted waveform data
#if defined (WIN) || defined (SAC)
  !!count nsubevent
  nsubevent = 0
  ot_tmp = ot_begin
  ttime_maxval(i) = 0.0_fp
  do j = 1, nev_master
    do i = 1, nsta
      if(ttime_maxval(i) .lt. ttime(i, j)) ttime_maxval(i) = ttime(i, j)
    enddo
  enddo
  count_loop: do
    if(ot_tmp .gt. ot_end) exit
    do i = 1, nsta
      if(int((ot_tmp - begin(i) + ttime_maxval(i)) + rms_tw) / sampling(i) + 0.5_fp) .gt. npts(i)) then
        exit count_loop
      endif
    enddo
    ot_tmp = ot_tmp + ot_shift
    nsubevent = nsubevent + 1
  enddo count_loop
  write(0, '(a, i0)') "nsubevent = ", nsubevent
  !!calculate rms amplitude of each subevent
  allocate(obsamp_sub(1 : nsta, 1 : nsubevent))
  allocate(data_residual(1 : nsubevent), nsta_use_tmp(1 : nsubevent))

  !!velocity and Qinv at master event location
  do i = 1, nev_master
    lon_index = int((evlon_master(i) - lon_str_w) / dlon_str) + 1
    lat_index = int((evlat_master(i) - lat_str_s) / dlat_str) + 1
    z_index   = int((evdp_master(i)  - z_str_min) / dz_str) + 1
    xgrid(1) = lon_str_w + real(lon_index - 1, kind = fp) * dlon_str; xgrid(2) = xgrid(1) + dlon_str
    ygrid(1) = lat_str_s + real(lat_index - 1, kind = fp) * dlat_str; ygrid(2) = ygrid(1) + dlat_str
    zgrid(1) = z_str_min + real(z_index - 1, kind = fp)   * dz_str;   zgrid(2) = zgrid(1) + dz_str
    val_3d(1 : 2, 1 : 2, 1 : 2) &
    &  = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1, wavetype)
    call linear_interpolation_3d(evlon_master(i), evlat_master(i), evdp_master(i), &
    &                            xgrid, ygrid, zgrid, val_3d, velocity_interpolate(i))
    val_3d(1 : 2, 1 : 2, 1 : 2) &
    &  = qinv(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1, wavetype)
    call block_interpolation_3d(evlon_master(i), evlat_master(i), evdp_master(i), &
    &                           xgrid, ygrid, zgrid, val_3d, qinv_interpolate(i))
  enddo

  !!Set up the observation vector and inversion matrix
#if defined (DAMPED)
  allocate(inversion_matrix(1 : nsta_use * nnearest_evmaster + 4, 1 : 4 + (nnearest_evmaster - 1)), &
  &        obsvector(1 : nsta_use, 1 : nsubevent), &
  &        obsvector_tmp(1 : nsta_use + 4, 1 : nlon_cor, 1 : nlat_cor, 1 : nz_cor))
#else
  allocate(inversion_matrix(1 : nsta_use * nnearest_evmaster, 1 : 4 + (nnearest_evmaster - 1)), &
  &        obsvector(1 : nsta_use, 1 : nsubevent), &
  &        obsvector_tmp(1 : nsta_use, 1 : nlon_cor, 1 : nlat_cor, 1 : nz_cor))
#endif
  allocate(inversion_matrix_copy(1 : ubound(inversion_matrix, 1), 1 : ubound(inversion_matrix, 2)), &
  &        obsvector_copy(1 : ubound(obsvector_tmp, 1)))


  do k = 1, nsubevent
    ot_tmp = ot_begin + ot_shift * real(k - 1, kind = fp)

!!loop for location
    do kk = 1, nz_cor
      do jj = 1, nlat_cor
        do ii = 1, nlon_cor
          obsvector_tmp(1 : ubound(obsvector_tmp, 1), ii, jj, kk) = 0.0_fp
          residual_loc_cor(ii, jj, kk) = huge
 
          do j = 1, nsta
            obsamp_sub(j, k) = 0.0_dp
            icount = 0
            do i = 1, int(rms_tw / sampling(j) + 0.5_fp)
              time_index = int((ot_tmp - begin(j) + ttime(j, nearest_evmaster_index(1, ii, jj, kk)) &
              &                                   + ttime_loc_cor(j, ii, jj, kk)) / sampling(j) + 0.5_fp) + i
              if(time_index .ge. 1 .and. time_index .le. npts(j)) then
                obsamp_sub(j, k) = obsamp_sub(j, k) + waveform_obs(time_index, j) ** 2
                icount = icount + 1
              endif
            enddo
            obsamp_sub(j, k) = sqrt(obsamp_sub(j, k) / real(icount, kind = dp))
          enddo
#endif /* (-DWIN || -DSAC) or not */


          !!check sn ratio
          use_flag_tmp(1 : nsta) = .false.
          nsta_use_tmp(k) = 0

          icount = 0
          do i = 1, nsta
            if(use_flag(i) .eqv. .false.) then
              use_flag_tmp(i) = .false.
              cycle
            else
              if(obsamp_sub(i, k) .gt. snratio_accept * obsamp_noise(i)) then
                use_flag_tmp(i) = .true.
                icount = icount + 1
              else
                use_flag_tmp(i) = .false.
              endif
            endif
          enddo
          if(icount .lt. nsta_use_minimum) cycle

          inversion_matrix(1 : ubound(inversion_matrix, 1), 1 : ubound(inversion_matrix, 2)) = 0.0_fp
          icount = 1
          do i = 1, nsta
            if(use_flag(i) .eqv. .false.) cycle
            if(obsamp_sub(i, k) .eq. 0.0_fp) use_flag_tmp(i) = .false.

            if(use_flag_tmp(i) .eqv. .true.) then
              do k2 = 1, nnearest_evmaster
                nsta_use_tmp(k) = nsta_use_tmp(k) + 1
                obsvector_tmp(icount, ii, jj, kk) &
                &  = log(obsamp_sub(i, k) / obsamp_master(i, nearest_evmaster_index(k2, ii, jj, kk)))
                normal_vector(1 : 3) = [sin(ray_azinc(2, i, nearest_evmaster_index(k2, ii, jj, kk))) &
                &                     * cos(ray_azinc(1, i, nearest_evmaster_index(k2, ii, jj, kk))), &
                &                       sin(ray_azinc(2, i, nearest_evmaster_index(k2, ii, jj, kk))) &
                &                     * sin(ray_azinc(1, i, nearest_evmaster_index(k2, ii, jj, kk))), &
                &                       cos(ray_azinc(2, i, nearest_evmaster_index(k2, ii, jj, kk)))]
                inversion_matrix(icount, 1) = 1.0_fp
                matrix_const = pi * freq * qinv_interpolate(nearest_evmaster_index(k2, ii, jj, kk)) &
                &                        / velocity_interpolate(nearest_evmaster_index(k2, ii, jj, kk)) &
                &                        + 1.0_fp / hypodist(i, nearest_evmaster_index(k2, ii, jj, kk))
                inversion_matrix(icount, 2 : 4) = -1.0_fp * matrix_const * normal_vector(1 : 3)
              enddo
            else
              obsvector_tmp(icount, ii, jj, kk) = 0.0_fp
              inversion_matrix(icount, 1 : 4) = 0.0_fp
            endif
            icount = icount + 1
          enddo

#if defined (DAMPED)
          do i = 1, 4
            inversion_matrix(nsta_use_tmp(k) + i, i) = damp(i)
          enddo
#endif

          !!copy observation vector and inversion matrix
          obsvector_copy(1 : ubound(obsvector_tmp, 1)) = obsvector_tmp(1 : ubound(obsvector_tmp, 1), ii, jj, kk)
          inversion_matrix_copy(1 : ubound(inversion_matrix, 1), 1 : ubound(inversion_matrix, 2)) &
          &  = inversion_matrix(1 : ubound(inversion_matrix, 1), 1 : ubound(inversion_matrix, 2))

          !!calculate least-squares solution
#if defined (MKL)
          call gels(inversion_matrix, obsvector_tmp(:, ii, jj, kk))
#else
          call la_gels(inversion_matrix, obsvector_tmp(:, ii, jj, kk))
#endif

          !!calculate mean data residual
          residual_loc_cor(ii, jj, kk) = 0.0_fp
          icount = 0
          do i = 1, nsta_use
            if(dot_product(inversion_matrix_copy(i, 1 : 4), obsvector_tmp(1 : 4, ii, jj, kk)) .eq. 0.0_fp) cycle
            icount = icount + 1
            residual_loc_cor(ii, jj, kk) = residual_loc_cor(ii, jj, kk) &
            &  + (obsvector_copy(i) - dot_product(inversion_matrix_copy(i, 1 : 4), obsvector_tmp(1 : 4, ii, jj, kk))) ** 2
          enddo
          residual_loc_cor(ii, jj , kk) = residual_loc_cor(ii, jj, kk) / real(icount, kind = dp)
          if(icount .lt. nsta) residual_loc_cor(ii, jj, kk) = real(huge, kind = dp)
          print *, ii, jj, kk, residual_loc_cor(ii, jj, kk), icount
        enddo  !!lon
      enddo    !!lat
    enddo      !!dep
    min_residual = minloc(residual_loc_cor)
    data_residual(k) = minval(residual_loc_cor)
    lon_tmp = evlon_master + dlon_cor * real(min_residual(1) - nlon_cor - 1, kind = fp)
    lat_tmp = evlat_master + dlat_cor * real(min_residual(2) - nlat_cor - 1, kind = fp)
    depth_tmp = evdp_master + dz_cor * real(min_residual(3) - nz_cor - 1, kind = fp)
    print *, min_residual(1) - nlon_cor - 1, min_residual(2) - nlat_cor - 1, min_residual(3) - nz_cor - 1
    print '(4(f8.4, 1x))', lon_tmp, lat_tmp, depth_tmp, data_residual(k)
    obsvector(1 : 4, k) &
    &  = obsvector_tmp(1 : 4, min_residual(1) - 1 - nlon_cor, min_residual(2) - 1 - nlat_cor, min_residual(3) - 1 - nz_cor)
  enddo        !!nsubevent


  !!output result
  open(unit = 10, file = trim(resultfile))

#if defined (WIN) || defined (SAC)
  write(10, '(a)') "# amp_ratio sigma_ampratio longitude sigma_lon latitude sigma_lat depth sigma_depth residual_sum nsta ot_tmp"
  do i = 1, nsubevent
    ot_tmp = ot_begin + ot_shift * real(i - 1, kind = fp)
    delta_lat = (obsvector(2, i) / (r_earth - evdp_master)) * rad2deg
    delta_lon = (obsvector(3, i) / ((r_earth - evdp_master) * sin(pi / 2.0_fp - evlat_master * deg2rad))) * rad2deg
    delta_depth = obsvector(4, i)

    !sigma_amp = exp(obsvector(4 * (i - 1) + 1)) * sqrt(error_matrix(4 * (i - 1) + 1, 4 * (i - 1) + 1)) * 2.0_fp
    !sigma_lat = sqrt(error_matrix(4 * (i - 1) + 2, 4 * (i - 1) + 2)) * rad2deg / (r_earth - evdp_master) * 2.0_fp
    !sigma_lon = sqrt(error_matrix(4 * (i - 1) + 3, 4 * (i - 1) + 3)) &
    !&         * sin(pi / 2.0_fp - evlat_master * deg2rad) * rad2deg / (r_earth - evdp_master) * 2.0_fp
    !sigma_depth = sqrt(error_matrix(4 * (i - 1) + 4, 4 * (i - 1) + 4)) * 2.0_fp
    sigma_amp = 0.0_fp
    sigma_lat = 0.0_fp
    sigma_lon = 0.0_fp
    sigma_depth = 0.0_fp

    write(10, '(9(e15.8, 1x), i0, 1x, f7.1)') &
    &          exp(obsvector(1, i)), sigma_amp, &
    &          evlon_master - delta_lon, sigma_lon, &
    &          evlat_master - delta_lat, sigma_lat, &
    &          evdp_master - delta_depth, sigma_depth, data_residual(i), nsta_use_tmp(i), ot_tmp


    write(0, '(a, i0, a, f8.1)')  "subevent index = ", i, " ot_tmp = ", ot_tmp
    write(0, '(a, 2(e14.7, 1x))') "amp_ratio and sigma_amp = ", exp(obsvector(1, i)), sigma_amp
    write(0, '(a, 2(e14.7, 1x))') "longitude and sigma_lon = ", evlon_master - delta_lon, sigma_lon
    write(0, '(a, 2(e14.7, 1x))') "latitude and sigma_lat = ", evlat_master - delta_lat, sigma_lat
    write(0, '(a, 2(e14.7, 1x))') "depth and sigma_depth = ", evdp_master - delta_depth, sigma_depth
  enddo
  close(10)
#else
  write(10, '(a)') "# amp_ratio sigma_ampratio longitude sigma_lon latitude sigma_lat depth sigma_depth residual_sum nsta"
  do i = 1, nsubevent
    delta_lat = (obsvector(2, i) / (r_earth - evdp_master)) * rad2deg
    delta_lon = (obsvector(3, i) / ((r_earth - evdp_master) * sin(pi / 2.0_fp - evlat_master * deg2rad))) * rad2deg
    delta_depth = obsvector(4, i)

    !sigma_amp = exp(obsvector(4 * (i - 1) + 1)) * sqrt(error_matrix(4 * (i - 1) + 1, 4 * (i - 1) + 1)) * 2.0_fp
    !sigma_lat = sqrt(error_matrix(4 * (i - 1) + 2, 4 * (i - 1) + 2)) * rad2deg / (r_earth - evdp_master) * 2.0_fp
    !sigma_lon = sqrt(error_matrix(4 * (i - 1) + 3, 4 * (i - 1) + 3)) &
    !&         * sin(pi / 2.0_fp - evlat_master * deg2rad) * rad2deg / (r_earth - evdp_master) * 2.0_fp
    !sigma_depth = sqrt(error_matrix(4 * (i - 1) + 4, 4 * (i - 1) + 4)) * 2.0_fp
    sigma_amp = 0.0_fp
    sigma_lat = 0.0_fp
    sigma_lon = 0.0_fp
    sigma_depth = 0.0_fp

    write(10, '(9(e15.8, 1x), i0)') &
    &          exp(obsvector(1, i)), sigma_amp, &
    &          evlon_master - delta_lon, sigma_lon, &
    &          evlat_master - delta_lat, sigma_lat, &
    &          evdp_master - delta_depth, sigma_depth, data_residual(i), nsta_use_tmp(i)

    write(0, '(a, i0)')           "subevent index = ", i
    write(0, '(a, 2(e14.7, 1x))') "amp_ratio and sigma_amp = ", exp(obsvector(1, i)), sigma_amp
    write(0, '(a, 2(e14.7, 1x))') "longitude and sigma_lon = ", evlon_master - delta_lon, sigma_lon
    write(0, '(a, 2(e14.7, 1x))') "latitude and sigma_lat = ", evlat_master - delta_lat, sigma_lat
    write(0, '(a, 2(e14.7, 1x))') "depth and sigma_depth = ", evdp_master - delta_depth, sigma_depth
  enddo
  close(10)
#endif
    

  stop
end program AmplitudeSourceLocation_masterevent

