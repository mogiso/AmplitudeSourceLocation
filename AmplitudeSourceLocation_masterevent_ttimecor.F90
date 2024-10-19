program AmplitudeSourceLocation_masterevent_ttimecor
  !!Relative amplitude source location (master event) method using multiple reference events
  !!using depth-dependent 1D velocity structure, 3D heterogeneous attenuation structure
  !!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2024
  !!License  : MIT License (https://opensource.org/licenses/MIT)

  !!-DWIN: use win-formatted waveform file for input
  !!-DSAC: use sac-formatted(binary) waveform file  for input
  !!default: txt format

  !!calculate traveltime correction table for determining events which
  !interevent distance are long

  use nrtype,               only : fp, dp
  use constants,            only : rad2deg, deg2rad, pi, r_earth
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

  type location
    real(kind = fp) :: lon, lat, xeast, ynorth, depth, vel, qinv, sigma_x, sigma_y, sigma_depth
  end type location
  type evsourceamp
    real(kind = fp) :: sourceamp, sigma_sourceamp
    integer         :: evindex
  end type evsourceamp

  integer,            parameter :: wavetype = 2          !!1 for P-wave, 2 for S-wave
  integer,            parameter :: nsta_use_minimum = 5
  real(kind = fp),    parameter :: snratio_accept = 0.0_fp
  !!Range for velocity and attenuation structure
  !!whole Japan
  real(kind = fp),    parameter :: lon_str_w = 122.0_fp, lon_str_e = 150.0_fp
  real(kind = fp),    parameter :: lat_str_s =  23.0_fp, lat_str_n = 46.0_fp
  real(kind = fp),    parameter :: z_str_min =  -3.0_fp, z_str_max = 100.0_fp
  real(kind = fp),    parameter :: dlon_str  =   0.1_fp, dlat_str  = 0.1_fp, dz_str = 1.0_fp
  integer,            parameter :: nlon_str  = int((lon_str_e - lon_str_w) / dlon_str) + 1
  integer,            parameter :: nlat_str  = int((lat_str_n - lat_str_s) / dlat_str) + 1
  integer,            parameter :: nz_str    = int((z_str_max - z_str_min) / dz_str)   + 1

  !!traveltime correction table
  real(kind = fp),    parameter :: lon_cor_w = 130.4_fp, lon_cor_e = 131.7_fp
  real(kind = fp),    parameter :: lat_cor_s = 32.4_fp,  lat_cor_n = 33.5_fp
  real(kind = fp),    parameter :: z_cor_min = 2.0_fp,   z_cor_max = 20.0_fp
  real(kind = fp),    parameter :: dlon_cor = 0.04_fp, dlat_cor = 0.04_fp, dz_cor = 2.0_fp
  integer,            parameter :: nlon_cor = int((lon_cor_e - lon_cor_w) / dlon_cor) + 1
  integer,            parameter :: nlat_cor = int((lat_cor_n - lat_cor_s) / dlat_cor) + 1
  integer,            parameter :: nz_cor   = int((z_cor_max - z_cor_min) / dz_cor)   + 1

#if defined (WIN) /* use win-format waveform file for input waveforms */
  real(kind = dp),    parameter   :: order = 1.0e+6_dp
  real(kind = fp),    allocatable :: sampling(:), begin(:)
  real(kind = dp),    allocatable :: waveform_obs(:, :)
  integer,            allocatable :: ikey(:), sampling_int(:), waveform_obs_int(:, :), npts_win(:, :), npts(:)
  character(len = 4), allocatable :: st_winch(:)
  type(winch__hdr),   allocatable :: chtbl(:)
  integer                         :: nsec, tim
  character(len = 129)            :: win_filename, win_chfilename
  character(len = 10)             :: cmpnm
#endif

#if defined (SAC) /* use sac-format wavefome files (NVHDR=6) */
#if defined (TESTDATA)
  real(kind = dp),    parameter   :: order = 1.0_dp
#else
  real(kind = dp),    parameter   :: order = 1.0_dp
#endif
  character(len = 3), parameter   :: sacfile_extension = "sac"
  real(kind = dp),    allocatable :: waveform_obs(:, :)
  real(kind = fp),    allocatable :: begin(:), sampling(:), stime(:)
  integer,            allocatable :: npts(:)
  character(len = 129)            :: sacfile, sacfile_index
  character(len = 10)             :: cmpnm
#endif

  real(kind = fp),    allocatable :: stlon(:), stlat(:), stdp(:), ttime(:, :), ttime_cor(:, :), ttime_loc_cor(:, :, :, :)
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

  real(kind = fp),    parameter :: alt_to_depth = -1.0e-3_fp
  real(kind = dp),    parameter :: huge = 1.0e+6_dp

  type(location),     allocatable :: evloc_master(:), evloc_sub(:)
  type(evsourceamp),  allocatable :: evamp_sub(:, :)

  real(kind = fp)               :: velocity(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                qinv(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                val_2d(1 : 2, 1 : 2), val_3d(1 : 2, 1 : 2, 1 : 2), &
  &                                xgrid(1 : 2), ygrid(1 : 2), zgrid(1 : 2), normal_vector(1 : 3) 
  real(kind = dp),  allocatable :: topography(:, :), lon_topo(:), lat_topo(:)
  real(kind = fp),  allocatable :: obsamp_master(:, :), obsamp_sub(:, :), obsamp_noise(:), &
  &                                hypodist(:, :), ray_azinc(:, :, :), evdist_master(:), &
  &                                obsvector(:), obsvector_org(:), obsvector_min(:), &
  &                                inversion_matrix(:, :), inversion_matrix_org(:, :), inversion_matrix_min(:, :), &
  &                                data_residual(:), error_matrix(:, :)
  integer,          allocatable :: ipiv(:), nsta_use_min(:), evindex_master(:, :, :, :)
  real(kind = fp)               :: epdist, epdelta, mean_lon, mean_lat, mean_depth, lon_tmp, lat_tmp, depth_tmp, dist_tmp, &
  &                                ttime_tmp, ttime_maxval, matrix_const, &
  &                                siteamp_tmp, residual_sum, mean_residual, sigma_residual, sigma_lat, sigma_lon
  real(kind = dp)               :: topography_interpolate, dlon_topo, dlat_topo, freq
  integer                       :: nlon_topo, nlat_topo, nsta, nev_master_max, nsubevent, lon_index, lat_index, z_index, &
  &                                i, j, k, ii, jj, kk, k2, icount, nsta_use, ios, time_index, nev_master, ndata, nsta_use_tmp
  character(len = 129)          :: topo_grd, station_param, masterevent_param, subevent_param, resultfile, resultfile_tmp
  character(len = 20)           :: cfmt, nsta_c
  character(len = 3)            :: nev_master_t

  !!Psuedobending parameters
  integer,                parameter :: ndiv_raypath = 10
  integer,                parameter :: nraypath_ini = 4
  real(kind = fp)                   :: raypath_lon((nraypath_ini - 1) * 2 ** ndiv_raypath + 1), &
  &                                    raypath_lat((nraypath_ini - 1) * 2 ** ndiv_raypath + 1), &
  &                                    raypath_dep((nraypath_ini - 1) * 2 ** ndiv_raypath + 1)
  integer                           :: nraypath

  icount = iargc()

#if defined (WIN)
  if(icount .ne. 15) then
    write(0, '(a)', advance="no") "usage: ./asl_masterevent "
    write(0, '(a)', advance="no") "(topography_grd) (station_param_file) (masterevent_param_file) (win_waveform) (win_chfile) "
    write(0, '(a)', advance="no") "(component name) (fl) (fh) (fs) (ot_begin) (ot_end) (ot_shift) (rms_time_window_length) "
    write(0, '(a)')               "(nev_master) (result_file)"
    error stop
  endif
  call getarg(1, topo_grd)
  call getarg(2, station_param)
  call getarg(3, masterevent_param)
  call getarg(4, win_filename)
  call getarg(5, win_chfilename)
  call getarg(6, cmpnm)
  call getarg(7, fl_t)         ; read(fl_t, *) fl
  call getarg(8, fh_t)         ; read(fh_t, *) fh
  call getarg(9, fs_t)         ; read(fs_t, *) fs
  call getarg(10, ot_begin_t)  ; read(ot_begin_t, *) ot_begin
  call getarg(11, ot_end_t)    ; read(ot_end_t, *) ot_end
  call getarg(12, ot_shift_t)  ; read(ot_shift_t, *) ot_shift
  call getarg(13, rms_tw_t)    ; read(rms_tw_t, *) rms_tw
  call getarg(14, nev_master_t); read(nev_master_t, *) nev_master
  call getarg(15, resultfile)
#elif defined (SAC)
  if(icount .ne. 14) then
    write(0, '(a)', advance="no") "usage: ./asl_masterevent "
    write(0, '(a)', advance="no") "(topography_grd) (station_param_file) (masterevent_param_file) (sacfile_index) "
    write(0, '(a)', advance="no") "(component_name) (fl) (fh) (fs) (ot_begin) (ot_end) (ot_shift) (rms_time_window_length) "
    write(0, '(a)')               "(nev_master) (result_file)"
    error stop
  endif
  call getarg(1, topo_grd)
  call getarg(2, station_param)
  call getarg(3, masterevent_param)
  call getarg(4, sacfile_index)
  call getarg(5, cmpnm)
  call getarg(6, fl_t)         ; read(fl_t, *) fl
  call getarg(7, fh_t)         ; read(fh_t, *) fh
  call getarg(8, fs_t)         ; read(fs_t, *) fs
  call getarg(9, ot_begin_t)   ; read(ot_begin_t, *) ot_begin
  call getarg(10, ot_end_t)    ; read(ot_end_t, *) ot_end
  call getarg(11, ot_shift_t)  ; read(ot_shift_t, *) ot_shift
  call getarg(12, rms_tw_t)    ; read(rms_tw_t, *) rms_tw
  call getarg(13, nev_master_t); read(nev_master_t, *) nev_master
  call getarg(14, resultfile)
#else
  if(icount .ne. 7) then
    write(0, '(a)', advance="no") "usage: ./asl_masterevent "
    write(0, '(a)', advance="no") "(topography_grd) (station_param_file) (masterevent_param_file) (subevent_param_file) "
    write(0, '(a)')               "(frequency) (nev_master) (result_file)"
    error stop
  endif
  call getarg(1, topo_grd)
  call getarg(2, station_param)
  call getarg(3, masterevent_param)
  call getarg(4, subevent_param)
  call getarg(5, freq_t)      ; read(freq_t, *) freq
  call getarg(6, nev_master_t); read(nev_master_t, *) nev_master
  call getarg(7, resultfile)
#endif

#if defined (DAMPED)
  write(0, '(a)') "-DDAMPED set"
#endif

#if defined (WIN) || defined (SAC)
  write(0, '(a, 3(f5.2, 1x))') "Bandpass filter parameter fl, fh, fs (Hz) = ", fl, fh, fs
  freq = (fl + fh) * 0.5_dp
#endif

  !!set velocity/attenuation structure
  call set_velocity(z_str_min, dz_str, velocity, qinv)

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
    read(10, *) stlon(i), stlat(i), stdp(i), stname(i), use_flag(i), &
    &           ttime_cor(i, 1), ttime_cor(i, 2), siteamp_tmp, obsamp_noise(i)
    write(0, '(a, i0, a, f9.4, a, f8.4, a, f6.3, 1x, a7, l2)') &
    &     "station(", i, ") lon(deg) = ", stlon(i), &
    &                     " lat(deg) = ", stlat(i), " depth(km) = ", stdp(i), &
    &                       trim(stname(i)), use_flag(i)
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
  nev_master_max = 0
  read(10, *)
  do
    read(10, *, iostat = ios)
    if(ios .ne. 0) exit
    nev_master_max = nev_master_max + 1
  enddo
  nev_master_max = nev_master_max / 2
  allocate(evloc_master(1 : nev_master_max))
  allocate(obsamp_master(1 : nsta, 1 : nev_master_max))
  rewind(10)
  read(10, *)
  mean_lon = 0.0_fp
  mean_lat = 0.0_fp
  mean_depth = 0.0_fp
  do j = 1, nev_master_max
    read(10, *) evloc_master(j)%lon, evloc_master(j)%lat, evloc_master(j)%depth
    read(10, *) (obsamp_master(i, j), i = 1, nsta)
    mean_lon = mean_lon + evloc_master(j)%lon
    mean_lat = mean_lat + evloc_master(j)%lat
    mean_depth = mean_depth + evloc_master(j)%depth
  enddo
  close(10)
  mean_lon = mean_lon / real(nev_master_max, kind = fp)
  mean_lat = mean_lat / real(nev_master_max, kind = fp)
  mean_depth = mean_depth / real(nev_master_max, kind = fp)
  write(0, '(a, 3(f9.4, 1x))') "mean lon, lat, depth = ", mean_lon, mean_lat, mean_depth
  do i = 1, nev_master_max
    evloc_master(i)%xeast = (evloc_master(i)%lon - mean_lon) * deg2rad &
    &                     * (r_earth - mean_depth) * sin(pi * 0.5_fp - mean_lat * deg2rad)
    evloc_master(i)%ynorth = (evloc_master(i)%lat - mean_lat) * deg2rad * (r_earth - mean_depth)

    !!velocity and Qinv at master event location
    lon_index = int((evloc_master(i)%lon   - lon_str_w) / dlon_str) + 1
    lat_index = int((evloc_master(i)%lat   - lat_str_s) / dlat_str) + 1
    z_index   = int((evloc_master(i)%depth - z_str_min) / dz_str)   + 1
    xgrid(1) = lon_str_w + real(lon_index - 1, kind = fp) * dlon_str; xgrid(2) = xgrid(1) + dlon_str
    ygrid(1) = lat_str_s + real(lat_index - 1, kind = fp) * dlat_str; ygrid(2) = ygrid(1) + dlat_str
    zgrid(1) = z_str_min + real(z_index   - 1, kind = fp) * dz_str  ; zgrid(2) = zgrid(1) + dz_str
    val_3d(1 : 2, 1 : 2, 1 : 2) &
    &  = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1, wavetype)
    call linear_interpolation_3d(evloc_master(i)%lon, evloc_master(i)%lat, evloc_master(i)%depth, &
    &                            xgrid, ygrid, zgrid, val_3d, evloc_master(i)%vel)
    val_3d(1 : 2, 1 : 2, 1 : 2) &
    &  = qinv(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1, wavetype)
    call block_interpolation_3d(evloc_master(i)%lon, evloc_master(i)%lat, evloc_master(i)%depth, &
    &                           xgrid, ygrid, zgrid, val_3d, evloc_master(i)%qinv)
  enddo 

  if(nev_master .gt. nev_master_max) nev_master = nev_master_max
  allocate(evindex_master(1 : nev_master, 1 : nlon_cor, 1 : nlat_cor, 1 : nz_cor))
  allocate(evdist_master(1 : nev_master))
  allocate(ttime(1 : nsta, 1 : nev_master_max))

#if defined (WIN) || defined (SAC)
  allocate(sampling(1 : nsta), npts(1 : nsta), begin(1 : nsta))

#if defined (WIN)
  allocate(st_winch(nsta), sampling_int(1 : nsta))
  !!read channel table
  call winch__read_tbl(trim(win_chfilename), chtbl)
  !!find chid from given stname(nsta) and cmpnm
  allocate(ikey(nsta))
  do i = 1, nsta
    call winch__st2chid(chtbl, stname(i), trim(cmpnm), st_winch(i), ikey(i))
    if(ikey(i) .eq. 0) then
      write(0, '(5(a, 1x))') "station/comp =", trim(stname(i)), trim(cmpnm), &
      &                      "does not exist in", trim(win_chfilename)
      error stop
    endif
    write(0, '(6a)') "station name = ", trim(stname(i)), " comp = ", trim(cmpnm), " chid = ", st_winch(i)
  enddo
  !!read all waveform data from win-formatted file
  call win__read_file(trim(win_filename), st_winch, sampling_int, nsec, tim, waveform_obs_int, npts_win)
  sampling(1 : nsta) = 1.0_fp / real(sampling_int(1 : nsta), kind = fp)
  npts(1 : nsta)     = sampling_int(1 : nsta) * nsec
  begin(1 : nsta)    = 0.0_fp
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
    sacfile = trim(sacfile_index) // "." // trim(stname(i)) // "." // trim(cmpnm) // "." // sacfile_extension
    call read_sachdr(sacfile, delta=sampling(i), npts=npts(i), begin=begin(i), t0=stime(i))
  enddo
  allocate(waveform_obs(maxval(npts), nsta))
  do i = 1, nsta
    sacfile = trim(sacfile_index) // "." // trim(stname(i)) // "." // trim(cmpnm) // "." // sacfile_extension
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
#endif
#endif /* -DWIN || -DSAC */



  !!calculate ray length, pulse width, unit vector of ray incident
  write(0, '(a)') "calculate ray length, pulse width, and ray incident vector for master events"

  !!check whether the depth of master event location is lower than the topo
  do i = 1, nev_master_max
    lon_index = int((evloc_master(i)%lon - lon_topo(1)) / dlon_topo) + 1
    lat_index = int((evloc_master(i)%lat - lat_topo(1)) / dlat_topo) + 1
    xgrid(1 : 2) = [lon_topo(lon_index), lon_topo(lon_index + 1)]
    ygrid(1 : 2) = [lat_topo(lat_index), lat_topo(lat_index + 1)]
    val_2d(1 : 2, 1 : 2) = topography(lon_index : lon_index + 1, lat_index : lat_index + 1)
    call linear_interpolation_2d(evloc_master(i)%lon, evloc_master(i)%lat, xgrid, ygrid, val_2d, topography_interpolate)
    if(evloc_master(i)%depth .lt. topography_interpolate) then
      write(0, '(a, f5.2, a)') "Depth of master event ", evloc_master(i)%depth, " is higher than the altitude there."
      error stop
    endif
  enddo

  allocate(hypodist(1 : nsta, 1 : nev_master_max), ray_azinc(1 : 2, 1 : nsta, 1 : nev_master_max))
  
  !!calculate traveltime from each master event to each station        
  masterevent_loop: do kk = 1, nev_master_max
    station_loop: do jj = 1, nsta
      !!calculate azimuth and hypocentral distance
      call greatcircle_dist(evloc_master(kk)%lat, evloc_master(kk)%lon, stlat(jj), stlon(jj), &
      &                     distance  = epdist, &
      &                     azimuth   = ray_azinc(1, jj, kk), &
      &                     delta_out = epdelta)
      lon_index = int((stlon(jj) - lon_str_w) / dlon_str) + 1
      lat_index = int((stlat(jj) - lat_str_s) / dlat_str) + 1
      z_index   = int((stdp(jj)  - z_str_min) / dz_str)   + 1
      !print *, lon_sta(jj), lon_w + real(lon_index - 1) * dlon
      !print *, lat_sta(jj), lat_s + real(lat_index - 1) * dlat
      hypodist(jj, kk) = sqrt((r_earth - evloc_master(kk)%depth) ** 2 + (r_earth - stdp(jj)) ** 2 &
      &            - 2.0_fp * (r_earth - evloc_master(kk)%depth) *      (r_earth - stdp(jj)) * cos(epdelta))

#if defined (V_CONST)
      !!homogeneous structure: ray incident angle is calculated using cosine function (assuming cartesian coordinate)
      ray_azinc(2, jj, kk) = acos(epdist / hypodist(jj, kk)) + pi * 0.5_fp
#else
      !!do ray tracing with pseudobending scheme
      raypath_lon(1) = evloc_master(kk)%lon
      raypath_lat(1) = evloc_master(kk)%lat
      raypath_dep(1) = evloc_master(kk)%depth
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
      ttime(jj, kk) = ttime_tmp + ttime_cor(jj, wavetype)

      !write(0, '(a, 4(f8.4, 1x))') "masterevent lon, lat, depth = ", &
      !&                             evloc_master(kk)%lon, evloc_master(kk)%lat, evloc_master(kk)%depth
      !write(0, '(a, 3(f8.4, 1x))') "station lon, lat, depth = ", stlon(jj), stlat(jj), stdp(jj)
      !write(0, '(a, 2(f8.4, 1x))') "ray azimuth and inc_angle (deg) = ", &
      !&                             ray_azinc(1, jj, kk) * rad2deg, ray_azinc(2, jj, kk) * rad2deg
      !write(0, '(a, f5.2)') "traveltime (s) = ", ttime(jj, kk)

#endif /* -DV_CONST or not */
    enddo station_loop
  enddo masterevent_loop

  !!make traveltime correction table
  write(0, '(a)') "calculate traveltime correction table"
  allocate(ttime_loc_cor(1 : nsta, 1 : nlon_cor, 1 : nlat_cor, 1 : nz_cor))
  do k = 1, nz_cor
    depth_tmp = z_cor_min + dz_cor * real(k, kind = fp)
    do j = 1, nlat_cor
      lat_tmp = lat_cor_s + dlat_cor * real(j, kind = fp)
      do i = 1, nlon_cor
        lon_tmp = lon_cor_w + dlon_cor * real(i, kind = fp)
        evdist_master(1 : nev_master) = huge
        evindex_master(1 : nev_master, i, j, k) = 1
        do kk = 1, nev_master_max
          call greatcircle_dist(evloc_master(kk)%lat, evloc_master(kk)%lon, lat_tmp, lon_tmp, delta_out = epdelta)
          dist_tmp = sqrt((r_earth - depth_tmp) ** 2      + (r_earth - evloc_master(kk)%depth) ** 2 &
          &             - (r_earth - depth_tmp) *  2.0_fp * (r_earth - evloc_master(kk)%depth) * cos(epdelta))
          do jj = 1, nev_master
            if(dist_tmp .le. evdist_master(jj)) then
              do ii = nev_master, jj + 1, -1
                evdist_master(ii) = evdist_master(ii - 1)
                evindex_master(ii, i, j, k) = evindex_master(ii - 1, i, j, k)
              enddo
              evdist_master(jj) = dist_tmp
              evindex_master(jj, i, j, k) = kk
              exit
            endif
          enddo
        enddo
          
        !!calculate travel time from each grid to each station
        do jj = 1, nsta
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
          !!use travel time of closest master event
          !ttime_loc_cor(jj, i, j, k) = ttime_loc_cor(jj, i, j, k) - ttime(jj, evindex_master(1, i, j, k))
        enddo  !!station loop
      enddo    !!longitude loop
    enddo      !!latitude loop
  enddo        !!depth loop

  open(unit = 10, file = "ttime_cor_table.txt")
  do kk = 1, nz_cor
    depth_tmp = z_cor_min + dz_cor * real(kk - 1, kind = fp)
    do jj = 1, nlat_cor
      lat_tmp = lat_cor_s + dlat_cor * real(jj - 1, kind = fp)
      do ii = 1, nlon_cor
        lon_tmp = lon_cor_w + dlon_cor * real(ii - 1, kind = fp)
        write(10, '(4(f8.4, 1x))') lon_tmp, lat_tmp, depth_tmp, ttime_loc_cor(1, ii, jj, kk)
      enddo
      write(10, *) ""
    enddo
  enddo
  close(10)
  !stop

#if defined (WITHOUT_TTIME)
  ttime(1 : nsta, 1 : nev_master_max) = 0.0_fp
  ttime_loc_cor(1 : nsta, 1 : nlon_cor, 1 : nlat_cor, 1 : nz_cor) = 0.0_fp
#endif

#if defined (WIN) || defined (SAC)
  !!make amplitude data from win- or sac-formatted waveform data
  !!count nsubevent
  nsubevent = 0
  ot_tmp = ot_begin
  count_loop: do
    if(ot_tmp .gt. ot_end) exit
    do i = 1, nsta
      ttime_maxval = maxval(ttime_loc_cor(i, :, :, :))
      if(int((ot_tmp - begin(i) + ttime_maxval + rms_tw) / sampling(i) + 0.5_fp) .gt. npts(i)) then
        exit count_loop
      endif
    enddo
    ot_tmp = ot_tmp + ot_shift
    nsubevent = nsubevent + 1
  enddo count_loop
  write(0, '(a, i0)') "nsubevent = ", nsubevent
  allocate(obsamp_sub(1 : nsta, 1 : nsubevent))
  allocate(data_residual(1 : nsubevent), nsta_use_min(1 : nsubevent))
  do j = 1, nsta
    obsamp_noise(j) = 0.0_fp
    icount = 0
    do i = int(rms_tw / sampling(j)) + 1, int(rms_tw / sampling(j)) * 2 
      obsamp_noise(j) = obsamp_noise(j) + waveform_obs(i, j) ** 2
      icount = icount + 1
    enddo
    obsamp_noise(j) = sqrt(obsamp_noise(j) / real(icount, kind = fp))
  enddo
#else
  !!read amplitudes of subevents from text file
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
  allocate(obsamp_sub(1 : nsta, 1 : nsubevent))
  allocate(data_residual(1 : nsubevent), nsta_use_min(1 : nsubevent))
  do j = 1, nsubevent
    read(10, *) (obsamp_sub(i, j), i = 1, nsta)
  enddo
  close(10)
  obsamp_noise(1 : nsta) = 0.0_fp
#endif /* -DWIN || -DSAC */

  allocate(evamp_sub(1 : nev_master, 1 : nsubevent), evloc_sub(1 : nsubevent))

  !!Set up the observation vector and inversion matrix
#if defined (DAMPED)
  allocate(inversion_matrix(1 : nsta_use * nev_master + 3 * nev_master, 1 : 3 + nev_master), &
  &               obsvector(1 : nsta_use * nev_master + 3 * nev_master))
#else
  allocate(inversion_matrix(1 : nsta_use * nev_master, 1 : 3 + nev_master), &
  &               obsvector(1 : nsta_use * nev_master))
#endif
  allocate(inversion_matrix_min(1 : ubound(inversion_matrix, 1), 1 : ubound(inversion_matrix, 2)), &
  &        inversion_matrix_org(1 : ubound(inversion_matrix, 1), 1 : ubound(inversion_matrix, 2)), &
  &        obsvector_min(1 : ubound(obsvector, 1)), obsvector_org(1 : ubound(obsvector, 1)))


  do k = 1, nsubevent
    data_residual(k) = huge
    nsta_use_min(k) = 0

    write(0, '(a, i0)') "k = ", k
#if defined (WIN) || defined (SAC)
    ot_tmp = ot_begin + ot_shift * real(k - 1, kind = fp)
    write(0, '(a, f8.1)') "origin time = ", ot_tmp
#endif

!!loop for location
    do kk = 1, nz_cor
      do jj = 1, nlat_cor
        do ii = 1, nlon_cor
 
#if defined (WIN) || defined (SAC)
          do j = 1, nsta
            obsamp_sub(j, k) = 0.0_dp
            icount = 0
            do i = 1, int(rms_tw / sampling(j) + 0.5_fp)
              time_index = int((ot_tmp - begin(j) &
              !&               + ttime(j, evindex_master(1, ii, jj, kk)) &
              &               + ttime_loc_cor(j, ii, jj, kk)) / sampling(j) + 0.5_fp) + i
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
          obsvector(1 : ubound(obsvector, 1)) = 0.0_fp
          icount = 1
          nsta_use_tmp = 0
          do j = 1, nsta
            if(use_flag(j) .eqv. .false.) cycle
            nsta_use_tmp = nsta_use_tmp + 1

            if(use_flag_tmp(j) .eqv. .true.) then
              do i = 1, nev_master
                normal_vector(1 : 3) = [sin(ray_azinc(2, j, evindex_master(i, ii, jj, kk)))  &
                &                     * cos(ray_azinc(1, j, evindex_master(i, ii, jj, kk))), &  !!latitude (ynorth)
                &                       sin(ray_azinc(2, j, evindex_master(i, ii, jj, kk)))  &
                &                     * sin(ray_azinc(1, j, evindex_master(i, ii, jj, kk))), &  !!longitude (xeast)
                &                       cos(ray_azinc(2, j, evindex_master(i, ii, jj, kk)))]    !!depth (down+)
                matrix_const = pi * freq * evloc_master(evindex_master(i, ii, jj, kk))%qinv &
                &                        / evloc_master(evindex_master(i, ii, jj, kk))%vel  &
                &                        + 1.0_fp / hypodist(j, evindex_master(i, ii, jj, kk))
                inversion_matrix(icount, 1) = 1.0_fp * matrix_const * normal_vector(1)
                inversion_matrix(icount, 2) = 1.0_fp * matrix_const * normal_vector(2)
                inversion_matrix(icount, 3) = 1.0_fp * matrix_const * normal_vector(3)
                inversion_matrix(icount, 3 + i) = 1.0_fp

                obsvector(icount) &
                &  = log(obsamp_sub(j, k) / obsamp_master(j, evindex_master(i, ii, jj, kk))) &
                &  +     matrix_const * (evloc_master(evindex_master(i, ii, jj, kk))%ynorth * normal_vector(1) &
                &                     +  evloc_master(evindex_master(i, ii, jj, kk))%xeast  * normal_vector(2) &
                &                     +  evloc_master(evindex_master(i, ii, jj, kk))%depth  * normal_vector(3))

                icount = icount + 1
              enddo
            endif
          enddo
          ndata = icount - 1

#if defined (DAMPED)
          do j = 1, nev_master
            do i = 1, 3
              if(i .eq. 3) then
                inversion_matrix(icount, i) = 1.0_fp
                obsvector(icount) = evloc_master(evindex_master(j, ii, jj, kk))%depth
              endif
              icount = icount + 1
            enddo
          enddo
#endif

          !!copy observation vector and inversion matrix
          obsvector_org(1 : ubound(obsvector, 1)) = obsvector(1 : ubound(obsvector, 1))
          inversion_matrix_org(1 : ubound(inversion_matrix, 1), 1 : ubound(inversion_matrix, 2)) &
          &  = inversion_matrix(1 : ubound(inversion_matrix, 1), 1 : ubound(inversion_matrix, 2))

          !!calculate least-squares solution
#if defined (MKL)
          call gels(inversion_matrix(1 : icount - 1, :) , obsvector(1 : icount - 1))
#else
          call la_gels(inversion_matrix(1 : icount - 1, :), obsvector(1 : icount - 1))
#endif

          !!calculate mean data residual
          residual_sum = 0.0_fp
          icount = 0
          do i = 1, ndata
            if(dot_product(inversion_matrix_org(i, 1 : 3 + nev_master), obsvector(1 : 3 + nev_master)) .eq. 0.0_fp) &
            &  cycle
            icount = icount + 1
            residual_sum = residual_sum &
            &  + (obsvector_org(i) &
            &  -  dot_product(inversion_matrix_org(i, 1 : 3 + nev_master), obsvector(1 : 3 + nev_master))) ** 2
          enddo
          residual_sum = residual_sum / real(icount, kind = dp)
          if(icount .lt. nsta) residual_sum = real(huge, kind = dp)
          !print *, ii, jj, kk, residual_sum, icount
          if(residual_sum .lt. data_residual(k)) then
            data_residual(k) = residual_sum
            nsta_use_min(k)  = nsta_use_tmp
            inversion_matrix_min(1 : ubound(inversion_matrix, 1), 1 : ubound(inversion_matrix, 2)) &
            &  = inversion_matrix_org(1 : ubound(inversion_matrix, 1), 1 : ubound(inversion_matrix, 2))
            obsvector_min(1 : ubound(obsvector, 1)) = obsvector(1 : ubound(obsvector, 1))
            evloc_sub(k)%ynorth = obsvector(1)
            evloc_sub(k)%xeast  = obsvector(2)
            evloc_sub(k)%depth  = obsvector(3)
            do i = 1, nev_master
              evamp_sub(i, k)%sourceamp = exp(obsvector(3 + i))
              evamp_sub(i, k)%evindex   = evindex_master(i, ii, jj, kk)
            enddo
          endif
        enddo  !!lon
      enddo    !!lat
    enddo      !!dep
    if(nsta_use_min(k) .lt. nsta_use_minimum) cycle

    evloc_sub(k)%lat = mean_lat + (evloc_sub(k)%ynorth / (r_earth - mean_depth) * rad2deg)
    evloc_sub(k)%lon = mean_lon &
    &                + ((evloc_sub(k)%xeast / ((r_earth - mean_depth) * sin(pi * 0.5_fp - mean_lat * deg2rad))) * rad2deg)

    !print '(4(f8.4, 1x))', lon_tmp, lat_tmp, depth_tmp, data_residual(k)

    !!calculate residual and its variance
    mean_residual = 0.0_fp
    do i = 1, nsta_use_min(k) * nev_master
      mean_residual = mean_residual + obsvector_org(i) &
      &             - dot_product(inversion_matrix_min(i, 1 : 3 + nev_master), obsvector_min(1 : 3 + nev_master))
    enddo
    mean_residual = mean_residual / real(nsta_use_min(k) * nev_master, kind = fp)
    sigma_residual = 0.0_fp
    do i = 1, nsta_use_min(k) * nev_master
      sigma_residual = sigma_residual &
      &  + (obsvector_org(i) - dot_product(inversion_matrix_min(i, 1 : 3 + nev_master), obsvector_min(1 : 3 + nev_master)) &
      &  - mean_residual) ** 2
    enddo
    sigma_residual = sigma_residual / real(nsta_use_min(k) * nev_master - 1, kind = fp)
   
    allocate(error_matrix(1 : 3 + nev_master, 1 : 3 + nev_master))
    allocate(ipiv(1 : ubound(error_matrix, 1)))
    error_matrix(1 : 3 + nev_master, 1 : 3 + nev_master) = matmul(transpose(inversion_matrix_min), inversion_matrix_min)
#if defined (MKL)
    call getrf(error_matrix, ipiv)
    call getri(error_matrix, ipiv)
#else
    call la_getrf(error_matrix, ipiv)
    call la_getri(error_matrix, ipiv)
#endif
    evloc_sub(k)%sigma_x     = sqrt(sigma_residual * error_matrix(1, 1))
    evloc_sub(k)%sigma_y     = sqrt(sigma_residual * error_matrix(2, 2))
    evloc_sub(k)%sigma_depth = sqrt(sigma_residual * error_matrix(3, 3))
    do i = 1, nev_master
      evamp_sub(i, k)%sigma_sourceamp = sqrt(sigma_residual * error_matrix(3 + i, 3 + i))
    enddo

    deallocate(error_matrix, ipiv)
  enddo        !!nsubevent

  !!output result
  resultfile_tmp = trim(resultfile) // "_evloc_sub.txt"
  open(unit = 10, file = trim(resultfile_tmp))

#if defined (WIN) || defined (SAC)
  write(10, '(a)') "# longitude sigma_lon latitude sigma_lat depth sigma_depth residual_sum nsta ot_tmp"
  do i = 1, nsubevent
    ot_tmp = ot_begin + ot_shift * real(i - 1, kind = fp)

    sigma_lon = evloc_sub(i)%sigma_x / ((r_earth - mean_depth) * sin(pi * 0.5_fp - mean_lat * deg2rad)) * rad2deg
    sigma_lat = evloc_sub(i)%sigma_y / (r_earth - mean_depth) * rad2deg

    write(10, '(7(e15.8, 1x), i0, 1x, f7.1)') &
    &          evloc_sub(i)%lon, sigma_lon, &
    &          evloc_sub(i)%lat, sigma_lat, &
    &          evloc_sub(i)%depth, evloc_sub(i)%sigma_depth, data_residual(i), nsta_use_min(i), ot_tmp

    write(0, '(a, i0, a, f8.1)')  "subevent index = ", i, " ot_tmp = ", ot_tmp
    write(0, '(a, 2(e14.7, 1x))') "longitude and sigma_lon = ", evloc_sub(i)%lon, sigma_lon
    write(0, '(a, 2(e14.7, 1x))') "latitude and sigma_lat = ",  evloc_sub(i)%lat, sigma_lat
    write(0, '(a, 2(e14.7, 1x))') "depth and sigma_depth = ",   evloc_sub(i)%depth, evloc_sub(i)%sigma_depth
  enddo
  close(10)
#else
  write(10, '(a)') "# longitude sigma_lon latitude sigma_lat depth sigma_depth residual_sum nsta"
  do i = 1, nsubevent
    sigma_lon = evloc_sub(i)%sigma_x / ((r_earth - mean_depth) * sin(pi * 0.5_fp - mean_lat * deg2rad)) * rad2deg
    sigma_lat = evloc_sub(i)%sigma_y / (r_earth - mean_depth) * rad2deg
    write(10, '(7(e15.8, 1x), i0)') &
    &          evloc_sub(i)%lon, sigma_lon, &
    &          evloc_sub(i)%lat, sigma_lat, &
    &          evloc_sub(i)%depth, evloc_sub(i)%sigma_depth, data_residual(i), nsta_use_min(i)
    write(0, '(a, i0)')  "subevent index = ", i
    write(0, '(a, 2(e14.7, 1x))') "longitude and sigma_lon = ", evloc_sub(i)%lon, sigma_lon
    write(0, '(a, 2(e14.7, 1x))') "latitude and sigma_lat = ",  evloc_sub(i)%lat, sigma_lat
    write(0, '(a, 2(e14.7, 1x))') "depth and sigma_depth = ",   evloc_sub(i)%depth, evloc_sub(i)%sigma_depth
    write(0, '(a, e14.7)') "xeast = ", evloc_sub(i)%xeast
    write(0, '(a, e14.7)') "ynorth = ", evloc_sub(i)%ynorth
  enddo
  close(10)
#endif

  resultfile_tmp = trim(resultfile) // "_evamp_sub.txt"
  open(unit = 10, file = trim(resultfile_tmp))
  write(10, '(a)') "# masterevent index, amplitude ratio, and error of amplitude ratio"
  do j = 1, nsubevent
#if defined (WIN) || defined (SAC)
    ot_tmp = ot_begin + ot_shift * real(j - 1, kind = fp)
    cfmt = '(x(i0, 1x), f8.1)'
    write(cfmt(2 : 2), '(i1)') nev_master
    write(10, cfmt) (evamp_sub(i, j)%evindex, i = 1, nev_master), ot_tmp
#else
    cfmt = '(x(i0, 1x))'
    write(cfmt(2 : 2), '(i1)') nev_master
    write(10, cfmt) (evamp_sub(i, j)%evindex, i = 1, nev_master)
#endif
    cfmt = '(x(e15.7, 1x))'
    write(cfmt(2 : 2), '(i1)') nev_master
    write(10, cfmt) (evamp_sub(i, j)%sourceamp, i = 1, nev_master)
    write(10, cfmt) (evamp_sub(i, j)%sigma_sourceamp, i = 1, nev_master)
  enddo
  close(10)

  stop
end program AmplitudeSourceLocation_masterevent_ttimecor

