program AmplitudeSourceLocation_masterevent
  !!Relative amplitude source location (master event)
  !!using depth-dependent 1D velocity structure, 3D heterogeneous attenuation structure
  !!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2020
  !!License  : MIT License (https://opensource.org/licenses/MIT)

  use nrtype,               only : fp, dp
  use constants,            only : rad2deg, deg2rad, pi, r_earth
  use rayshooting,          only : rayshooting3D
  use set_velocity_model,   only : set_velocity
  use linear_interpolation, only : linear_interpolation_1d, linear_interpolation_2d, linear_interpolation_3d, &
  &                                block_interpolation_3d
  use greatcircle,          only : greatcircle_dist
  use grdfile_io,           only : read_grdfile_2d
  use read_shmdump,         only : read_shmdump_win
  use m_winch
#ifdef MKL
  use lapack95
#else
  use f95_lapack
#endif

  implicit none

  integer,            parameter :: wavetype = 2          !!1 for P-wave, 2 for S-wave
  integer,            parameter :: nsec_wavebuf = 20
  integer,            parameter :: nsample_max = 200;
  integer,            parameter :: nsec_read = 1 
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

  real(kind = dp),    parameter   :: order = 1.0e+6_dp
  real(kind = fp),    allocatable :: sampling(:), begin(:), ttime(:), conv(:)
  real(kind = dp),    allocatable :: waveform_obs(:, :)
  integer,            allocatable :: ndata_sec(:, :), yr(:), mo(:), dy(:), hh(:), mm(:), ss(:)
  character(len = 4), allocatable :: st_winch(:)
  type(winch__hdr),   allocatable :: chtbl(:)
  character(len = 129)            :: win_chfilename
  character(len = 10)             :: cmpnm

  real(kind = fp),    allocatable :: stlon(:), stlat(:), stdp(:), ttime_cor(:, :)
  character(len = 6), allocatable :: stname(:)
  logical,            allocatable :: use_flag(:)
  character(len = 10)             :: freq_t

  real(kind = fp)                 :: ot_shift, ot_tmp, rms_tw
  character(len = 10)             :: ot_shift_t, rms_tw_t

  !!bandpass filter
  real(kind = dp),    parameter   :: ap = 0.5_dp, as = 5.0_dp
  real(kind = dp),    allocatable :: h(:, :), uv(:, :), gn(:)
  integer,            allocatable :: m(:)
  logical,            allocatable :: filtered(:, :)
  character(len = 6)              :: fl_t, fh_t, fs_t
  real(kind = dp)                 :: fl, fh, fs, c, sampling_tmp
  integer                         :: n

#ifdef DAMPED
  real(kind = fp),    parameter :: damp(4) = [0.0_fp, 0.0_fp, 0.0_fp, 1.0_fp]  !![amp, dx, dy, dz]
#endif

  real(kind = fp),    parameter :: alt_to_depth = -1.0e-3_fp
  real(kind = dp),    parameter :: huge = 1.0e+10_dp

  real(kind = fp)               :: velocity(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                qinv(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                val_1d(1 : 2), val_2d(1 : 2, 1 : 2), val_3d(1 : 2, 1 : 2, 1 : 2), &
  &                                xgrid(1 : 2), ygrid(1 : 2), zgrid(1 : 2), inc_angle_ini_min(0 : nrayshoot), &
  &                                normal_vector(1 : 3)
  real(kind = dp),  allocatable :: topography(:, :), lon_topo(:), lat_topo(:)
  real(kind = fp),  allocatable :: obsamp_master(:), obsamp_sub(:, :), hypodist(:), ray_azinc(:, :), dist_min(:), &
  &                                obsvector(:), obsvector_copy(:), &
  &                                inversion_matrix(:, :), inversion_matrix_copy(:, :), &
  &                                sigma_inv_data(:, :), error_matrix(:, :)
  integer,          allocatable :: ipiv(:)
  real(kind = fp)               :: evlon_master, evlat_master, evdp_master, &
  &                                epdist, epdelta, az_ini, dinc_angle_org, dinc_angle, inc_angle_ini, &
  &                                lon_tmp, lat_tmp, depth_tmp, az_tmp, inc_angle_tmp, dist_tmp, velocity_interpolate, &
  &                                dvdz, lon_new, lat_new, depth_new, az_new, inc_angle_new, ttime_tmp, matrix_const, &
  &                                lon_min, lat_min, depth_min, delta_depth, delta_lon, delta_lat, &
  &                                data_residual, data_variance, sigma_lon, sigma_lat, sigma_depth, sigma_amp, &
  &                                depth_max_tmp, depth_max, data_residual_square
  real(kind = dp)               :: topography_interpolate, dlon_topo, dlat_topo, qinv_interpolate, freq
  integer                       :: nlon_topo, nlat_topo, nsta, lon_index, lat_index, z_index, &
  &                                i, j, k, ii, jj, kk, icount, nsta_use, ios, index_tmp(1), ikey
  character(len = 129)          :: topo_grd, station_param, masterevent_param, subevent_param, resultfile
  character(len = 20)           :: cfmt, nsta_c


  icount = iargc()
  if(icount .ne. 10) then
    write(0, '(a)', advance="no") "usage: shmdump -tq xxx | asl_masterevent_shmdump (topography_grd) "
    write(0, '(a)', advance="no") "(station_param_file) (masterevent_param_file) (win_chfile) (component name) (fl) (fh) (fs) "
    write(0, '(a)')               "(ot_shift) (length of rms time window)"
    error stop
  endif
  call getarg(1, topo_grd)
  call getarg(2, station_param)
  call getarg(3, masterevent_param)
  call getarg(4, win_chfilename)
  call getarg(5, cmpnm)
  call getarg(6, fl_t); read(fl_t, *) fl
  call getarg(7, fh_t); read(fh_t, *) fh
  call getarg(8, fs_t); read(fs_t, *) fs
  call getarg(9, ot_shift_t); read(ot_shift_t, *) ot_shift
  call getarg(10, rms_tw_t); read(rms_tw_t, *) rms_tw

#ifdef DAMPED
  write(0, '(a)') "-DDAMPED set"
  write(0, '(a, 4(f5.2, 1x), a)') "Damped matrix = [ ", (damp(i), i = 1, 4), "]"
#endif

  write(0, '(a, 3(f5.2, 1x))') "Bandpass filter parameter fl, fh, fs (Hz) = ", fl, fh, fs
  freq = (fl + fh) * 0.5_dp

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
  !!allocation related to station parameters
  allocate(stlon(1 : nsta), stlat(1 : nsta), stdp(1 : nsta), stname(1 : nsta), ttime_cor(1 : nsta, 1 : 2), use_flag(1 : nsta))
  
  nsta_use = 0
  do i = 1, nsta
    read(10, *) stlon(i), stlat(i), stdp(i), stname(i), use_flag(i), ttime_cor(i, 1), ttime_cor(i, 2)
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
  allocate(obsamp_master(nsta))
  open(unit = 10, file = masterevent_param)
  read(10, *)
  read(10, *) evlon_master, evlat_master, evdp_master
  read(10, *) (obsamp_master(i), i = 1, nsta)
  close(10)


  !!allocation related to channel table
  allocate(st_winch(nsta), conv(nsta))
  !!read win channel table
  call winch__read_tbl(trim(win_chfilename), chtbl)
  !!find chid from given stname(nsta) and cmpnm
  do i = 1, nsta
    call winch__st2chid(chtbl, stname(i), trim(cmpnm), st_winch(i), ikey)
    if(ikey .eq. 0) then
      write(0, '(5(a, 1x))') "station/comp =", trim(stname(i)), trim(cmpnm), "does not exist in", trim(win_chfilename)
      error stop
    endif
    write(0, '(6a)') "station name = ", trim(stname(i)), " comp = ", trim(cmpnm), " chid = ", st_winch(i)
    conv(i) = chtbl(ikey)%conv
  enddo

  !!bandpass filter
  !write(0, '(a)') "Applying bandpass filter"
  !do j = 1, nsta
  !  call calc_bpf_order(fl, fh, fs, ap, as, real(sampling(j), kind = dp), m, n, c)
  !  allocate(h(4 * m))
  !  call calc_bpf_coef(fl, fh, real(sampling(j), kind = dp), m, n, h, c, gn)
  !  call tandem1(waveform_obs(:, j), waveform_obs(:, j), npts(j), h, m, 1)
  !  waveform_obs(1 : npts(j), j) = waveform_obs(1 : npts(j), j) * gn
  !  deallocate(h)
  !enddo

  !!set velocity/attenuation structure
  call set_velocity(z_str_min, dz_str, velocity, qinv)

  !!calculate ray length, pulse width, unit vector of ray incident
  write(0, '(a)') "calculate ray length, pulse width, and ray incident vector for master event"

  !!check whether the depth of master event location is lower than the topo
  lon_index = int((evlon_master - lon_topo(1)) / dlon_topo) + 1
  lat_index = int((evlat_master - lat_topo(1)) / dlat_topo) + 1
  xgrid(1 : 2) = [lon_topo(lon_index), lon_topo(lon_index + 1)]
  ygrid(1 : 2) = [lat_topo(lat_index), lat_topo(lat_index + 1)]
  val_2d(1 : 2, 1 : 2) = topography(lon_index : lon_index + 1, lat_index : lat_index + 1)
  call linear_interpolation_2d(evlon_master, evlat_master, xgrid, ygrid, val_2d, topography_interpolate)
  if(evdp_master .lt. topography_interpolate) then
    write(0, '(a, f5.2, a)') "Depth of master event ", evdp_master, " is higher than the altitude there."
    error stop
  endif

  !!allocation related to ray tracing
  allocate(hypodist(1 : nsta), ray_azinc(1 : 2, 1 : nsta), dist_min(1 : nsta), ttime(1 : nsta))
  
  !!calculate traveltime to station        
  station_loop: do jj = 1, nsta
    !!calculate azimuth and hypocentral distance
    call greatcircle_dist(evlat_master, evlon_master, stlat(jj), stlon(jj), &
    &                     distance = epdist, azimuth = az_ini, delta_out = epdelta)
    lon_index = int((stlon(jj) - lon_str_w) / dlon_str) + 1
    lat_index = int((stlat(jj) - lat_str_s) / dlat_str) + 1
    z_index   = int((stdp(jj) - z_str_min) / dz_str) + 1
    !print *, lon_sta(jj), lon_w + real(lon_index - 1) * dlon
    !print *, lat_sta(jj), lat_s + real(lat_index - 1) * dlat
    hypodist(jj) = sqrt((r_earth - evdp_master) ** 2 + (r_earth - stdp(jj)) ** 2 &
    &            - 2.0_fp * (r_earth - evdp_master) * (r_earth - stdp(jj)) * cos(epdelta))
    ray_azinc(1, jj) = az_ini

#ifdef V_CONST
    !!homogeneous structure: ray incident angle is calculated using cosine function (assuming cartesian coordinate)
    ray_azinc(2, jj) = acos(epdist / hypodist(jj)) + pi / 2.0_fp

#else
          
    !!do ray shooting
    dist_min(jj) = real(huge, kind = fp)
    incangle_loop2: do kk = 1, nrayshoot
      if(kk .eq. 1) then
        dinc_angle_org = pi / 2.0_fp
      else
        dinc_angle_org = dinc_angle
      endif
      dinc_angle = 2.0_fp * dinc_angle_org / real(ninc_angle, kind = fp)
      inc_angle_ini_min(0) = dinc_angle_org
      !print *, "dinc_angle = ", dinc_angle * rad2deg, inc_angle_ini_min(kk - 1) * rad2deg

      depth_max_tmp = evdp_master
      incangle_loop: do ii = 1, ninc_angle
        inc_angle_ini = (inc_angle_ini_min(kk - 1) - dinc_angle_org) + real(ii, kind = fp) * dinc_angle
        !print '(2(i0, 1x), a, e15.7)', ii, kk, "inc_angle_ini = ", inc_angle_ini * rad2deg
               
        lon_tmp = evlon_master
        lat_tmp = evlat_master
        depth_tmp = evdp_master
        az_tmp = az_ini
        inc_angle_tmp = inc_angle_ini

        !!loop until ray arrives at surface/boundary
        ttime_tmp = 0.0_fp
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

          !!exit if ray approaches to the boundary
          if(lon_index .lt. 1 .or. lon_index .gt. nlon_str - 1      &
          &  .or. lat_index .lt. 1 .or. lat_index .gt. nlat_str - 1 &
          &  .or. z_index .lt. 1 .or. z_index .gt. nz_str - 1) then
            exit shooting_loop
          endif

          xgrid(1) = lon_str_w + real(lon_index - 1, kind = fp) * dlon_str; xgrid(2) = xgrid(1) + dlon_str
          ygrid(1) = lat_str_s + real(lat_index - 1, kind = fp) * dlat_str; ygrid(2) = ygrid(1) + dlat_str
          zgrid(1) = z_str_min + real(z_index - 1, kind = fp) * dz_str;     zgrid(2) = zgrid(1) + dz_str

          !!calculate distance between ray and station
          call greatcircle_dist(lat_tmp, lon_tmp, stlat(jj), stlon(jj), delta_out = epdelta)
          dist_tmp = sqrt((r_earth - depth_tmp) ** 2 + (r_earth - stdp(jj)) ** 2 &
          &        - 2.0_fp * (r_earth - depth_tmp) * (r_earth - stdp(jj)) * cos(epdelta))
          !print *, real(ii, kind = fp) * dinc_angle, lat_tmp, lon_tmp, depth_tmp
          !print *, jj, lat_sta(jj), lon_sta(jj), z_sta(jj), dist_tmp

          if(dist_tmp .lt. dist_min(jj)) then
            dist_min(jj) = dist_tmp
            lon_min = lon_tmp
            lat_min = lat_tmp
            depth_min = depth_tmp
            inc_angle_ini_min(kk) = inc_angle_ini
            ray_azinc(2, jj) = inc_angle_ini
            depth_max = depth_max_tmp
            ttime(jj) = ttime_tmp
          endif

          !!shooting the ray
          val_1d(1 : 2) = velocity(lon_index, lat_index, z_index : z_index + 1, wavetype)
          call linear_interpolation_1d(depth_tmp, zgrid, val_1d, velocity_interpolate)
          dvdz = (val_1d(2) - val_1d(1)) / dz_str
          call rayshooting3D(lon_tmp, lat_tmp, depth_tmp, az_tmp, inc_angle_tmp, time_step, velocity_interpolate, &
          &                  dvdlon, dvdlat, dvdz, lon_new, lat_new, depth_new, az_new, inc_angle_new)

          lon_tmp = lon_new
          lat_tmp = lat_new
          depth_tmp = depth_new
          az_tmp = az_new
          inc_angle_tmp = inc_angle_new
          ttime_tmp = ttime_tmp + time_step
          if(depth_tmp .ge. depth_max_tmp) depth_max_tmp = depth_tmp
        enddo shooting_loop

      enddo incangle_loop
    enddo incangle_loop2
    write(0, '(a, 4(f8.4, 1x))') "masterevent lon, lat, depth, az_ini = ", &
    &                            evlon_master, evlat_master, evdp_master, az_ini * rad2deg
    write(0, '(a, 3(f8.4, 1x))') "station lon, lat, depth = ", stlon(jj), stlat(jj), stdp(jj)
    write(0, '(a, 4(f8.4, 1x))') "rayshoot lon, lat, depth, dist = ", lon_min, lat_min, depth_min, dist_min(jj)
    write(0, '(a, 2(f8.4, 1x))') "inc_angle (deg), ray turning depth (km) = ", &
    &                             inc_angle_ini_min(nrayshoot) * rad2deg, depth_max
    write(0, '(a, f5.2)') "traveltime (s) = ", ttime(jj)
#endif
  enddo station_loop

  ttime(1 : nsta) = ttime(1 : nsta) + ttime_cor(1 : nsta, wavetype)
#ifdef WITHOUT_TTIME
  ttime(1 : nsta) = 0.0_fp
#endif

  !!velocity and Qinv at master event location
  lon_index = int((evlon_master - lon_str_w) / dlon_str) + 1
  lat_index = int((evlat_master - lat_str_s) / dlat_str) + 1
  z_index   = int((evdp_master - z_str_min) / dz_str) + 1
  xgrid(1) = lon_str_w + real(lon_index - 1, kind = fp) * dlon_str; xgrid(2) = xgrid(1) + dlon_str
  ygrid(1) = lat_str_s + real(lat_index - 1, kind = fp) * dlat_str; ygrid(2) = ygrid(1) + dlat_str
  zgrid(1) = z_str_min + real(z_index - 1, kind = fp) * dz_str;     zgrid(2) = zgrid(1) + dz_str
  val_3d(1 : 2, 1 : 2, 1 : 2) &
  &  = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1, wavetype)
  call linear_interpolation_3d(evlon_master, evlat_master, evdp_master, xgrid, ygrid, zgrid, val_3d, velocity_interpolate)
  val_3d(1 : 2, 1 : 2, 1 : 2) &
  &  = qinv(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1, wavetype)
  call block_interpolation_3d(evlon_master, evlat_master, evdp_master, xgrid, ygrid, zgrid, val_3d, qinv_interpolate)

  !!Set up the observation vector and inversion matrix
#ifdef DAMPED
  allocate(inversion_matrix(1 : nsta_use + 4, 1 : 4), obsvector(1 : nsta_use + 4))
  inversion_matrix(1 : nsta_use + 4, 1 : 4) = 0.0_fp
  obsvector(1 : nsta_use + 4) = 0.0_fp
#else
  allocate(inversion_matrix(1 : nsta_use, 1 : 4), obsvector(1 : nsta_use))
  inversion_matrix(1 : nsta_use, 1 : 4) = 0.0_fp
  obsvector(1 : nsta_use) = 0.0_fp
#endif
  allocate(inversion_matrix_copy(1 : nsta_use, 1 : 4), obsvector_copy(1 : nsta_use))

  allocate(yr(nsec_wavebuf), mo(nsec_wavebuf), dy(nsec_wavebuf), hh(nsec_wavebuf), mm(nsec_wavebuf), ss(nsec_wavebuf))
  allocate(waveform_obs(nsec_wavebuf * nsample_max, nsta), ndata_sec(nsec_wavebuf, nsta))
  allocate(m(nsta), filtered(nsec_wavebuf, nsta))
  filtered(1 : nsec_wavebuf, 1 : nsta) = .true.
  waveform_obs(1 : nsec_wavebuf * nsample_max, 1 : nsta) = 0.0_dp
  ndata_sec(1 : nsec_wavebuf, 1 : nsta) = 0
  do 
    !!read waveform from stdin
    do i = 1, nsec_read
      !call read_shmdump_win(st_winch, conv, yr, mo, dy, hh, mm, ss, waveform_obs, ndata_sec)
    enddo
    j = sum(ndata_sec(1 : nsec_wavebuf, 1))
    do i = 1, j
      print *, waveform_obs(i, 1)
    enddo
    if(allocated(h) .eqv. .false.) then
      do i = 1, nsta
        sampling_tmp = 1.0_dp / real(ndata_sec(1, i), kind = dp)
        call calc_bpf_order(fl, fh, fs, ap, as, sampling_tmp, m(i), n, c)
      enddo
      allocate(h(1 : 4 * maxval(m), 1 : nsta), uv(1 : 4 * maxval(m), 1 : nsta), gn(1 : nsta))
      h(1 : 4 * maxval(m), 1 : nsta) = 0.0_dp
      uv(1 : 4 * maxval(m), 1 : nsta) = 0.0_dp
      gn(1 : nsta) = 0.0_dp
      do i = 1, nsta
        sampling_tmp = 1.0_dp / real(ndata_sec(1, i))
        call calc_bpf_coef(fl, fh, sampling_tmp, m(i), n, h(1 : 4 * m(i), i), c, gn(i))
      enddo
    endif
    
    !do j = 1, nsta
    !  ndata_tmp = 0
    !  sec_loop: do i = 1, nsec_read
    !    if(filtered(i, j) .eqv. .true.) cycle
    !    call tandem2(waveform_obs

        


!    icount = 1
!    do i = 1, nsta
!      if(use_flag(i) .eqv. .false.) cycle
!      obsvector(icount) = log(obsamp_sub(i, j) / obsamp_master(i))
!      normal_vector(1 : 3) = [sin(ray_azinc(2, i)) * cos(ray_azinc(1, i)), &
!      &                       sin(ray_azinc(2, i)) * sin(ray_azinc(1, i)), &
!      &                       cos(ray_azinc(2, i))]
!      inversion_matrix(icount, 1) = 1.0_fp
!      matrix_const = pi * freq * qinv_interpolate / velocity_interpolate + 1.0_fp / hypodist(i)
!      do ii = 1, 3
!        inversion_matrix(icount, 1 + ii) = (-1.0_fp) * matrix_const * normal_vector(ii)
!      enddo
!      icount = icount + 1
!    enddo
!#ifdef DAMPED
!    do i = 1, 4
!      inversion_matrix(nsta_use + i, i) = damp(i)
!    enddo
!#endif
!
!    !!copy observation vector and inversion matrix
!    obsvector_copy(1 : nsta_use) = obsvector(1 : nsta_use)
!    inversion_matrix_copy(1 : nsta_use, 1 : 4) = inversion_matrix(1 : nsta_use, 1 : 4)
!
!    !!calculate least-squares solution
!#ifdef MKL
!    call gels(inversion_matrix, obsvector)
!#else
!    call la_gels(inversion_matrix, obsvector)
!#endif
!
!    !!calculate mean data residual
!    data_residual = 0.0_fp
!    data_residual_square = 0.0_fp
!    icount = 0
!    do i = 1, nsta_use
!      if(dot_product(inversion_matrix_copy(i, 1 : 4), obsvector(1 : 4)) .eq. 0.0_fp) cycle
!      icount = icount + 1
!      data_residual = data_residual &
!      &             + (obsvector_copy(i) - dot_product(inversion_matrix_copy(i, 1 : 4), obsvector(1 : 4)))
!      data_residual_square = data_residual_square &
!      &             + (obsvector_copy(i) - dot_product(inversion_matrix_copy(i, 1 : 4), obsvector(1 : 4))) ** 2
!    enddo
!    data_residual_square = data_residual_square / real(icount, kind = fp)
!    data_residual = data_residual / real(icount, kind = fp)
!
!    !!estimate errors of inverted model parameters
!    allocate(error_matrix(1 : 4, 1 : 4), sigma_inv_data(1 : nsta_use, 1 : nsta_use))
!    error_matrix(1 : 4, 1 : 4) = 0.0_fp
!    sigma_inv_data(1 : nsta_use, 1 : nsta_use) = 0.0_fp  
!    !!calculate variance
!    data_variance = 0.0_fp
!    do i = 1, nsta_use
!      data_variance = data_variance + (data_residual &
!      &              - (obsvector_copy(i) -  dot_product(inversion_matrix_copy(i, 1 : 4), obsvector(1 : 4)))) ** 2
!    enddo
!    if(nsta_use - 1 .ne. 0) data_variance = data_variance / real(nsta_use - 1, kind = fp)
!    do i = 1, nsta_use
!      sigma_inv_data(i, i) = 1.0_fp / data_variance
!    enddo
!    error_matrix = matmul(matmul(transpose(inversion_matrix_copy), sigma_inv_data), inversion_matrix_copy)
!    allocate(ipiv(1 : size(error_matrix, 1)))
!#ifdef MKL
!    call getrf(error_matrix, ipiv)
!    call getri(error_matrix, ipiv)
!#else
!    call la_getrf(error_matrix, ipiv)
!    call la_getri(error_matrix, ipiv)
!#endif
!    deallocate(ipiv, sigma_inv_data)
!
!  !!output result into stdout
!    delta_lat = (obsvector(4 * (i - 1) + 2) / (r_earth - evdp_master)) * rad2deg
!    delta_lon = (obsvector(4 * (i - 1) + 3) / ((r_earth - evdp_master) * sin(pi / 2.0_fp - evlat_master * deg2rad))) * rad2deg
!    delta_depth = obsvector(4 * (i - 1) + 4)
!
!    sigma_amp = exp(obsvector(4 * (i - 1) + 1)) * sqrt(error_matrix(4 * (i - 1) + 1, 4 * (i - 1) + 1)) * 2.0_fp
!    sigma_lat = sqrt(error_matrix(4 * (i - 1) + 2, 4 * (i - 1) + 2)) * rad2deg / (r_earth - evdp_master) * 2.0_fp
!    sigma_lon = sqrt(error_matrix(4 * (i - 1) + 3, 4 * (i - 1) + 3)) &
!    &         * sin(pi / 2.0_fp - evlat_master * deg2rad) * rad2deg / (r_earth - evdp_master) * 2.0_fp
!    sigma_depth = sqrt(error_matrix(4 * (i - 1) + 4, 4 * (i - 1) + 4)) * 2.0_fp
!
!    write(6, '(9(e15.8, 1x), f7.1)') &
!    &          exp(obsvector(4 * (i - 1) + 1)), sigma_amp, &
!    &          evlon_master - delta_lon, sigma_lon, &
!    &          evlat_master - delta_lat, sigma_lat, &
!    &          evdp_master - delta_depth, sigma_depth, data_residual_square, ot_tmp

  enddo

  stop
end program AmplitudeSourceLocation_masterevent

