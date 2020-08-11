program TraveltimeSourceLocation_masterevent
  !!Relative traveltime source location (master event)
  !!using depth-dependent 1D velocity structure, 3D heterogeneous attenuation structure
  !!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2020
  !!License  : MIT License (https://opensource.org/licenses/MIT)

  !!velocity(1 : nlon, 1 : nlat, 1 : nz) defines wave type; if P-wave,
  !!observation should be P-wave travel time, while S-wave travel time when
  !!S-wave velocity is given in velocity(1 : nlon, 1 : nlat, 1 : nz).
  !!In this source, I do not check whether P- or S-wave velocity is given

  use nrtype,               only : fp, dp
  use constants,            only : rad2deg, deg2rad, pi, r_earth
  use rayshooting,          only : rayshooting3D
  use set_velocity_model,   only : set_velocity
  use linear_interpolation, only : linear_interpolation_1d, linear_interpolation_2d, linear_interpolation_3d, &
  &                                block_interpolation_3d
  use greatcircle,          only : greatcircle_dist
  use grdfile_io,           only : read_grdfile_2d
  !$ use omp_lib
#ifdef MKL
  use lapack95
#else
  use f95_lapack
#endif

  implicit none

  !!Range for velocity and attenuation structure
  real(kind = fp),    parameter :: lon_str_w = 143.5_fp, lon_str_e = 144.1_fp
  real(kind = fp),    parameter :: lat_str_s = 43.0_fp, lat_str_n = 43.50_fp
  real(kind = fp),    parameter :: z_str_min = -1.5_fp, z_str_max = 10.0_fp
  real(kind = fp),    parameter :: dlon_str = 0.01_fp, dlat_str = 0.01_fp, dz_str = 0.1_fp
  integer,            parameter :: nlon_str = int((lon_str_e - lon_str_w) / dlon_str) + 2
  integer,            parameter :: nlat_str = int((lat_str_n - lat_str_s) / dlat_str) + 2
  integer,            parameter :: nz_str = int((z_str_max - z_str_min) / dz_str) + 2
  !!Ray shooting
  real(kind = fp),    parameter :: dvdlon = 0.0_fp, dvdlat = 0.0_fp         !!assume 1D structure
  integer,            parameter :: ninc_angle = 180                         !!grid search in incident angle
  integer,            parameter :: nrayshoot = 2                            !!number of grid search
  real(kind = fp),    parameter :: time_step = 0.01_fp
  real(kind = fp),    parameter :: rayshoot_dist_thr = 0.05_fp
  !!Use station

  real(kind = fp),    parameter :: alt_to_depth = -1.0e-3_fp
  real(kind = dp),    parameter :: freq = 7.5_dp
  real(kind = dp),    parameter :: huge = 1.0e+10_dp

  real(kind = fp)               :: velocity(1 : nlon_str, 1 : nlat_str, 1 : nz_str), &
  &                                qinv(1 : nlon_str, 1 : nlat_str, 1 : nz_str), &
  &                                val_1d(1 : 2), val_2d(1 : 2, 1 : 2), val_3d(1 : 2, 1 : 2, 1 : 2), &
  &                                xgrid(1 : 2), ygrid(1 : 2), zgrid(1 : 2), inc_angle_ini_min(0 : nrayshoot), &
  &                                normal_vector(1 : 3)
  real(kind = dp),  allocatable :: topography(:, :), lon_topo(:), lat_topo(:)
  real(kind = fp),  allocatable :: stlon(:), stlat(:), stdp(:), traveltime_master(:), traveltime_sub(:, :), &
  &                                hypodist(:), ray_azinc(:, :), dist_min(:), &
  &                                obsvector(:), obsvector_copy(:), &
  &                                inversion_matrix(:, :), inversion_matrix_copy(:, :), &
  &                                sigma_inv_data(:, :), error_matrix(:, :)
  integer,          allocatable :: ipiv(:)
  
  real(kind = fp)               :: evlon_master, evlat_master, evdp_master, &
  &                                epdist, epdelta, az_ini, dinc_angle_org, dinc_angle, inc_angle_ini, &
  &                                lon_tmp, lat_tmp, depth_tmp, az_tmp, inc_angle_tmp, dist_tmp, velocity_interpolate, &
  &                                dvdz, lon_new, lat_new, depth_new, az_new, inc_angle_new, matrix_const, &
  &                                lon_min, lat_min, depth_min, delta_depth, delta_lon, delta_lat, &
  &                                data_residual, data_variance, sigma_lon, sigma_lat, sigma_depth, sigma_amp

  real(kind = dp)               :: topography_interpolate, dlon_topo, dlat_topo, qinv_interpolate

  integer                       :: nlon_topo, nlat_topo, nsta, nsubevent, lon_index, lat_index, z_index, &
  &                                i, j, ii, jj, kk, icount

  character(len = 129)          :: topo_grd, station_param, masterevent_param, subevent_param, resultfile

  !!OpenMP variable
  !$ integer                    :: omp_thread

  icount = iargc()
  if(icount .ne. 5) then
    write(0, '(a)', advance="no") "usage: ./asl_masterevent "
    write(0, '(a)', advance="no") "(topography_grd) (station_param_file) (masterevent_param_file) (subevent_param_file) "
    write(0, '(a)')               "(result_file)"
    error stop
  endif
  call getarg(1, topo_grd)
  call getarg(2, station_param)
  call getarg(3, masterevent_param)
  call getarg(4, subevent_param)
  call getarg(5, resultfile)

  !!read topography file (netcdf grd format)
  call read_grdfile_2d(topo_grd, lon_topo, lat_topo, topography)
  nlon_topo = ubound(lon_topo, 1)
  nlat_topo = ubound(lat_topo, 1)
  dlon_topo = lon_topo(2) - lon_topo(1)
  dlat_topo = lat_topo(2) - lat_topo(1)
  topography(1 : nlon_topo, 1 : nlat_topo) = topography(1 : nlon_topo, 1 : nlat_topo) * alt_to_depth

  !!read station parameter
  open(unit = 10, file = station_param)
  read(10, *) nsta
  if(nsta .lt. 4) then
    close(10)
    write(0, '(a)') "Number of station nsta should be larger than 4"
    error stop
  endif
  allocate(stlon(nsta), stlat(nsta), stdp(nsta))
  do i = 1, nsta
    read(10, *) stlon(i), stlat(i), stdp(i)
  enddo
  close(10)
  !!read masterevent parameter
  allocate(traveltime_master(nsta))
  open(unit = 10, file = masterevent_param)
  read(10, *)
  read(10, *) evlon_master, evlat_master, evdp_master
  read(10, *) (traveltime_master(i), i = 1, nsta)
  close(10)
  !!read subevent paramter
  open(unit = 10, file = subevent_param)
  read(10, *)
  read(10, *) nsubevent
  allocate(traveltime_sub(1 : nsta, 1 : nsubevent))
  do j = 1, nsubevent
    read(10, *) (traveltime_sub(i, j), i = 1, nsta)
  enddo
  close(10)

  !!set velocity/attenuation structure
  call set_velocity(z_str_min, dz_str, velocity, qinv)

  !!calculate ray length, pulse width, unit vector of ray incident
  write(0, '(a)') "calculate ray length and ray incident vector for master event"

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

  allocate(hypodist(1 : nsta), ray_azinc(1 : 2, 1 : nsta), dist_min(1 : nsta))
  
  !$omp parallel default(none), &
  !$omp&         shared(nsta, evlon_master, evlat_master, evdp_master, stlon, stlat, stdp, hypodist, ray_azinc, &
  !$omp&                dist_min, lon_topo, lat_topo, dlon_topo, dlat_topo, topography, nlon_topo, nlat_topo, velocity, qinv), &
  !$omp&         private(epdist, az_ini, epdelta, lon_index, lat_index, z_index, dinc_angle_org, dinc_angle, &
  !$omp&                 inc_angle_ini_min, inc_angle_ini, lon_tmp, lat_tmp, depth_tmp, az_tmp, inc_angle_tmp, &
  !$omp&                 xgrid, ygrid, zgrid, val_1d, val_2d, val_3d, topography_interpolate, &
  !$omp&                 dist_tmp, lon_min, lat_min, depth_min, velocity_interpolate, dvdz, &
  !$omp&                 lon_new, lat_new, depth_new, az_new, inc_angle_new, ii, kk, omp_thread)

  !$ omp_thread = omp_get_thread_num()

  !!calculate traveltime to station        
  !$omp do
  station_loop: do jj = 1, nsta
    !!calculate azimuth and hypocentral distance
    call greatcircle_dist(evlat_master, evlon_master, stlat(jj), stlon(jj), &
    &                     distance = epdist, azimuth = az_ini, delta_out = epdelta)
    !lon_index = int((stlon(jj) - lon_w) / dlon) + 1
    !lat_index = int((stlat(jj) - lat_s) / dlat) + 1
    !z_index   = int((stdp(jj) - z_min) / dz) + 1
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

      incangle_loop: do ii = 1, ninc_angle
        inc_angle_ini = (inc_angle_ini_min(kk - 1) - dinc_angle_org) + real(ii, kind = fp) * dinc_angle
        !print '(2(i0, 1x), a, e15.7)', ii, kk, "inc_angle_ini = ", inc_angle_ini * rad2deg
               
        lon_tmp = evlon_master
        lat_tmp = evlat_master
        depth_tmp = evdp_master
        az_tmp = az_ini
        inc_angle_tmp = inc_angle_ini

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
            !print '(a, 4(f8.4, 1x))', "rayshoot_tmp lon, lat, depth, dist_min = ", lon_min, lat_min, depth_min, &
            !&                                                                      dist_min
          endif

          !!shooting the ray
          val_1d(1 : 2) = velocity(lon_index, lat_index, z_index : z_index + 1)
          call linear_interpolation_1d(depth_tmp, zgrid, val_1d, velocity_interpolate)
          dvdz = (val_1d(2) - val_1d(1)) / dz_str
          call rayshooting3D(lon_tmp, lat_tmp, depth_tmp, az_tmp, inc_angle_tmp, time_step, velocity_interpolate, &
          &                  dvdlon, dvdlat, dvdz, lon_new, lat_new, depth_new, az_new, inc_angle_new)

          lon_tmp = lon_new
          lat_tmp = lat_new
          depth_tmp = depth_new
          az_tmp = az_new
          inc_angle_tmp = inc_angle_new
        enddo shooting_loop

      enddo incangle_loop
    enddo incangle_loop2
    print '(a, 4(f8.4, 1x))', "masterevent lon, lat, depth, az_ini = ", &
    &                          evlon_master, evlat_master, evdp_master, az_ini * rad2deg
    print '(a, 3(f8.4, 1x))', "station lon, lat, depth = ", stlon(jj), stlat(jj), stdp(jj)
    print '(a, 4(f8.4, 1x))', "rayshoot lon, lat, depth, dist = ", lon_min, lat_min, depth_min, dist_min(jj)
#endif
  enddo station_loop
  !$omp end do
  !$omp end parallel

  !!velocity and Qinv at master event location
  lon_index = int((evlon_master - lon_str_w) / dlon_str) + 1
  lat_index = int((evlat_master - lat_str_s) / dlat_str) + 1
  z_index   = int((evdp_master - z_str_min) / dz_str) + 1
  xgrid(1) = lon_str_w + real(lon_index - 1, kind = fp) * dlon_str; xgrid(2) = xgrid(1) + dlon_str
  ygrid(1) = lat_str_s + real(lat_index - 1, kind = fp) * dlat_str; ygrid(2) = ygrid(1) + dlat_str
  zgrid(1) = z_str_min + real(z_index - 1, kind = fp) * dz_str;     zgrid(2) = zgrid(1) + dz_str
  val_3d(1 : 2, 1 : 2, 1 : 2) = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
  call linear_interpolation_3d(evlon_master, evlat_master, evdp_master, xgrid, ygrid, zgrid, val_3d, velocity_interpolate)

  !!Set up the observation vector and inversion matrix
  allocate(obsvector(1 : nsta * nsubevent), obsvector_copy(1 : nsta * nsubevent))
  allocate(inversion_matrix(1 : nsta * nsubevent, 1 : 4 * nsubevent), &
  &        inversion_matrix_copy(1 : nsta * nsubevent, 1 : 4 * nsubevent))
  inversion_matrix(1 : nsta * nsubevent, 1 : 4 * nsubevent) = 0.0_fp
  do j = 1, nsubevent
    do i = 1, nsta
      obsvector(nsta * (j - 1) + i) = traveltime_master(i) - traveltime_sub(i, j)
      normal_vector(1 : 3) = [sin(ray_azinc(2, i)) * cos(ray_azinc(1, i)), &
      &                       sin(ray_azinc(2, i)) * sin(ray_azinc(1, i)), &
      &                       cos(ray_azinc(2, i))]
      matrix_const = 1.0_fp / velocity_interpolate
      do ii = 1, 3
        inversion_matrix(nsta * (j - 1) + i, 4 * (j - 1) + ii) = matrix_const * normal_vector(ii)
      enddo
      inversion_matrix(nsta * (j - 1) + i, 4 * (j - 1) + 4) = 1.0_fp
    enddo
  enddo

  !!copy observation vector and inversion matrix
  obsvector_copy(1 : nsta * nsubevent) = obsvector(1 : nsta * nsubevent)
  inversion_matrix_copy(1 : nsta * nsubevent, 1 : 4 * nsubevent) &
  &  = inversion_matrix(1 : nsta * nsubevent, 1 : 4 * nsubevent)

  !!calculate least-squares solution
#ifdef MKL
  call gels(inversion_matrix_copy, obsvector)
#else
  call la_gels(inversion_matrix_copy, obsvector)
#endif

  !!calculate mean data residual
  data_residual = 0.0_fp
  do i = 1, nsta * nsubevent
    data_residual = data_residual &
    &             + (obsvector_copy(i) - dot_product(inversion_matrix(i, 1 : 4 * nsubevent), obsvector(1 : 4 * nsubevent)))
  enddo
  data_residual = data_residual / real(nsta * nsubevent, kind = fp)
  !!calculate variance
  data_variance = 0.0_fp
  do i = 1, nsta * nsubevent
    data_variance = data_variance + (data_residual &
    &             - (obsvector_copy(i) - dot_product(inversion_matrix(i, 1 : 4 * nsubevent), obsvector(1 : 4 * nsubevent)))) ** 2
  enddo
  data_variance = data_variance / real(nsta * nsubevent - 1)

  !!estimate error of inverted model parameters
  allocate(sigma_inv_data(1 : nsta * nsubevent, 1 : nsta * nsubevent))
  allocate(error_matrix(1 : nsta * nsubevent, 1 : nsta * nsubevent))
  sigma_inv_data(1 : nsta * nsubevent, 1 : nsta * nsubevent) = 0.0_fp
  do i = 1, nsta * nsubevent
    sigma_inv_data(i, i) = 1.0_fp / data_variance
    !sigma_inv_data(i, i) = &
    !&  abs(obsvector_copy(i) - dot_product(inversion_matrix(i, 1 : 4 * nsubevent), obsvector(1 : 4 * nsubevent)))
    !sigma_inv_data(i, i) = 1.0_fp / sigma_inv_data(i, i)
  enddo
  error_matrix = matmul(matmul(transpose(inversion_matrix), sigma_inv_data), inversion_matrix)
  allocate(ipiv(1 : size(error_matrix, 1)))
#ifdef MKL
  call getrf(error_matrix, ipiv)
  call getri(error_matrix, ipiv)
#else
  call la_getrf(error_matrix, ipiv)
  call la_getri(error_matrix, ipiv)
#endif
  

  !!output result
  open(unit = 10, file = trim(resultfile))
  write(10, '(a)') "# OTdiff sigma_OTdiff longitude sigma_lon latitude sigma_lat depth sigma_depth"
  do i = 1, nsubevent
    delta_lat = (obsvector(4 * (i - 1) + 1) / (r_earth - evdp_master)) * rad2deg
    delta_lon = (obsvector(4 * (i - 1) + 2) / ((r_earth - evdp_master) * sin(pi / 2.0_fp - evlat_master * deg2rad))) * rad2deg
    delta_depth = obsvector(4 * (i - 1) + 3)
    !print *, obsvector(4 * (i - 1) + 1), delta_lon, delta_lat, delta_depth
    sigma_lat = sqrt(error_matrix(4 * (i - 1) + 1, 4 * (i - 1) + 1)) * rad2deg / (r_earth - evdp_master) * 2.0_fp
    sigma_lon = sqrt(error_matrix(4 * (i - 1) + 2, 4 * (i - 1) + 2)) &
    &         * sin(pi / 2.0_fp - evlat_master * deg2rad) * rad2deg / (r_earth - evdp_master) * 2.0_fp
    sigma_depth = sqrt(error_matrix(4 * (i - 1) + 3, 4 * (i - 1) + 3)) * 2.0_fp

    write(10, '(8(e14.7, 1x))') &
    &          obsvector(4 * (i - 1) + 4), sqrt(error_matrix(4 * (i - 1) + 4, 4 * (i - 1) + 4)) * 2.0_fp, 
    &          evlon_master + delta_lon, sigma_lon, &
    &          evlat_master + delta_lat, sigma_lat, &
    &          evdp_master + delta_depth, sigma_depth


    write(0, '(a, i0)')           "subevent index = ", i
    write(0, '(a, 2(e14.7, 1x))') "longitude and sigma_lon = ", evlon_master + delta_lon, sigma_lon
    write(0, '(a, 2(e14.7, 1x))') "latitude and sigma_lat = ", evlat_master + delta_lat, sigma_lat
    write(0, '(a, 2(e14.7, 1x))') "depth and sigma_depth = ", evdp_master + delta_depth, sigma_depth
    write(0, '(a, 2(e14.7, 1x))') "time_diff and sigma_time_diff = ", &
    &                             obsvector(4 * (i - 1) + 4), &
    &                             sqrt(error_matrix(4 * (i - 1) + 4, 4 * (i - 1) + 4)) * 2.0_fp
  enddo
  close(10)
    

  stop
end program TraveltimeSourceLocation_masterevent
