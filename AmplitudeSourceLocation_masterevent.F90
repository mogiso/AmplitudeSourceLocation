program AmplitudeSourceLocation_masterevent
  !!Relative amplitude source location (master event)
  !!using depth-dependent 1D velocity structure, 3D heterogeneous attenuation structure
  !!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2020
  !!License  : MIT License (https://opensource.org/licenses/MIT)

  use nrtype,               only : fp, sp, dp
  use constants,            only : rad2deg, deg2rad, pi, r_earth
  use rayshooting,          only : rayshooting3D
  use set_velocity_model,   only : set_velocity
  use linear_interpolation, only : linear_interpolation_1d, linear_interpolation_2d, block_interpolation_3d
  use greatcircle,          only : greatcircle_dist
  use itoa,                 only : int_to_char
  use grdfile_io,           only : read_grdfile_2d
  !$ use omp_lib

  implicit none

  !!Range for velocity and attenuation structure
  real(kind = fp),    parameter :: lon_w = 143.90_fp, lon_e = 144.05_fp
  real(kind = fp),    parameter :: lat_s = 43.35_fp, lat_n = 43.410_fp
  real(kind = fp),    parameter :: z_min = -1.5_fp, z_max = 5.0_fp
  !!Ray shooting
  real(kind = fp),    parameter :: dvdlon = 0.0_fp, dvdlat = 0.0_fp         !!assume 1D structure
  integer,            parameter :: ninc_angle = 200                         !!grid search in incident angle
  integer,            parameter :: nrayshoot = 2                            !!number of grid search
  real(kind = fp),    parameter :: time_step = 0.01_fp
  real(kind = fp),    parameter :: rayshoot_dist_thr = 0.05_fp
  !!Use station
  integer,            parameter :: nsta = 4

  real(kind = fp),    parameter :: alt_to_depth = -1.0e-3_fp
  real(kind = dp),    parameter :: freq = 7.5_dp

  real(kind = fp)               :: velocity(1 : nlon, 1 : nlat, 1 : nz), qinv(1 : nlon, 1 : nlat, 1 : nz), &
  &                                sampling(1 : nsta), begin(1 : nsta), &
  &                                val_1d(1 : 2), val_2d(1 : 2, 1 : 2), val_3d(1 : 2, 1 : 2, 1 : 2), &
  &                                xgrid(1 : 2), ygrid(1 : 2), zgrid(1 : 2), &
  &                                lon_sta(1 : nsta), lat_sta(1 : nsta), z_sta(1 : nsta), &
  &                                width_min(1 : nsta, 1 : nlon, 1 : nlat, 1 : nz), &
  &                                ttime_min(1 : nsta, 1 : nlon, 1 : nlat, 1 : nz), &
  &                                hypodist(1 : nsta, 1 : nlon, 1 : nlat, 1 : nz), inc_angle_ini_min(0 : nrayshoot)
  real(kind = dp), allocatable  :: topography(:, :), lon_topo(:), lat_topo(:)
  
  real(kind = fp)               :: ttime_tmp, width_tmp, velocity_interpolate, qinv_interpolate, &
  &                                az_tmp, inc_angle_tmp, inc_angle_min, az_new, inc_angle_new, az_ini, &
  &                                lon_tmp, lat_tmp, depth_tmp, lon_new, lat_new, depth_new, origintime, &
  &                                lon_grid, lat_grid, depth_grid, dist_min, dist_tmp, dvdz, epdelta, &
  &                                ot_begin, ot_end, ot_shift, rms_tw, lon_min, lat_min, depth_min, &
  &                                dinc_angle, dinc_angle_org, inc_angle_ini
  real(kind = dp)               :: topography_interpolate, dlon_topo, dlat_topo
  
  integer                       :: i, j, k, ii, jj, kk, icount, wave_index, time_count, lon_index, lat_index, z_index, &
  &                                npts_max, nlon_topo, nlat_topo
  character(len = 129)          :: dem_file, sacfile, sacfile_index, ot_begin_t, ot_end_t, rms_tw_t, ot_shift_t, &
  &                                grdfile, resultfile, resultdir
  character(len = maxlen)       :: time_count_char

  !!OpenMP variable
  !$ integer                    :: omp_thread

  icount = iargc()
  if(icount .ne. 4) then
    write(0, '(a)', advance="no") "usage: ./asl_masterevent "
    write(0, '(a)')               "(topography_grd) (station_param_file) (masterevent_param_file) (subevent_param_file)"
    error stop
  endif
  call getarg(1, topo_grd)
  call getarg(2, station_param)
  call getarg(3, masterevent_param)
  call getarg(4, subevent_param)

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
  allocate(stlon(nsta), stlat(nsta), stdp(nsta))
  do i = 1, nsta
    read(10, *) stlon(i), stlat(i), stdp(i)
  enddo
  close(10)
  !!read masterevent parameter
  allocate(obsamp_master(nsta))
  open(unit = 10, file = masterevent_param)
  read(10, *) evlon_master, evlat_master, evdp_master
  read(10, *) sourceamp_master
  read(10, *) (obsamp_master(i), i = 1, nsta)
  close(10)
  !!read subevent paramter
  open(unit = 10, file = subevent_param_file)
  read(10, *) nsubevent
  allocate(obsamp_sub(1 : nsta, 1 : nsubevent)
  do j = 1, nsubevent
    read(10, *) (obsamp_sub(i, j), i = 1, nsta)
  enddo
  close(10)

  !!set velocity/attenuation structure
  call set_velocity(z_min, dz, velocity, qinv)

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

  allocate(hypodist(1 : nsta), ray_azinc(1 : 2, 1 : nsta))
  
  !$omp parallel default(none), &
  !$omp&         shared(topography, lat_sta, lon_sta, z_sta, velocity, qinv, ttime_min, width_min, hypodist, &
  !$omp&                lon_topo, lat_topo, dlon_topo, dlat_topo), &
  !$omp&         private(i, j, k, ii, jj, kk, lon_grid, lat_grid, depth_grid, az_ini, epdelta, lon_index, lat_index, z_index, &
  !$omp&                dist_min, inc_angle_tmp, inc_angle_min, lon_tmp, lat_tmp, depth_tmp, az_tmp, val_2d, &
  !$omp&                topography_interpolate, ttime_tmp, width_tmp, xgrid, ygrid, zgrid, dist_tmp, val_1d, &
  !$omp&                velocity_interpolate, val_3d, qinv_interpolate, dvdz, lon_new, lat_new, depth_new, &
  !$omp&                az_new, inc_angle_new, omp_thread, lon_min, lat_min, depth_min, dinc_angle, dinc_angle_org, &
  !$omp&                inc_angle_ini, inc_angle_ini_min)

  !$ omp_thread = omp_get_thread_num()

  !!calculate traveltime to station        
  !$omp do
  station_loop: do jj = 1, nsta
    !!calculate azimuth and hypocentral distance
    call greatcircle_dist(evlat_master, evlon_master, stlat(jj), stlon(jj), &
    &                     distance = epdist, azimuth = az_ini, delta_out = epdelta)
    lon_index = int((stlon(jj) - lon_w) / dlon) + 1
    lat_index = int((stlat(jj) - lat_s) / dlat) + 1
    z_index   = int((stdp(jj) - z_min / dz) + 1
    !print *, lon_sta(jj), lon_w + real(lon_index - 1) * dlon
    !print *, lat_sta(jj), lat_s + real(lat_index - 1) * dlat
    hypodist(jj) = sqrt((r_earth - evdp_master) ** 2 + (r_earth - stdp(jj)) ** 2 &
    &            - 2.0_fp * (r_earth - evdp_master) * (r_earth - stdp(jj)) * cos(epdelta))
    ray_azinc(1, jj) = az_ini

#ifdef V_CONST
    !!homogeneous structure: using velocity/qinv at the grid
    ttime_min = hypodist(jj) / velocity(lon_index, lat_index, z_index)
    width_min(jj) = ttime_min(jj) * qinv(lon_index, lat_index, z_index)
    ray_azinc(2, jj) = acos(epdist / hypodist(jj)) + pi / 2.0_fp

#else
          
    !!do ray shooting
    dist_min = huge
    ttime_min = real(huge, kind = fp)
    width_min(jj) = real(huge, kind = fP)
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
          xgrid(1 : 2) = [lon_topo(lon_index), lon_topo(lon_index + 1)]
          ygrid(1 : 2) = [lat_topo(lat_index), lat_topo(lat_index + 1)]
          val_2d(1 : 2, 1 : 2) = topography(lon_index : lon_index + 1, lat_index : lat_index + 1)
          call linear_interpolation_2d(lon_tmp, lat_tmp, xgrid, ygrid, val_2d, topography_interpolate)
          if(depth_tmp .lt. topography_interpolate) then
            !print '(a, 3(f9.4, 1x))', "ray surface arrived, lon/lat = ", lon_tmp, lat_tmp, depth_tmp
            exit shooting_loop
          endif

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
          call greatcircle_dist(lat_tmp, lon_tmp, stlat(jj), stlon(jj), delta_out = epdelta)
          dist_tmp = sqrt((r_earth - depth_tmp) ** 2 + (r_earth - stdp(jj)) ** 2 &
          &        - 2.0_fp * (r_earth - depth_tmp) * (r_earth - stdp(jj)) * cos(epdelta))
          !print *, real(ii, kind = fp) * dinc_angle, lat_tmp, lon_tmp, depth_tmp
          !print *, jj, lat_sta(jj), lon_sta(jj), z_sta(jj), dist_tmp

          if(dist_tmp .lt. dist_min) then
            dist_min = dist_tmp
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
          val_3d(1 : 2, 1 : 2, 1 : 2) = qinv(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
          dvdz = (val_1d(2) - val_1d(1)) / dz
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
    !print '(a, 4(f8.4, 1x))', "grid lon, lat, depth, az_ini = ", lon_grid, lat_grid, depth_grid, az_ini * rad2deg
    !print '(a, 3(f8.4, 1x))', "station lon, lat, depth = ", lon_sta(jj), lat_sta(jj), z_sta(jj)
    !print '(a, 4(f8.4, 1x))', "rayshoot lon, lat, depth, inc_angle = ", lon_min, lat_min, depth_min, &
    !&                                                                   inc_angle_ini_min(nrayshoot) * rad2deg
    !print '(a, 3(f8.4, 1x))', "dist_min, ttime, width = ", dist_min, ttime_min(jj, i, j, k), width_min(jj, i, j, k)
    if(dist_min .gt. rayshoot_dist_thr) then
      width_min(jj, i, j, k) = hu
    endif
#endif
  enddo station_loop
  !$omp end do
  !$omp end parallel

  !!

  stop
end program AmplitudeSourceLocation_masterevent

