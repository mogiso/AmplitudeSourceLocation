program AmplitudeSourceLocation_synthwave
  !!Making synthetic waveform using depth-dependent 1D velocity structure, 3D heterogeneous attenuation structure
  !!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2020
  !!License  : MIT License (https://opensource.org/licenses/MIT)

  use nrtype,               only : fp, sp, dp
  use constants,            only : rad2deg, deg2rad, pi, r_earth
  use rayshooting,          only : rayshooting3D
  use read_sacfile,         only : read_sachdr
  use set_velocity_model,   only : set_velocity
  use linear_interpolation, only : linear_interpolation_1d, linear_interpolation_2d, block_interpolation_3d
  use greatcircle,          only : greatcircle_dist
  use grdfile_io,           only : read_grdfile_2d
  use xorshift1024star
  !$ use omp_lib

  implicit none

  integer,            parameter :: wavetype = 2      !!1 for P-wave, 2 for S-wave
  !!Structure range
  real(kind = fp),    parameter :: lon_str_w = 143.5_fp, lon_str_e = 144.1_fp
  real(kind = fp),    parameter :: lat_str_s = 43.0_fp, lat_str_n = 43.5_fp
  real(kind = fp),    parameter :: z_str_min = -1.5_fp, z_str_max = 10.0_fp
  real(kind = fp),    parameter :: dlon_str = 0.001_fp, dlat_str = 0.001_fp, dz_str = 0.1_fp
  !!Ray shooting
  real(kind = fp),    parameter :: dvdlon = 0.0_fp, dvdlat = 0.0_fp         !!assume 1D structure
  integer,            parameter :: ninc_angle = 200                         !!grid search in incident angle
  integer,            parameter :: nrayshoot = 2                            !!number of grid search
  real(kind = fp),    parameter :: time_step = 0.01_fp
  real(kind = fp),    parameter :: rayshoot_dist_thr = 0.05_fp
  !!assumed hypocenter
  integer,            parameter :: nhypo = 1
  real(kind = fp),    parameter :: origintime(1 : nhypo) = [2.0_fp]
  real(kind = fp),    parameter :: lon_hypo(1 : nhypo)   = [144.0183_fp]
  real(kind = fp),    parameter :: lat_hypo(1 : nhypo)   = [43.3902_fp]
  real(kind = fp),    parameter :: depth_hypo(1 : nhypo) = [0.1_fp]
  real(kind = dp),    parameter :: amp_hypo(1 : nhypo)   = [1.3_dp]
  real(kind = dp),    parameter :: characteristicfreq    = 0.5_dp
  real(kind = dp),    parameter :: asl_freq              = 7.5_dp

  !!assumed data length, sampling frequency (s)
  real(kind = fp),    parameter :: waveform_length = 100.0_fp
  real(kind = fp),    parameter :: wavelet_length = 15.0_fp
  real(kind = dp),    parameter :: sampling = 0.01_dp
  integer,            parameter :: npts_waveform = int(waveform_length / real(sampling, kind = fp) + 0.5_fp)
  integer,            parameter :: npts_wavelet = int(wavelet_length / real(sampling, kind = fp) + 0.5_fp)
  real(kind = dp),    parameter :: sigma_noise = 0.01_dp

  real(kind = fp),    parameter :: alt_to_depth = -1.0e-3_fp
  real(kind = dp),    parameter :: huge = 1.0e+5_dp

  !real(kind = fp),    parameter :: dinc_angle1 = pi / real(ninc_angle, kind = fp)
  !real(kind = fp),    parameter :: dinc_angle2 = 2.0_fp * dinc_angle1 / real(ninc_angle - 1, kind = fp)
  integer,            parameter :: nlon_str = int((lon_str_e - lon_str_w) / dlon_str) + 2
  integer,            parameter :: nlat_str = int((lat_str_n - lat_str_s) / dlat_str) + 2
  integer,            parameter :: nz_str   = int((z_str_max - z_str_min) / dz_str) + 2

  real(kind = fp)               :: velocity(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                qinv(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                val_1d(1 : 2), val_2d(1 : 2, 1 : 2), val_3d(1 : 2, 1 : 2, 1 : 2), &
  &                                xgrid(1 : 2), ygrid(1 : 2), zgrid(1 : 2), inc_angle_ini_min(0 : nrayshoot)
  real(kind = dp), allocatable  :: topography(:, :), lon_topo(:), lat_topo(:), waveform(:), wavelet(:)
  real(kind = fp), allocatable  :: lon_sta(:), lat_sta(:), z_sta(:), ttime_min(:, :), width_min(:, :), hypodist(:, :)
  
  real(kind = fp)               :: ttime_tmp, width_tmp, velocity_interpolate, qinv_interpolate, &
  &                                az_tmp, inc_angle_tmp, az_new, inc_angle_new, az_ini, &
  &                                lon_tmp, lat_tmp, depth_tmp, lon_new, lat_new, depth_new, &
  &                                dist_min, dist_tmp, dvdz, epdelta, lon_min, lat_min, depth_min, &
  &                                dinc_angle, dinc_angle_org, inc_angle_ini
  real(kind = dp)               :: topography_interpolate, dlon_topo, dlat_topo, random_number1, random_number2
  
  integer                       :: nsta, i, j, ii, jj, icount, wave_index, lon_index, lat_index, z_index, nlon_topo, nlat_topo
  character(len = 129)          :: dem_file
  character(len = 129), allocatable :: sacfile(:)
  character(len = 4),   allocatable :: sachdr(:, :)

  !!random number
  type(xorshift1024star_state)  :: random_state
  
  !!OpenMP variable
  !$ integer                    :: omp_thread

  icount = iargc()
  if(icount .le. 1) then
    write(0, '(a)') "usage: ./a.out dem_grdfile_name sacfile1 sacfile2 ..."
    error stop
  endif
  nsta = icount - 1
  allocate(sacfile(1 : nsta), sachdr(1 : 158, 1 : nsta), lon_sta(1 : nsta), lat_sta(1 : nsta), z_sta(1 : nsta), &
  &        ttime_min(1 : nhypo, 1 : nsta), width_min(1 : nhypo, 1 : nsta), hypodist(1 : nhypo, 1 : nsta))
  
  call getarg(1, dem_file)
  do i = 1, nsta
    call getarg(i + 1, sacfile(i))
  enddo

  write(0, '(a, 3(1x, f8.3))') "lon_str_w, lat_str_s, z_str_min =", lon_str_w, lat_str_s, z_str_min
  write(0, '(a, 3(1x, i0))') "nlon_str, nlat_str, nz_str =", nlon_str, nlat_str, nz_str

  !!read topography file (netcdf grd format)
  call read_grdfile_2d(dem_file, lon_topo, lat_topo, topography)
  nlon_topo = ubound(lon_topo, 1)
  nlat_topo = ubound(lat_topo, 1)
  dlon_topo = lon_topo(2) - lon_topo(1)
  dlat_topo = lat_topo(2) - lat_topo(1)
  topography(1 : nlon_topo, 1 : nlat_topo) = topography(1 : nlon_topo, 1 : nlat_topo) * alt_to_depth

  !!set velocity/attenuation structure
  call set_velocity(z_str_min, dz_str, velocity, qinv)

  !!read sac file
  do i = 1, nsta
    call read_sachdr(sacfile(i), header = sachdr(:, i), stlat = lat_sta(i), stlon = lon_sta(i), stdp = z_sta(i))
#ifdef STDP_COR
    z_sta(i) = z_sta(i) * (alt_to_depth)
#endif
  enddo

  !!make traveltime/pulse width table for each grid point
  write(0, '(a)') "calculate traveltime / pulse width ..."
  !$omp parallel default(none), &
  !$omp&         shared(topography, nsta, lat_sta, lon_sta, z_sta, velocity, qinv, &
  !$omp&                ttime_min, width_min, hypodist, lon_topo, lat_topo, dlon_topo, dlat_topo, nlon_topo, nlat_topo), &
  !$omp&         private(i, ii, jj, az_ini, epdelta, lon_index, lat_index, z_index, &
  !$omp&                dist_min, inc_angle_tmp, lon_tmp, lat_tmp, depth_tmp, az_tmp, val_2d, &
  !$omp&                topography_interpolate, ttime_tmp, width_tmp, xgrid, ygrid, zgrid, dist_tmp, val_1d, &
  !$omp&                velocity_interpolate, val_3d, qinv_interpolate, dvdz, lon_new, lat_new, depth_new, &
  !$omp&                az_new, inc_angle_new, omp_thread, lon_min, lat_min, depth_min, dinc_angle, dinc_angle_org, &
  !$omp&                inc_angle_ini, inc_angle_ini_min)

  !$ omp_thread = omp_get_thread_num()

  !$omp do schedule(guided)
  !!calculate traveltime to station        
  station_loop: do j = 1, nsta

    hypo_loop: do i = 1, nhypo
      ttime_min(i, j) = huge
      width_min(i, j) = huge

      !!check the grid is lower than the topo
      lon_index = int((lon_hypo(i) - lon_topo(1)) / dlon_topo) + 1
      lat_index = int((lat_hypo(i) - lat_topo(1)) / dlat_topo) + 1
      xgrid(1 : 2) = [lon_topo(lon_index), lon_topo(lon_index + 1)]
      ygrid(1 : 2) = [lat_topo(lat_index), lat_topo(lat_index + 1)]
      val_2d(1 : 2, 1 : 2) = topography(lon_index : lon_index + 1, lat_index : lat_index + 1)
      call linear_interpolation_2d(lon_hypo(i), lat_hypo(i), xgrid, ygrid, val_2d, topography_interpolate)
      if(depth_hypo(i) .lt. topography_interpolate) then
        write(0, '(a, i0)') "depth_hypo is higher than altitude there, hypo_index = ", i
        error stop
      endif

      !!calculate azimuth and hypocentral distance
      call greatcircle_dist(lat_hypo(i), lon_hypo(i), lat_sta(j), lon_sta(j), azimuth = az_ini, delta_out = epdelta)
      !lon_index = int((lon_sta(j) - lon_w) / dlon) + 1
      !lat_index = int((lat_sta(j) - lat_s) / dlat) + 1
      !print *, lon_sta(jj), lon_w + real(lon_index - 1) * dlon
      !print *, lat_sta(jj), lat_s + real(lat_index - 1) * dlat
      hypodist(i, j) = sqrt((r_earth - depth_hypo(i)) ** 2 + (r_earth - z_sta(j)) ** 2 &
      &              - 2.0_fp * (r_earth - depth_hypo(i)) * (r_earth - z_sta(j)) * cos(epdelta))

#ifdef V_CONST
      !!homogeneous structure: using velocity/qinv at the grid
      lon_index = int((lon_hypo(i) - lon_str_w) / dlon_str) + 1
      lat_index = int((lat_hypo(i) - lat_str_s) / dlat_str) + 1
      z_index = int((depth_hypo(i) - depth_min) / dz_str) + 1
      ttime_min(i, j) = hypodist(i, j) / velocity(lon_index, lat_index, z_index)
      width_min(i, j) = ttime_min(i, j) * qinv(lon_index, lat_index, z_index)

#else
          
      !!do ray shooting
      dist_min = huge
      incangle_loop2: do jj = 1, nrayshoot
        if(jj .eq. 1) then
          dinc_angle_org = pi / 2.0_fp
        else
          dinc_angle_org = dinc_angle
        endif
        dinc_angle = 2.0_fp * dinc_angle_org / real(ninc_angle, kind = fp)
        inc_angle_ini_min(0) = dinc_angle_org
        !print *, "dinc_angle = ", dinc_angle * rad2deg, inc_angle_ini_min(kk - 1) * rad2deg

        incangle_loop: do ii = 1, ninc_angle
          inc_angle_ini = (inc_angle_ini_min(jj - 1) - dinc_angle_org) + real(ii, kind = fp) * dinc_angle
          !print '(2(i0, 1x), a, e15.7)', ii, kk, "inc_angle_ini = ", inc_angle_ini * rad2deg
          
          lon_tmp = lon_hypo(i)
          lat_tmp = lat_hypo(i)
          depth_tmp = depth_hypo(i)
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
            if(lon_index .lt. 1 .or. lon_index .gt. nlon_topo - 1 &
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
            call greatcircle_dist(lat_tmp, lon_tmp, lat_sta(j), lon_sta(j), delta_out = epdelta)
            dist_tmp = sqrt((r_earth - depth_tmp) ** 2 + (r_earth - z_sta(j)) ** 2 &
            &        - 2.0_fp * (r_earth - depth_tmp) * (r_earth - z_sta(j)) * cos(epdelta))
            !print *, real(ii, kind = fp) * dinc_angle, lat_tmp, lon_tmp, depth_tmp
            !print *, jj, lat_sta(jj), lon_sta(jj), z_sta(jj), dist_tmp

            if(dist_tmp .lt. dist_min) then
              dist_min = dist_tmp
              ttime_min(i, j) = ttime_tmp
              width_min(i, j) = width_tmp
              lon_min = lon_tmp
              lat_min = lat_tmp
              depth_min = depth_tmp
              inc_angle_ini_min(jj) = inc_angle_ini
              !print '(a, 4(f8.4, 1x))', "rayshoot_tmp lon, lat, depth, dist_min = ", lon_min, lat_min, depth_min, &
              !&                                                                      dist_min
            endif
 
            !!shooting the ray
            val_1d(1 : 2) = velocity(lon_index, lat_index, z_index : z_index + 1, wavetype)
            call linear_interpolation_1d(depth_tmp, zgrid, val_1d, velocity_interpolate)
            val_3d(1 : 2, 1 : 2, 1 : 2) &
            &  = qinv(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1, wavetype)
            call block_interpolation_3d(lon_tmp, lat_tmp, depth_tmp, xgrid, ygrid, zgrid, val_3d, qinv_interpolate)

            dvdz = (val_1d(2) - val_1d(1)) / dz_str
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
      print '(a, 4(f8.4, 1x))', "hypo lon, lat, depth, az_ini = ", lon_hypo(i), lat_hypo(i), depth_hypo(i), az_ini * rad2deg
      print '(a, 3(f8.4, 1x))', "station lon, lat, depth = ", lon_sta(j), lat_sta(j), z_sta(j)
      print '(a, 4(f8.4, 1x))', "rayshoot lon, lat, depth, inc_angle = ", lon_min, lat_min, depth_min, &
      &                                                                   inc_angle_ini_min(nrayshoot) * rad2deg
      print '(a, 3(f8.4, 1x))', "dist_min, ttime, width = ", dist_min, ttime_min(i, j), width_min(i, j)
      if(dist_min .gt. rayshoot_dist_thr) then
        ttime_min(i, j) = huge
        width_min(i, j) = 0.0_fp
        write(0, '(a, i0)') "rayshooting failed, hypo_index = ", i
        error stop
      endif
#endif
    enddo hypo_loop
  enddo station_loop
  !$omp end do
  !$omp end parallel

  !!make waveform
  allocate(waveform(1 : npts_waveform))
  allocate(wavelet(1 : npts_wavelet))
  call state_init_self(random_state)
  do j = 1, nsta
    write(0, '(2a)') "waveform file = ", trim(sacfile(j))
    waveform(1 : npts_waveform) = 0.0_dp
    do i = 1, nhypo
      call ricker_wavelet(npts_wavelet, sampling, characteristicfreq, amp_hypo(i), 0.0_dp, wavelet)
      wave_index = int((origintime(i) + ttime_min(i, j)) / real(sampling, kind = fp) + 0.5_fp) + 1
      do ii = 1, npts_wavelet
        if(wave_index + ii - 1 .le. npts_waveform) then
          waveform(wave_index + ii - 1) = wavelet(ii) / hypodist(i, j) * exp(-pi * asl_freq * width_min(i, j))
        endif
      enddo
      write(0, '(a, i0, a, 2(e15.7, 1x))') "hypo index = ", i, " theo. amp. and travel time = ", &
      &                             amp_hypo(i) / hypodist(i, j) * exp(-pi * asl_freq * width_min(i, j)), &
      &                             ttime_min(i, j)
    enddo


    do i = 1, npts_waveform
      random_number1 = draw_uniform(random_state)
      random_number2 = draw_uniform(random_state)
      waveform(i) = waveform(i) &
      &           + sqrt(-log(random_number1 * random_number1)) * sin(2.0_dp * pi * random_number2) * sigma_noise
    enddo

    write(0, '(2a)') "writing waveform to ", trim(sacfile(j))
    open(unit = 10, file = trim(sacfile(j)), form = "unformatted", access = "direct", recl = 4)
    do i = 1, 158
      write(10, rec = i) sachdr(i, j)
    enddo
    write(10, rec = 1) real(sampling, kind = sp)
    !!reset origin time, travel time
    write(10, rec = 6) 0.0_sp
    write(10, rec = 7) real(sampling, kind = sp) * real(npts_waveform - 1, kind = sp)
    write(10, rec = 8) -12345.0_sp
    write(10, rec = 9) real(ttime_min(1, j) + origintime(1), kind = sp)
    write(10, rec = 11) -12345.0_sp
    write(10, rec = 60) 0.0_sp
    write(10, rec = 61) real(sampling, kind = sp) * real(npts_waveform - 1, kind = sp)
    write(10, rec = 62) real(minval(waveform), kind = sp)
    write(10, rec = 63) real(maxval(waveform), kind = sp)
    write(10, rec = 80) npts_waveform
    do i = 1, npts_waveform
      write(10, rec = 158 + i) real(waveform(i), kind = sp)
    enddo
    close(10)
    

  enddo


  stop
end program AmplitudeSourceLocation_synthwave

