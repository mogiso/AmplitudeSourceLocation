program AmplitudeSourceLocation_PulseWidth
  !!Amplitude Source Location using depth-dependent 1D velocity structure, 3D heterogeneous attenuation structure
  use nrtype,               only : fp, sp, dp
  use constants,            only : rad2deg, deg2rad, pi, r_earth
  use rayshooting,          only : rayshooting3D
  use set_velocity_model,   only : set_velocity
  use linear_interpolation, only : linear_interpolation_1d
  use greatcircle,          only : greatcircle_dist

  implicit none

  real(kind = fp),    parameter :: dlon = 0.001_fp, dlat = 0.001_fp, dz = 0.1_fp
  real(kind = fp),    parameter :: lon_w = 143.980_fp, lon_e = 144.041_fp
  real(kind = fp),    parameter :: lat_s = 43.365_fp, lat_n = 43.410_fp
  real(kind = fp),    parameter :: z_min = -1.5_fp, z_max = 3.0_fp
  real(kind = fp),    parameter :: dvdlon = 0.0_fp, dvdlat = 0.0_fp
  integer,            parameter :: ninc_angle = 200
  real(kind = fp),    parameter :: freq = 7.5_fp
  real(kind = fp),    parameter :: rms_tw = 10.0_fp
  integer,            parameter :: npts_max = 3600 * 200
  integer,            parameter :: nsta = 4
  character(len = 6), parameter :: stname(1 : 4) = (/"V.MEAB", "V.MEAA", "V.PMNS", "V.NSYM"/)
  real(kind = dp),    parameter :: siteamp(1 : 4) = (/1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp/)
  real(kind = fp),    parameter :: ot_begin = 0.0_fp, ot_shift = 5.0_fp
  real(kind = fp),    parameter :: time_step = 0.01_fp

  real(kind = fp),    parameter :: dinc_angle = pi / real(ninc_angle, kind = fp)
  integer,            parameter :: nlon = int((lon_e - lon_w) / dlon) + 1
  integer,            parameter :: nlat = int((lat_n - lat_s) / dlat) + 1
  integer,            parameter :: nz   = int((z_max - z_min) / dz) + 1

  real(kind = fp)               :: velocity(1 : nlon, 1 : nlat, 1 : nz), qinv(1 : nlon, 1 : nlat, 1 : nz), &
  &                                topography(1 : nlon, 1 : nlat), width_min(1 : nsta), &
  &                                sampling(1 : nsta), begin(1 : nsta), &
  &                                val_1d(1 : 2), val_2d(1 : 2, 1 : 2), val_3d(1 : 2, 1 : 2, 1 : 2), &
  &                                xgrid(1 : 2), ygrid(1 : 2), zgrid(1 : 2), lon_sta(1 : nsta), lat_sta(1 : nsta), &
  &                                z_sta(1 : nsta)
  real(kind = dp)               :: waveform_obs(1 : npts_max, 1 : nsta), rms_amp_obs(1 : nsta), hypodist(1 : nsta), &
  &                                residual(1 : nlon, 1 : nlat, 1 : nz)
  integer                       :: npts(1 : nsta), lon_index, lat_index, z_index
  
  real(kind = fp)               :: ttime_tmp, ttime_min, width_tmp, velocity_interpolate, qinv_interpolate, &
  &                                topography_interpolate, az_tmp, inc_angle_tmp, az_new, inc_angle_new, &
  &                                lon_tmp, lat_tmp, depth_tmp, lon_new, lat_new, depth_new, origintime, &
  &                                lon_grid, lat_grid, depth_grid, dist_min, dist_tmp, dvdz, epdelta
  real(kind = sp)               :: lon_r, lat_r, topo_r
  real(kind = dp)               :: source_amp, residual_normalize
  
  integer                       :: i, j, k, ii, jj, icount, wave_index
  character(len = 129)          :: dem_file, sacfile, sacfile_index
  

  !!read topography file
  !!pay attention to the order
  open(unit = 10, file = dem_file, form = "unformatted", access = "direct", recl = 12)
  icount = 1
  do j = nlat, 1, -1
    do i = 1, nlon
      read(10, rec = icount) lon_r, lat_r, topo_r
      topography(i, j) = real(topo_r, kind = fp)
      icount = icount + 1
    enddo
  enddo
  close(10)

  !!read sac file
  do j = 1, nsta
    sacfile = sacfile_index // "." // trim(stname(i)) // ".sac"
    call read_sachdr(sacfile, delta=sampling(j), stlat = lat_sta(j), stlon = lon_sta(j), npts = npts(j), begin = begin(j))
    if(j .ne. 1) then
      do i = 1, j - 1
        if(begin(j) .ne. begin(i)) then
          write(0, '(3a)') "beginning time is different: ", trim(stname(j)), trim(stname(i))
          stop
        endif
      enddo
    endif
    call read_sacdata(sacfile, npts_max, waveform_obs(:, j))
  enddo

  time_loop: do
    z_loop: do k = 1, nz - 1
      depth_grid = z_min + dz * real(k - 1, kind = fp)
      lat_loop: do j = 1, nlat - 1
        lat_grid = lat_s + dlat * real(j - 1, kind = fp)
        lon_loop: do i = 1, nlon - 1
          lon_grid = lon_w + dlon * real(i - 1, kind = fp)
  
          !!check the grid is lower than the topo
          if(depth_grid .lt. topography(i, j)) exit
  
          !!calculate traveltime to station        
          station_loop: do jj = 1, nsta
            !!calculate azimuth and hypocentral distance
            call greatcircle_dist(lat_grid, lon_grid, lat_sta(jj), lon_sta(jj), azimuth = az_tmp, delta_out = epdelta)
            lon_index = int((lon_sta(jj) - lon_w) / dlon) + 1
            lat_index = int((lat_sta(jj) - lat_s) / dlat) + 1
            xgrid(1) = lon_e + real(lon_index - 1, kind = fp) * dlon; xgrid(2) = xgrid(1) + dlon
            ygrid(1) = lat_s + real(lat_index - 1, kind = fp) * dlat; ygrid(2) = ygrid(1) + dlat
            val_2d(1 : 2, 1 : 2) = topography(lon_index : lon_index + 1, lat_index : lat_index + 1)
            call linear_interpolation_2d(lon_sta(jj), lat_sta(jj), xgrid, ygrid, val_2d, z_sta(jj))
            hypodist(jj) = sqrt((r_earth - depth_grid) ** 2 + (r_earth - z_sta(jj)) ** 2 &
            &                   - 2.0_fp * (r_earth - depth_grid) * (r_earth - z_sta(jj)) * cos(epdelta))
            
            !!do ray shooting
            lon_tmp = lon_grid
            lat_tmp = lat_grid
            depth_tmp = depth_grid
            ttime_min = 1.0e+10_fp
            width_min(jj) = 0.0_fp
  
            dist_min = 1.0e+10_fp
            incangle_loop: do ii = 1, ninc_angle
              inc_angle_tmp = real(ii - 1) * dinc_angle
              if(inc_angle_tmp .lt. 30.0_fp * deg2rad) then
                cycle
              endif
  
              ttime_tmp = 0.0_fp
              width_tmp = 0.0_fp
              !!loop until ray arrives at surface/boundary
              shooting_loop: do
                lon_index = int((lon_tmp - lon_w) / dlon) + 1
                lat_index = int((lat_tmp - lat_s) / dlat) + 1
                z_index   = int((depth_tmp - z_min) / dz) + 1
  
                if(lon_index .lt. 1 .or. lon_index .gt. nlon - 1      &
                &  .or. lat_index .lt. 1 .or. lat_index .gt. nlat - 1 &
                &  .or. z_index .lt. 1 .or. z_index .gt. nz - 1) then
                  dist_tmp = 1.0e+11_fp
                  exit shooting_loop
                endif
  
                xgrid(1) = lon_e + real(lon_index - 1, kind = fp) * dlon; xgrid(2) = xgrid(1) + dlon
                ygrid(1) = lat_s + real(lat_index - 1, kind = fp) * dlat; ygrid(2) = ygrid(1) + dlat
                zgrid(1) = z_min + real(z_index - 1, kind = fp) * dz;     zgrid(2) = zgrid(1) + dz
  
                !!exit if the ray arrived at the surface
                val_2d(1 : 2, 1 : 2) = topography(lon_index : lon_index + 1, lat_index : lat_index + 1)
                call linear_interpolation_2d(lon_tmp, lat_tmp, xgrid, ygrid, val_2d, topography_interpolate)
                if(depth_tmp .le. topography_interpolate) then
                  !!calculate distance between station and ray arrived
                  call greatcircle_dist(lat_tmp, lon_tmp, lat_sta(jj), lon_sta(jj), distance = dist_tmp)
                  exit shooting_loop
                endif
   
                !!otherwise shooting the ray
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
  
              if(dist_tmp .lt. dist_min) then
                dist_min = dist_tmp
                ttime_min = ttime_tmp
                width_min(jj) = width_tm
            enddo incangle_loop
  
            !!caluculate site-corrected amplitude
            !!caluculate index of waveform array
            wave_index = int((origintime + ttime_min) / sampling(jj) + 0.5_fp) + 1
            rms_amp_obs(jj) = 0.0_fp
            icount = 0
            do ii = wave_index, wave_index + int(rms_tw / sampling(jj) + 0.5_fp)
              if(ii .lt. npts(jj)) then
                rms_amp_obs(jj) = rms_amp_obs(jj) + waveform_obs(ii, jj) * waveform_obs(ii, jj)
                icount = icount + 1
              endif
            enddo
            rms_amp_obs(jj) = sqrt(rms_amp_obs(jj) / real(icount, kind = fp))
          enddo station_loop
  
          !!calculate source amplitude and residual
          source_amp = 0.0_dp
          do ii = 1, nsta
            source_amp = source_amp &
            &            + rms_amp_obs(ii) / siteamp(ii) * hypodist(ii) * exp(width_min(ii) * (pi * freq)) / real(nsta, kind = fp)
          enddo
          residual(i, j, k) = 0.0_dp
          residual_normalize = 0.0_dp
          do ii = 1, nsta
            residual(i, j, k) = residual(i, j, k) &
            &                 + (rms_amp_obs(ii) / siteamp(ii) - source_amp / hypodist(ii) * exp(-pi * freq * width_min(ii))) ** 2
            residual_normalize = residual_normalize + (rms_amp_obs(ii) / siteamp(ii)) ** 2
          enddo
          residual(i, j, k) = residual(i, j, k) / residual_normalize
  
        enddo lon_loop
      enddo lat_loop
    enddo z_loop
  enddo time_loop

  stop
end program AmplitudeSourceLocation_PulseWidth
