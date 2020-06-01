program AmplitudeSourceLocation_PulseWidth
  !!Amplitude Source Location using depth-dependent 1D velocity structure, 3D heterogeneous attenuation structure
  !$use omp
  use nrtype,               only : fp, sp, dp
  use constants,            only : rad2deg, deg2rad, pi, r_earth
  use rayshooting,          only : rayshooting3D
  use read_sacfile,         only : read_sachdr, read_sacdata
  use set_velocity_model,   only : set_velocity
  use linear_interpolation, only : linear_interpolation_1d, linear_interpolation_2d, block_interpolation_3d
  use greatcircle,          only : greatcircle_dist
  use itoa,                 only : int_to_char
  use GMT

  implicit none

  real(kind = fp),    parameter :: dlon = 0.001_fp, dlat = 0.001_fp, dz = 0.1_fp
  real(kind = fp),    parameter :: lon_w = 143.90_fp, lon_e = 144.05_fp
  real(kind = fp),    parameter :: lat_s = 43.35_fp, lat_n = 43.410_fp
  real(kind = fp),    parameter :: z_min = -1.5_fp, z_max = 3.0_fp
  real(kind = fp),    parameter :: dvdlon = 0.0_fp, dvdlat = 0.0_fp
  integer,            parameter :: ninc_angle = 200
  real(kind = fp),    parameter :: rms_tw = 10.0_fp
  integer,            parameter :: npts_max = 7200 * 200
  integer,            parameter :: nsta = 4
  character(len = 6), parameter :: stname(1 : nsta) = (/"V.MEAB", "V.MEAA", "V.PMNS", "V.NSYM"/)
  real(kind = dp),    parameter :: siteamp(1 : nsta) = (/1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp/)
  real(kind = fp),    parameter :: time_step = 0.01_fp
  real(kind = fp),    parameter :: freq = 7.5_fp
  real(kind = dp),    parameter :: fl = 5.0_dp, fh = 10.0_dp, fs = 12.0_dp
  real(kind = dp),    parameter :: ap = 0.5_dp, as = 5.0_dp
  real(kind = dp),    parameter :: order = 1.0e-3_dp
  integer,            parameter :: maxlen = 5

  real(kind = fp),    parameter :: dinc_angle = pi / real(ninc_angle, kind = fp)
  integer,            parameter :: nlon = int((lon_e - lon_w) / dlon) + 1
  integer,            parameter :: nlat = int((lat_n - lat_s) / dlat) + 1
  integer,            parameter :: nz   = int((z_max - z_min) / dz) + 1

  real(kind = fp)               :: velocity(1 : nlon, 1 : nlat, 1 : nz), qinv(1 : nlon, 1 : nlat, 1 : nz), &
  &                                topography(1 : nlon, 1 : nlat), width_min(1 : nsta), sampling(1 : nsta), begin(1 : nsta), &
  &                                val_1d(1 : 2), val_3d(1 : 2, 1 : 2, 1 : 2), xgrid(1 : 2), ygrid(1 : 2), zgrid(1 : 2), &
  &                                lon_sta(1 : nsta), lat_sta(1 : nsta), z_sta(1 : nsta)
  real(kind = dp)               :: waveform_obs(1 : npts_max, 1 : nsta), rms_amp_obs(1 : nsta), hypodist(1 : nsta), &
  &                                residual(1 : nlon, 1 : nlat, 1 : nz), xrange(1 : 2), yrange(1 : 2), spacing(1 : 2), &
  &                                source_amp(1 : nlon, 1 : nlat, 1 : nz)
  integer                       :: npts(1 : nsta), residual_minloc(3)
  real(kind = sp), allocatable  :: residual_grd(:, :)
  
  real(kind = fp)               :: ttime_tmp, ttime_min, width_tmp, velocity_interpolate, qinv_interpolate, &
  &                                az_tmp, inc_angle_tmp, az_new, inc_angle_new, az_ini, &
  &                                lon_tmp, lat_tmp, depth_tmp, lon_new, lat_new, depth_new, origintime, &
  &                                lon_grid, lat_grid, depth_grid, dist_min, dist_tmp, dvdz, epdelta, &
  &                                ot_begin, ot_end, ot_shift
  real(kind = sp)               :: lon_r, lat_r, topo_r
  real(kind = dp)               :: residual_normalize, amp_avg
  
  integer                       :: i, j, k, ii, jj, icount, wave_index, time_count, lon_index, lat_index, z_index, grd_status
  character(len = 129)          :: dem_file, sacfile, sacfile_index, ot_begin_t, ot_end_t, ot_shift_t, grdfile
  character(len = maxlen)       :: time_count_char

  real(kind = dp), allocatable  :: h(:), waveform_tmp(:)
  real(kind = dp)               :: gn, c
  integer                       :: m, n
  
  call getarg(1, sacfile_index)
  call getarg(2, dem_file)
  call getarg(3, ot_begin_t); read(ot_begin_t, *) ot_begin
  call getarg(4, ot_end_t)  ; read(ot_end_t, *) ot_end
  call getarg(5, ot_shift_t); read(ot_shift_t, *) ot_shift

  write(0, '(a, 3(f8.3, 1x))') "lon_w, lat_s, z_min = ", lon_w, lat_s, z_min
  write(0, '(a, 3(i0))') "nlon, nlat, nz = ", nlon, nlat, nz

  !!read topography file
  !!pay attention to the order
  open(unit = 10, file = dem_file, form = "unformatted", access = "direct", recl = 12)
  icount = 1
  do j = nlat, 1, -1
    do i = 1, nlon
      read(10, rec = icount) lon_r, lat_r, topo_r
      topography(i, j) = real(topo_r, kind = fp) * (-0.001_fp)
      icount = icount + 1
    enddo
  enddo
  close(10)

  !!set velocity/attenuation structure
  call set_velocity(z_min, dz, velocity, qinv)

  !!read sac file
  do j = 1, nsta
    sacfile = trim(sacfile_index) // trim(stname(j)) // "__U__.sac"
    call read_sachdr(sacfile, delta=sampling(j), stlat = lat_sta(j), stlon = lon_sta(j), stdp = z_sta(j), &
    &                npts = npts(j), begin = begin(j))
    z_sta(j) = z_sta(j) * 0.001_fp
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
      waveform_obs(i, j) = (waveform_obs(i, j) - amp_avg) * order
    enddo
  enddo

  !!bandpass filter
  do j = 1, nsta
    !!design
    call calc_bpf_order(fl, fh, fs, ap, as, sampling(j), m, n, c)
    allocate(h(4 * m), waveform_tmp(npts(j)))
    waveform_tmp(1 : npts(j)) = waveform_obs(1 : npts(j), j)
    call calc_bpf_coef(fl, fh, sampling(j), m, n, h, c, gn)
    call tandem1(waveform_tmp, waveform_tmp, npts(j), h, m, 1)
    waveform_obs(1 : npts(j), j) = waveform_tmp(1 : npts(j)) * gn
    deallocate(h, waveform_tmp)
  enddo

 
  time_count = 0
  residual(1 : nlon, 1 : nlat, 1 : nz) = 1.0e+10_dp
  time_loop: do
    origintime = ot_begin + ot_shift * real(time_count, kind = fp)
    if(origintime .gt. ot_end) exit time_loop
    write(0, '(a, f6.1)') "origin time = ", origintime
    z_loop: do k = 1, nz - 1
      depth_grid = z_min + dz * real(k - 1, kind = fp)
      lat_loop: do j = 1, nlat - 1
        lat_grid = lat_s + dlat * real(j - 1, kind = fp)
        lon_loop: do i = 1, nlon - 1
          lon_grid = lon_w + dlon * real(i - 1, kind = fp)


          !!check the grid is lower than the topo
          if(depth_grid .lt. topography(i, j)) cycle
          !write(0, '(a, (f7.3, 1x, f6.3, 1x, f6.2), a)') "(lon, lat, z) = (", lon_grid, lat_grid, depth_grid, ")"
  
          !!calculate traveltime to station        
          station_loop: do jj = 1, nsta
            if(origintime .gt. npts(jj) * sampling(jj)) exit time_loop
            !!calculate azimuth and hypocentral distance
            call greatcircle_dist(lat_grid, lon_grid, lat_sta(jj), lon_sta(jj), azimuth = az_ini, delta_out = epdelta)
            lon_index = int((lon_sta(jj) - lon_w) / dlon) + 1
            lat_index = int((lat_sta(jj) - lat_s) / dlat) + 1
            !print *, lon_sta(jj), lon_w + real(lon_index - 1) * dlon
            !print *, lat_sta(jj), lat_s + real(lat_index - 1) * dlat
            hypodist(jj) = sqrt((r_earth - depth_grid) ** 2 + (r_earth - z_sta(jj)) ** 2 &
            &                   - 2.0_fp * (r_earth - depth_grid) * (r_earth - z_sta(jj)) * cos(epdelta))
            
            !!do ray shooting
            dist_min = 1.0e+10_fp
            ttime_min = 1.0e+10_fp
            width_min(jj) = 0.0_fp
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
                !print *, real(ii, kind = fp) * dinc_angle, lat_tmp, lon_tmp, depth_tmp, topography_interpolate
                !print *, jj, lat_sta(jj), lon_sta(jj), z_sta(jj), dist_tmp

                if(dist_tmp .lt. dist_min) then
                  dist_min = dist_tmp
                  ttime_min = ttime_tmp
                  width_min(jj) = width_tmp
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
            !print *, "station lon, lat, dist = ", lon_sta(jj), lat_sta(jj), dist_min, hypodist(jj), ttime_min
  
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
          source_amp(i, j, k) = 0.0_dp
          do ii = 1, nsta
            source_amp(i, j, k) = source_amp(i, j, k) &
            &            + rms_amp_obs(ii) / siteamp(ii) * hypodist(ii) * exp(width_min(ii) * (pi * freq)) / real(nsta, kind = fp)
          enddo
          residual(i, j, k) = 0.0_dp
          residual_normalize = 0.0_dp
          do ii = 1, nsta
            residual(i, j, k) = residual(i, j, k) &
            &                 + (rms_amp_obs(ii) / siteamp(ii) &
            &                    - source_amp(i, j, k) / hypodist(ii) * exp(-pi * freq * width_min(ii))) ** 2
            residual_normalize = residual_normalize + (rms_amp_obs(ii) / siteamp(ii)) ** 2
          enddo
          residual(i, j, k) = residual(i, j, k) / residual_normalize
  
        enddo lon_loop
      enddo lat_loop
    enddo z_loop

    !!search minimum residual
    residual_minloc = minloc(residual)
    lon_grid = lon_w + real(residual_minloc(1) - 1, kind = fp) * dlon
    lat_grid = lat_s + real(residual_minloc(2) - 1, kind = fp) * dlat
    depth_grid = z_min + real(residual_minloc(3) - 1, kind = fp) * dz
    write(0, '(a, f6.1, a, f7.3, 1x, f6.3, 1x, f6.2, 2(a, e11.7))') &
    &                 "OT = ", origintime, " residual_minimum (lon, lat, dep) = (", lon_grid, lat_grid, depth_grid, ")", &
    &                 " source_amp = ", source_amp(residual_minloc(1), residual_minloc(2), residual_minloc(3)), &
    &                 " residual = ", residual(residual_minloc(1), residual_minloc(2), residual_minloc(3))

    !!output grd files
    call int_to_char(time_count, maxlen, time_count_char)
    !!output horizontal slice
    grdfile = "min_err_lon-lat_" // trim(time_count_char) // ".grd"
    allocate(residual_grd(nlon, nlat))
    xrange(1) = lon_w
    xrange(2) = lon_w + real(nlon - 1, kind = dp)
    yrange(1) = lat_s
    yrange(2) = lat_s + real(nlat - 1, kind = dp)
    spacing(1 : 2) = (/dlon, dlat/)
    residual_grd(1 : nlon, 1 : nlat) = real(residual(1 : nlon, 1 : nlat, residual_minloc(3)), kind = sp)
    grd_status = grd_create(grdfile, residual_grd, xrange, yrange, spacing, jscan = 1, overwrite = .true.)
    deallocate(residual_grd)

    !!output depth slice -- lon-dep
    grdfile = "min_err_lon-dep_" // trim(time_count_char) // ".grd"
    allocate(residual_grd(nlon, nz))
    xrange(1) = lon_w
    xrange(2) = lon_w + real(nlon - 1, kind = dp)
    yrange(1) = z_min
    yrange(2) = z_min + real(nz - 1, kind = dp)
    spacing(1 : 2) = (/dlon, dz/)
    residual_grd(1 : nlon, 1 : nz) = real(residual(1 : nlon, residual_minloc(2), 1 : nz), kind = sp)
    grd_status = grd_create(grdfile, residual_grd, xrange, yrange, spacing, jscan = 1, overwrite = .true.)
    deallocate(residual_grd)

    !!output depth slice -- dep-lat
    grdfile = "min_err_dep-lat_" // trim(time_count_char) // ".grd"
    allocate(residual_grd(nz, nlat))
    xrange(1) = z_min
    xrange(2) = z_min + real(nz - 1, kind = dp)
    yrange(1) = lat_s
    yrange(2) = lat_s + real(nlat - 1, kind = dp)
    spacing(1 : 2) = (/dz, dlat/)
    do j = 1, nz
      do i = 1, nlat
        residual_grd(j, i) = real(residual(residual_minloc(1), i, j), kind = sp)
      enddo
    enddo
    grd_status = grd_create(grdfile, residual_grd, xrange, yrange, spacing, jscan = 1, overwrite = .true.)
    deallocate(residual_grd)

    time_count = time_count + 1
  enddo time_loop

  stop
end program AmplitudeSourceLocation_PulseWidth
