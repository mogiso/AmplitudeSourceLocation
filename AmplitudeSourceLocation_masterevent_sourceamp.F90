program AmplitudeSourceLocation_masterevent_sourceamp
  !!Estimate relative source amplitudes among several master events
  !!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2024
  !!License  : MIT License (https://opensource.org/licenses/MIT)

  use nrtype,               only : fp, dp
  use constants,            only : rad2deg, deg2rad, pi, r_earth
  use set_velocity_model,   only : set_velocity
  use linear_interpolation, only : linear_interpolation_1d, linear_interpolation_2d, linear_interpolation_3d, &
  &                                block_interpolation_3d
  use greatcircle,          only : greatcircle_dist
  use grdfile_io,           only : read_grdfile_2d
  use raybending,           only : pseudobending3D

#if defined (MKL)
  use lapack95
#else
  use f95_lapack
#endif

  implicit none

  type location
    real(kind = fp) :: lon, lat, xeast, ynorth, depth, vel, qinv, sigma_x, sigma_y, sigma_depth
  end type location
  type eventsourceamp
    real(kind = fp) :: sourceamp, sigma_sourceamp
    integer         :: evindex
  end type eventsourceamp

  integer,            parameter :: wavetype = 2          !!1 for P-wave, 2 for S-wave
  integer,            parameter :: nsta_use_minimum = 5
  real(kind = fp),    parameter :: ratio_dist_max = 10.0_fp
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
  real(kind = fp),    parameter :: lat_cor_s = 32.4_fp,  lat_cor_n = 33.2_fp
  real(kind = fp),    parameter :: z_cor_min = 5.0_fp,   z_cor_max = 10.0_fp
  real(kind = fp),    parameter :: dlon_cor = 0.05_fp, dlat_cor = 0.05_fp, dz_cor = 1.0_fp
  integer,            parameter :: nlon_cor = int((lon_cor_e - lon_cor_w) / dlon_cor) + 1
  integer,            parameter :: nlat_cor = int((lat_cor_n - lat_cor_s) / dlat_cor) + 1
  integer,            parameter :: nz_cor   = int((z_cor_max - z_cor_min) / dz_cor)   + 1

  real(kind = fp),    allocatable :: stlon(:), stlat(:), stdp(:), ttime(:, :), ttime_cor(:, :)
  character(len = 6), allocatable :: stname(:)
  logical,            allocatable :: use_flag(:)

  real(kind = fp),    parameter :: alt_to_depth = -1.0e-3_fp
  real(kind = dp),    parameter :: huge = 1.0e+6_dp

  type(location),       allocatable :: evloc_master(:)
  type(eventsourceamp), allocatable :: evsourceamp(:)

  real(kind = fp)               :: velocity(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                qinv(1 : nlon_str, 1 : nlat_str, 1 : nz_str, 1 : 2), &
  &                                val_2d(1 : 2, 1 : 2), val_3d(1 : 2, 1 : 2, 1 : 2), &
  &                                xgrid(1 : 2), ygrid(1 : 2), zgrid(1 : 2), normal_vector(1 : 3) 
  real(kind = dp),  allocatable :: topography(:, :), lon_topo(:), lat_topo(:)
  real(kind = fp),  allocatable :: obsamp_master(:, :), obsamp_noise(:), &
  &                                hypodist(:, :), ray_azinc(:, :, :), obsvector(:), obsvector_org(:), &
  &                                inversion_matrix(:, :), inversion_matrix_org(:, :), error_matrix(:, :)
  integer,          allocatable :: ipiv(:)
  logical,          allocatable :: eventpair(:, :)
  real(kind = fp)               :: epdist, epdelta, mean_lon, mean_lat, mean_depth, dist_tmp, freq, vel_mean, qinv_mean, &
  &                                ttime_tmp, matrix_const, siteamp_tmp, mean_residual, sigma_residual, attenuationcoef, &
  &                                sigma_attenuationcoef, mean_vel
  real(kind = dp)               :: topography_interpolate, dlon_topo, dlat_topo
  integer                       :: nlon_topo, nlat_topo, nsta, lon_index, lat_index, z_index, i, j, k, ii, jj, kk, &
  &                                icount, nsta_use, ios, ref_evindex, nev_master, neventpair, neventpair_max
  character(len = 129)          :: topo_grd, station_param, masterevent_param, resultfile, outfile
  character(len = 20)           :: freq_t
  character(len = 3)            :: ref_evindex_t

  !!Psuedobending parameters
  integer,                parameter :: ndiv_raypath = 10
  integer,                parameter :: nraypath_ini = 4
  real(kind = fp)                   :: raypath_lon((nraypath_ini - 1) * 2 ** ndiv_raypath + 1), &
  &                                    raypath_lat((nraypath_ini - 1) * 2 ** ndiv_raypath + 1), &
  &                                    raypath_dep((nraypath_ini - 1) * 2 ** ndiv_raypath + 1)
  integer                           :: nraypath

  !!3-D Delaunay triangulation variables
  integer,           parameter :: bf_max = 10000, fc_max = 10000
  real(kind = fp), allocatable :: vcl(:, :)
  integer,         allocatable :: vm(:), ht(:)
  integer                      :: bf(1 : 3, 1 : bf_max), fc(1 : 7, 1 : fc_max)
  integer                      :: sizht, bf_num, nfc, nface, ntetra, ierr
  

  icount = iargc()

  if(icount .ne. 6) then
    write(0, '(a)', advance="no") "usage: ./asl_masterevent_sourceamp "
    write(0, '(a)', advance="no") "(topography_grd) (station_param_file) (masterevent_param_file) "
    write(0, '(a)')               "(frequency) (reference event index) (result_file)"
    error stop
  endif
  call getarg(1, topo_grd)
  call getarg(2, station_param)
  call getarg(3, masterevent_param)
  call getarg(4, freq_t); read(freq_t, *) freq
  call getarg(5, ref_evindex_t); read(ref_evindex_t, *) ref_evindex
  call getarg(6, resultfile)

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
  &        use_flag(1 : nsta), obsamp_noise(1 : nsta))
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
  nev_master = 0
  read(10, *)
  do
    read(10, *, iostat = ios)
    if(ios .ne. 0) exit
    nev_master = nev_master + 1
  enddo
  nev_master = nev_master / 2
  print *, "nev_master = ", nev_master
  allocate(evloc_master(1 : nev_master), obsamp_master(1 : nsta, 1 : nev_master))
  rewind(10)
  read(10, *)
  mean_lon = 0.0_fp
  mean_lat = 0.0_fp
  mean_depth = 0.0_fp
  do j = 1, nev_master
    read(10, *) evloc_master(j)%lon, evloc_master(j)%lat, evloc_master(j)%depth
    read(10, *) (obsamp_master(i, j), i = 1, nsta)
    mean_lon = mean_lon + evloc_master(j)%lon
    mean_lat = mean_lat + evloc_master(j)%lat
    mean_depth = mean_depth + evloc_master(j)%depth
  enddo
  close(10)
  mean_lon = mean_lon / real(nev_master, kind = fp)
  mean_lat = mean_lat / real(nev_master, kind = fp)
  mean_depth = mean_depth / real(nev_master, kind = fp)
  write(0, '(a, 3(f9.4, 1x))') "mean lon, lat, depth = ", mean_lon, mean_lat, mean_depth
  do i = 1, nev_master
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


  !!calculate ray length, pulse width, unit vector of ray incident
  write(0, '(a)') "calculate ray length, pulse width, and ray incident vector for master events"

  !!check whether the depth of master event location is lower than the topo
  do i = 1, nev_master
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

  !!calculate traveltime from each master event to each station        
  allocate(ttime(1 : nsta, 1 : nev_master))
  allocate(hypodist(1 : nsta, 1 : nev_master), ray_azinc(1 : 2, 1 : nsta, 1 : nev_master))
  masterevent_loop: do kk = 1, nev_master
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
#endif /* -DV_CONST or not */
    enddo station_loop
  enddo masterevent_loop

  !!3-D Delaunay triangulation
  sizht = 2 * nev_master
  allocate(vcl(1 : 3, 1 : nev_master), vm(1 : nev_master), ht(0 : sizht - 1))
  do i = 1, nev_master
    vcl(1, i) = evloc_master(i)%ynorth
    vcl(2, i) = evloc_master(i)%xeast
    vcl(3, i) = evloc_master(i)%depth
    vm(i) = i
  enddo
  call dtris3(nev_master, sizht, bf_max, fc_max, vcl, vm, bf_num, nfc, nface, ntetra, bf, fc, ht, ierr)
  print *, "nfc = ", nfc, ntetra
  !!check event pairs
  neventpair_max = 1
  do i = 1, nev_master - 1
    neventpair_max = neventpair_max * (nev_master - i)
  enddo
  allocate(eventpair(1 : nev_master, 1 : nev_master))
  eventpair(1 : nev_master, 1 : nev_master) = .false.
  neventpair = 0
  do k = 1, nfc
    do j = 1, 2
      if(fc(j, k) .le. 0) cycle
      do i = j + 1, 3
        if(fc(i, k) .le. 0) cycle
        dist_tmp = sqrt((evloc_master(vm(fc(i, k)))%ynorth - evloc_master(vm(fc(j, k)))%ynorth) ** 2 &
        &             + (evloc_master(vm(fc(i, k)))%xeast - evloc_master(vm(fc(j, k)))%xeast) ** 2 &
        &             + (evloc_master(vm(fc(i, k)))%depth - evloc_master(vm(fc(j, k)))%depth) ** 2)
        if(dist_tmp .gt. ratio_dist_max) cycle
        if(eventpair(vm(fc(j, k)), vm(fc(i, k))) .eqv. .false.) then
          eventpair(vm(fc(j, k)), vm(fc(i, k))) = .true.
          eventpair(vm(fc(i, k)), vm(fc(j, k))) = .true.
          neventpair = neventpair + 1
        endif
      enddo
    enddo
  enddo
  print *, "neventpair_max = ", neventpair_max, " neventpair = ", neventpair

  !!Set up the observation vector and inversion matrix
  allocate(inversion_matrix(1 : nsta_use * neventpair + 1, 1 : nev_master + 1), obsvector(1 : nsta_use * neventpair + 1))
  allocate(inversion_matrix_org(1 : ubound(inversion_matrix, 1), 1 : ubound(inversion_matrix, 2)), &
  &        obsvector_org(1 : ubound(obsvector, 1)))
  allocate(evsourceamp(1 : nev_master))

  inversion_matrix(1 : ubound(inversion_matrix, 1), 1 : ubound(inversion_matrix, 2)) = 0.0_fp
  obsvector(1 : ubound(obsvector, 1)) = 0.0_fp

  eventpair(1 : nev_master, 1 : nev_master) = .false.
  icount = 1
  do kk = 1, nfc
    do jj = 1, 2
      if(fc(jj, kk) .le. 0) cycle
      do ii = jj + 1, 3
        if(fc(ii, kk) .le. 0) cycle
        dist_tmp = sqrt((evloc_master(vm(fc(ii, kk)))%ynorth - evloc_master(vm(fc(jj, kk)))%ynorth) ** 2 &
        &             + (evloc_master(vm(fc(ii, kk)))%xeast  - evloc_master(vm(fc(jj, kk)))%xeast)  ** 2 &
        &             + (evloc_master(vm(fc(ii, kk)))%depth  - evloc_master(vm(fc(jj, kk)))%depth)  ** 2)
        if(dist_tmp .gt. ratio_dist_max) cycle
        if(eventpair(vm(fc(ii, kk)), vm(fc(jj, kk))) .eqv. .true.) cycle
        eventpair(vm(fc(ii, kk)), vm(fc(jj, kk))) = .true.
        eventpair(vm(fc(jj, kk)), vm(fc(ii, kk))) = .true.
        do i = 1, nsta
          normal_vector(1 : 3) = [sin(ray_azinc(2, i, vm(fc(jj, kk)))) * cos(ray_azinc(1, i, vm(fc(jj, kk)))), &  !!(ynorth)
          &                       sin(ray_azinc(2, i, vm(fc(jj, kk)))) * sin(ray_azinc(1, i, vm(fc(jj, kk)))), &  !!(xeast)
          &                       cos(ray_azinc(2, i, vm(fc(jj, kk))))]                               !!depth (down+)
          matrix_const = (evloc_master(vm(fc(jj, kk)))%ynorth - evloc_master(vm(fc(ii, kk)))%ynorth) * normal_vector(1) &
          &            + (evloc_master(vm(fc(jj, kk)))%xeast  - evloc_master(vm(fc(ii, kk)))%xeast ) * normal_vector(2) &
          &            + (evloc_master(vm(fc(jj, kk)))%depth  - evloc_master(vm(fc(ii, kk)))%depth ) * normal_vector(3)

          inversion_matrix(icount, vm(fc(ii, kk))) =  1.0_fp
          inversion_matrix(icount, vm(fc(jj, kk))) = -1.0_fp
          mean_vel = (evloc_master(vm(fc(jj, kk)))%vel + evloc_master(vm(fc(ii, kk)))%vel) * 0.5_fp
          inversion_matrix(icount, nev_master + 1) = -matrix_const * pi * freq / mean_vel
          obsvector(icount) = log(obsamp_master(i, vm(fc(ii, kk))) / obsamp_master(i, vm(fc(jj, kk)))) &
          &                 + matrix_const / hypodist(i, vm(fc(jj, kk))) 
          icount = icount + 1
        enddo
      enddo
    enddo
  enddo
  !!constraint
  obsvector(icount) = 0.0_fp
  inversion_matrix(icount, ref_evindex) = 1.0_fp

  !!copy observation vector and inversion matrix
  obsvector_org(1 : ubound(obsvector, 1)) = obsvector(1 : ubound(obsvector, 1))
  inversion_matrix_org(1 : ubound(inversion_matrix, 1), 1 : ubound(inversion_matrix, 2)) &
  &  = inversion_matrix(1 : ubound(inversion_matrix, 1), 1 : ubound(inversion_matrix, 2))

  !!calculate least-squares solution
#if defined (MKL)
  call gels(inversion_matrix(:, :) , obsvector(:))
#else
  call la_gels(inversion_matrix(:, :), obsvector(:))
#endif

  !!calculate residual and its variance
  mean_residual = 0.0_fp
  do i = 1, nsta_use * neventpair
    mean_residual = mean_residual + obsvector_org(i) &
    &             - dot_product(inversion_matrix_org(i, 1 : nev_master + 1), obsvector(1 : nev_master + 1))
  enddo
  mean_residual = mean_residual / real(nsta_use * neventpair, kind = fp)
  sigma_residual = 0.0_fp
  do i = 1, nsta_use * neventpair
    sigma_residual = sigma_residual &
    &  + (obsvector_org(i) &
    &  - dot_product(inversion_matrix_org(i, 1 : nev_master + 1), obsvector(1 : nev_master + 1)) &
    &  - mean_residual) ** 2
  enddo
  sigma_residual = sigma_residual / real(nsta_use * neventpair - 1, kind = fp)
   
  allocate(error_matrix(1 : nev_master + 1, 1 : nev_master + 1))
  allocate(ipiv(1 : ubound(error_matrix, 1)))
  error_matrix(:, :) = matmul(transpose(inversion_matrix_org), inversion_matrix_org)
#if defined (MKL)
  call getrf(error_matrix, ipiv)
  call getri(error_matrix, ipiv)
#else
  call la_getrf(error_matrix, ipiv)
  call la_getri(error_matrix, ipiv)
#endif

  do i = 1, nev_master
    evsourceamp(i)%sourceamp = exp(obsvector(i))
    evsourceamp(i)%sigma_sourceamp = sigma_residual * error_matrix(i, i)
  enddo
  attenuationcoef = obsvector(nev_master + 1)
  sigma_attenuationcoef = sigma_residual * error_matrix(nev_master + 1, nev_master + 1)

  !!output result
  outfile = trim(resultfile) // "_sourceamp.txt"
  open(unit = 10, file = outfile)
  write(10, '(a, 2(e15.7, 1x))') "Attenuation param = ", attenuationcoef, sigma_attenuationcoef
  do i = 1, nev_master
    write(10, '(i0, 2(1x, e15.7))') i, evsourceamp(i)%sourceamp, evsourceamp(i)%sigma_sourceamp
  enddo
  close(10)

  deallocate(error_matrix, ipiv)
  stop
end program AmplitudeSourceLocation_masterevent_sourceamp

