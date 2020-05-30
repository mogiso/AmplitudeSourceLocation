module set_velocity_model
  private
  public :: set_velocity

  contains
 
  subroutine set_velocity(dz, velocity)
    use nrtype, only : sp, fp
    implicit none
    real(kind = fp), intent(in) :: dz                  !!size of depth direction
    real(kind = fp), intent(out) :: velocity(:, :, :, :)

    integer :: i, nlon, nlat, nz
    real(kind = fp) :: depth
#ifdef EJIVSM
    integer :: j, k, irec
    real(kind = sp) :: lon_tmp, lat_tmp, dep_tmp, vel_tmp
#endif

    nlon = ubound(velocity, 1)
    nlat = ubound(velocity, 2)
    nz = ubound(velocity, 3)
  
#ifdef EJIVSM
    write(0, '(a)') "velocity = eJIVSM"
    open(unit = 30, file = "velocity_ejivsm.dat", form = "unformatted", access = "direct", recl = 4 * 4)
    irec = 1
    do k = 1, nz
      do j = 1, nlat
        do i = 1, nlon
          read(30, rec = irec) lon_tmp, lat_tmp, dep_tmp, vel_tmp
          irec = irec + 1
          velocity(i, j, k, 2) = real(vel_tmp, kind = fp)
          velocity(i, j, k, 1) = velocity(i, j, k, 2) * sqrt(3.0_fp)
        enddo
      enddo
    enddo
    close(30)
#else
    do i = 1, nz
      depth = dz * real(i - 1, kind = fp)
#if defined(HOKUDAI)
      if(i .eq. 0) write(0, '(a)') "velocity = HOKUDAI routine"
      !!1D structure used by Hokkaido Univ routine
      if(depth .le. 30.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 2) &
        &  = ((6.8_fp - 5.8_fp) / (30.0_fp - 0.0_fp) * depth + 5.8_fp) / sqrt(3.0_fp)
      elseif(depth .gt. 30.0_fp .and. depth .le. 35.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 2) &
        &  = ((7.6_fp - 6.8_fp) / (35.0_fp - 30.0_fp) * (depth - 30.0_fp) + 6.8_fp) / sqrt(3.0_fp)
      elseif(depth .gt. 35.0_fp) then
        velocity(1 : nlon, 1 : nlat1, i, 2) &
        &  = ((8.0_fp - 7.6_fp) / (200.0_fp - 35.0_fp) * (depth - 35.0_fp) + 7.6_fp) / sqrt(3.0_fp)
      endif
      velocity(1 : nlon, 1 : nlat, i, 1) = velocity(1 : nlon, 1 : nlat, i, 2) * sqrt(3.0_fp)
#elif defined(JMA2001)
      !!1D structure similar to the JMA2001 velocity model
      if(i .eq. 0) write(0, '(a)') "velocity = JMA2001"
      if(depth .le. 4.5_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = ((5.780_fp - 4.800_fp) / 4.5_fp) * depth + 4.800_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = ((3.409_fp - 2.844_fp) / 4.5_fp) * depth + 2.844_fp
      elseif(depth .gt. 4.5_fp .and. depth .le. 21.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = &
        &  ((6.502_fp - 5.780_fp) / (21.0_fp - 4.5_fp)) * (depth - 4.5_fp) + 5.780_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = &
        &  ((3.782_fp - 3.409_fp) / (21.0_fp - 4.5_fp)) * (depth - 4.5_fp) + 3.409_fp
      elseif(depth .gt. 21.0_fp .and. depth .le. 45.5_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = &
        &  ((7.638_fp - 6.502_fp) / (45.5_fp - 21.0_fp)) * (depth - 21.0_fp) + 6.502_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = &
        &  ((4.362_fp - 3.782_fp) / (45.5_fp - 21.0_fp)) * (depth - 21.0_fp) + 3.782_fp
      elseif(depth .gt. 45.5_fp .and. depth .le. 83.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = &
        &  ((7.900_fp - 7.638_fp) / (83.0_fp - 45.5_fp)) * (depth - 45.5_fp) + 7.638_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = &
        &  ((4.424_fp - 4.362_fp) / (83.0_fp - 45.5_fp)) * (depth - 45.5_fp) + 4.362_fp
      elseif(depth .gt. 83.0_fp .and. depth .le. 175.5_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = &
        &  ((8.184_fp - 7.900_fp) / (175.5_fp - 83.0_fp)) * (depth - 83.0_fp) + 7.900_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = &
        &  ((4.558_fp - 4.424_fp) / (175.5_fp - 83.0_fp)) * (depth - 83.0_fp) + 4.424_fp
      elseif(depth .gt. 175.5_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = &
        &  ((8.269_fp - 8.184_fp) / (200.0_fp - 175.5_fp)) * (depth - 175.5_fp) + 8.184_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = &
        &  ((4.602_fp - 4.558_fp) / (200.0_fp - 175.5_fp)) * (depth - 175.5_fp) + 4.558_fp
      endif
#elif defined (UKAWA)
      !!1D structure of Kanto region deribed by Ukawa et al. (1984, NIED research report)
      if(i .eq. 0) write(0, '(a)') "velocity = UKAWA"
      if(depth .le. 10.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = ((5.98_fp - 5.50_fp) / 10.0_fp) * depth + 5.50_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = ((3.49_fp - 3.25_fp) / 10.0_fp) * depth + 3.25_fp
      elseif(depth .gt. 10.0_fp .and. depth .le. 20.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = ((6.51_fp - 5.98_fp) / (20.0_fp - 10.0_fp)) * (depth - 10.0_fp) + 5.98_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = ((3.74_fp - 3.49_fp) / (20.0_fp - 10.0_fp)) * (depth - 10.0_fp) + 3.49_fp
      elseif(depth .gt. 20.0_fp .and. depth .lt. 32.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = ((7.20_fp - 6.51_fp) / (32.0_fp - 20.0_fp)) * (depth - 20.0_fp) + 6.51_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = ((4.07_fp - 3.74_fp) / (32.0_fp - 20.0_fp)) * (depth - 20.0_fp) + 3.74_fp
      elseif(depth .eq. 32.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = 7.50_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = 4.24_fp
      elseif(depth .gt. 32.0_fp .and. depth .le. 50.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = ((7.85_fp - 7.80_fp) / (50.0_fp - 32.0_fp)) * (depth - 32.0_fp) + 7.80_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = ((4.43_fp - 4.41_fp) / (50.0_fp - 32.0_fp)) * (depth - 32.0_fp) + 4.41_fp
      elseif(depth .gt. 50.0_fp .and. depth .le. 100.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = ((8.00_fp - 7.85_fp) / (100.0_fp - 50.0_fp)) * (depth - 50.0_fp) + 7.85_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = ((4.50_fp - 4.43_fp) / (100.0_fp - 50.0_fp)) * (depth - 50.0_fp) + 4.43_fp
      elseif(depth .gt. 100.0_fp .and. depth .le. 150.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = ((8.12_fp - 8.00_fp) / (150.0_fp - 100.0_fp)) * (depth - 100.0_fp) + 8.0_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = ((4.53_fp - 4.50_fp) / (150.0_fp - 100.0_fp)) * (depth - 100.0_fp) + 4.5_fp
      elseif(depth .gt. 150.0_fp .and. depth .le. 200.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = ((8.26_fp - 8.12_fp) / (200.0_fp - 150.0_fp)) * (depth - 150.0_fp) + 8.12_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = ((4.60_fp - 4.53_fp) / (200.0_fp - 150.0_fp)) * (depth - 150.0_fp) + 4.53_fp
      elseif(depth .gt. 200.0_fp .and. depth .le. 250.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = ((8.42_fp - 8.26_fp) / (250.0_fp - 200.0_fp)) * (depth - 200.0_fp) + 8.26_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = ((4.76_fp - 4.60_fp) / (250.0_fp - 200.0_fp)) * (depth - 200.0_fp) + 4.60_fp
      elseif(depth .gt. 250.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = ((8.58_fp - 8.42_fp) / (300.0_fp - 250.0_fp)) * (depth - 250.0_fp) + 8.42_fp
        velocity(1 : nlon, 1 : nlat, i, 2) = ((4.76_fp - 4.68_fp) / (300.0_fp - 250.0_fp)) * (depth - 250.0_fp) + 4.68_fp
      endif
#else
      if(i .eq. 0) write(0, '(a)') "velocity = CONST"
      velocity(1 : nlon, 1 : nlat, i, 2) = 3.5_fp
      velocity(1 : nlon, 1 : nlat, i, 1) = velocity(1 : nlon, 1 : nlat, i, 2) * sqrt(3.0_fp)
#endif
    enddo
#endif
    return
  end subroutine set_velocity
end module set_velocity_model



