module set_velocity_model
  !!gives velocity and attenuation structure for AmplitudeSourceLocation_PulseWidth.F90
  !!Author: Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2020
  !!License  : MIT License https://opensource.org/licenses/MIT
  private
  public :: set_velocity

  contains
 
  subroutine set_velocity(z_min, dz, velocity, qinv)
    use nrtype, only : fp
    implicit none
    real(kind = fp), intent(in)  :: z_min, dz                  !!size of depth direction
    real(kind = fp), intent(out) :: velocity(:, :, :, :), qinv(:, :, :, :)

    real(kind = fp), parameter   :: VpVs = sqrt(3.0_fp)
    integer :: i, nlon, nlat, nz
    real(kind = fp) :: depth

    nlon = ubound(velocity, 1)
    nlat = ubound(velocity, 2)
    nz = ubound(velocity, 3)
  
    qinv(1 : nlon, 1 : nlat, 1 : nz, 1 : 2) = 0.0_fp
    velocity(1 : nlon, 1 : nlat, 1 : nz, 1 : 2) = 0.0_fp
    do i = 1, nz
      depth = z_min + dz * real(i - 1, kind = fp)

#if defined (V_MEA1D)
      if(i .eq. 1) write(0, '(a)') "velocity/Qinv = meakan 1d structure"
      if(depth .lt. 0.5_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = (3.2_fp - 3.0_fp) / (0.5_fp - 0.0_fp) * depth + 3.0_fp
      elseif(depth .ge. 0.5_fp .and. depth .lt. 1.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = (4.5_fp - 3.2_fp) / (1.0_fp - 0.5_fp) * (depth - 0.5_fp) + 3.2_fp
      elseif(depth .ge. 1.0_fp .and. depth .lt. 1.5_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = (5.0_fp - 4.5_fp) / (1.5_fp - 1.0_fp) * (depth - 1.0_fp) + 4.5_fp
      elseif(depth .ge. 1.5_fp .and. depth .lt. 2.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = (5.133_fp - 5.0_fp) / (2.0_fp - 1.5_fp) * (depth - 1.5_fp) + 5.0_fp
      elseif(depth .ge. 2.0_fp .and. depth .lt. 2.5_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = (5.262_fp - 5.133_fp) / (2.5_fp - 2.0_fp) * (depth - 2.0_fp) + 5.133_fp
      elseif(depth .ge. 2.5_fp .and. depth .lt. 3.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = (5.383_fp - 5.262_fp) / (3.0_fp - 2.5_fp) * (depth - 2.5_fp) + 5.262_fp
      elseif(depth .ge. 3.0_fp) then
        velocity(1 : nlon, 1 : nlat, i, 1) = (5.494_fp - 5.383_fp) / (3.5_fp - 3.0_fp) * (depth - 3.0_fp) + 5.383_fp
      endif
      !!Vp to Vs
      velocity(1 : nlon, 1 : nlat, i, 2) = velocity(1 : nlon, 1 : nlat, i, 1) / VpVs

      !!the values of Qinv are taken from Kumagai et al (2019): case of Nevado
      !del Ruiz volcano (layer thickness is different)
      if(depth .lt. 1.0_fp) then
        qinv(1 : nlon, 1 : nlat, i, 2) = 1.0_fp / 40.0_fp
      else
        qinv(1 : nlon, 1 : nlat, i, 2) = 1.0_fp / 180.0_fp
      endif
      qinv(1 : nlon, 1 : nlat, i, 1) = qinv(1 : nlon, 1 : nlat, i, 2) / (9.0_fp / 4.0_fp)

#elif defined (V_CONST)
      if(i .eq. 1) write(0, '(a)') "Velocity/Qinv = constant"
      velocity(1 : nlon, 1 : nlat, i, 1) = 2.5_fp
      !!Vp to Vs
      velocity(1 : nlon, 1 : nlat, i, 2) = velocity(1 : nlon, 1 : nlat, i, 1) / VpVs
      qinv(1 : nlon, 1 : nlat, i, 2) = 1.0_fp / 50.0_fp
      qinv(1 : nlon, 1 : nlat, i, 1) = qinv(1 : nlon, 1 : nlat, i, 2) / (9.0_fp / 4.0_fp)

#elif defined (JMA2001)
      if(i .eq. 1) write(0, '(a)') "velocity = (approx.) JMA2001"
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
      qinv(1 : nlon, 1 : nlat, i, 2) = 1.0_fp / 200.0_fp
      qinv(1 : nlon, 1 : nlat, i, 1) = qinv(1 : nlon, 1 : nlat, i, 2) / (9.0_fp / 4.0_fp)

#else
      write(0, *) "please set -DV_CONST or appropriate definition"
      stop
#endif


    enddo
    return
  end subroutine set_velocity

end module set_velocity_model



