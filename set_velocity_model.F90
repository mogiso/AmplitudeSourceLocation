module set_velocity_model
  private
  public :: set_velocity

  contains
 
  subroutine set_velocity(z_min, dz, velocity, qinv)
    use nrtype, only : fp
    implicit none
    real(kind = fp), intent(in)  :: z_min, dz                  !!size of depth direction
    real(kind = fp), intent(out) :: velocity(:, :, :), qinv(:, :, :)

    real(kind = fp), paramter    :: VpVs = sqrt(3.0_fp)
    integer :: i, nlon, nlat, nz
    real(kind = fp) :: depth

    nlon = ubound(velocity, 1)
    nlat = ubound(velocity, 2)
    nz = ubound(velocity, 3)
  
    do i = 1, nz
      depth = z_min + dz * real(i - 1, kind = fp)
#ifdef MEA_1D
      if(i .eq. 1) write(0, '(a)') "velocity/Qinv = meakan 1d structure"
      if(depth .lt. 0.5_fp) then
        velocity(1 : nlon, 1 : nlat, i) = (3.2_fp - 3.0_fp) / (0.5_fp - 0.0_fp) * depth + 3.0_fp
      elseif(depth .ge. 0.5_fp .and depth .lt. 1.0_fp) then
        velocity(1 : nlon, 1 : nlat, i) = (4.5_fp - 3.2_fp) / (1.0_fp - 0.5_fp) * (depth - 0.5_fp) + 3.2_fp
      elseif(depth .ge. 1.0_fp .and depth .lt. 1.5_fp) then
        velocity(1 : nlon, 1 : nlat, i) = (5.0_fp - 4.5_fp) / (1.5_fp - 1.0_fp) * (depth - 1.0_fp) + 4.5_fp
      elseif(depth .ge. 1.5_fp .and depth .lt. 2.0_fp) then
        velocity(1 : nlon, 1 : nlat, i) = (5.133_fp - 5.0_fp) / (2.0_fp - 1.5_fp) * (depth - 1.5_fp) + 5.0_fp
      elseif(depth .ge. 2.0_fp .and depth .lt. 2.5_fp) then
        velocity(1 : nlon, 1 : nlat, i) = (5.262_fp - 5.133_fp) / (2.5_fp - 2.0_fp) * (depth - 2.0_fp) + 5.133_fp
      elseif(depth .ge. 2.5_fp .and depth .lt. 3.0_fp) then
        velocity(1 : nlon, 1 : nlat, i) = (5.383_fp - 5.262_fp) / (3.0_fp - 2.5_fp) * (depth - 2.5_fp) + 5.262_fp
      elseif(depth .ge. 3.0_fp) then
        velocity(1 : nlon, 1 : nlat, i) = (5.494_fp - 5.383_fp) / (3.5_fp - 3.0_fp) * (depth - 3.0_fp) + 5.383_fp
      endif

      !!the values of Qinv are taken from Kumagai et al (2019): case of Nevado
      !del Ruiz volcano (layer thickness is different)
      if(depth .lt. 1.0_fp) then
        qinv(1 : nlon, 1 : nlat, i) = 1.0_fp / 40.0_fp
      else
        qinv(1 : nlon, 1 : nlat, i) = 1.0_fp / 180.0_fp
      endif
#else
      if(i .eq. 1) write(0, '(a)') "Velocity/Qinv = constant"
      velocity(1 : nlon, 1 : nlat, i) = 2.5_fp
      qinv(1 : nlon, 1 : nlat, i) = 1.0_fp / 50.0_fp
#endif

      !!Vp to Vs
      velocity(1 : nlon, 1 : nlat, i) = velocity(1 : nlon, 1 : nlat, i) / VpVs

    enddo
    return
  end subroutine set_velocity

end module set_velocity_model



