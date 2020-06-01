module set_velocity_model
  private
  public :: set_velocity

  contains
 
  subroutine set_velocity(z_min, dz, velocity, qinv)
    use nrtype, only : fp
    implicit none
    real(kind = fp), intent(in)  :: z_min, dz                  !!size of depth direction
    real(kind = fp), intent(out) :: velocity(:, :, :), qinv(:, :, :)

    integer :: i, nlon, nlat, nz
    real(kind = fp) :: depth

    nlon = ubound(velocity, 1)
    nlat = ubound(velocity, 2)
    nz = ubound(velocity, 3)
  
    do i = 1, nz
      depth = z_min + dz * real(i - 1, kind = fp)
#ifdef VEL_1D
      if(i .eq. 0) write(0, '(a)') "velocity = 1d structure"
      if(depth .lt. 0.0_fp) then
        velocity(1 : nlon, 1 : nlat, i) = (3.0_fp - 1.5_fp) / (0.0_fp + 1.5_fp) * (depth + 1.5_fp) + 1.5_fp
      else
        velocity(1 : nlon, 1 : nlat, i) = 3.0_fp
      endif
#else
      if(i .eq. 0) write(0, '(a)') "velocity = CONST"
      velocity(1 : nlon, 1 : nlat, i) = 2.5_fp
#endif

#ifdef QINV_1D
      if(i .eq. 0) write(0, '(a)') "qinv = 1d_structure"
      if(depth .lt. 0.0_fp) then
        qinv(1 : nlon, 1 : nlat, i) = 0.05_fp
      else
        qinv(1 : nlon, 1 : nlat, i) = 0.01_fp
      endif
#else
      if(i .eq. 0) write(0, '(a)') "qinv = CONST"
      qinv(1 : nlon, 1 : nlat, i) = 0.02_fp
#endif

    enddo
    return
  end subroutine set_velocity

end module set_velocity_model



