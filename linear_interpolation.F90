module linear_interpolation
  !!provide linear interpolation function (1D/2D/3D) and block interpolation function
  use nrtype, only : fp, dp, sp

  implicit none
  private

  public :: linear_interpolation_1d, linear_interpolation_dvdz_1d, &
  &         linear_interpolation_2d, linear_interpolation_dxdy_2d, &
  &         linear_interpolation_3d, linear_interpolation_dxdydz_3d, &
  &         block_interpolation_2d, block_interpolation_3d

contains

!!1D linear interpolation
!!value_interpolate = value(1) + (x - xgrid(1)) / (xgrid(2) - xgrid(1)) * (value(2) - value(1))
  subroutine linear_interpolation_1d(x, xgrid, value, value_interpolate)
    real(fp), intent(in) :: x, xgrid(2), value(2)
    real(fp), intent(out) :: value_interpolate

    if(xgrid(2) .ne. xgrid(1)) then
      value_interpolate = value(1) + (x - xgrid(1)) / (xgrid(2) - xgrid(1)) * (value(2) - value(1))
    else
      value_interpolate = value(1)
    endif

    return
  end subroutine linear_interpolation_1d

!!1D linear interpolation (derivative)
  subroutine linear_interpolation_dvdz_1d(z, zgrid, v, dvdz)
    real(fp), intent(in) :: z, zgrid(3), v(3)
    real(fp), intent(out) :: dvdz
  
    if(z .ne. zgrid(2)) then
      !!z locates between zgrid(1)とzgrid(2)
      dvdz = (v(2) - v(1)) / (zgrid(2) - zgrid(1))
    else
      !!z locates zgrid(2)
      dvdz = 0.5_fp * ((v(2) - v(1)) / (zgrid(2) - zgrid(1)) + (v(3) - v(2)) / (zgrid(3) - zgrid(2)))
    endif
    return
  end subroutine linear_interpolation_dvdz_1d

!!2D linear interpolation
!!value_interpolate = sigma_x sigma_y (value(xgrid(i), ygrid(j)) * (1 - abs((x - xgrid(i))/(xgrid(2) - xgrid(1))))
!!                                                               * (1 - abs((y - ygrid(j))/(ygrid(2) - ygrid(1))))
  subroutine linear_interpolation_2d(x, y, xgrid, ygrid, value, value_interpolate)
    !!input: (x, y): location where interpolated value are required
    !!      xgrid(2), ygrid(2): location of adjacent grid
    !!      value(2, 2)       : values at grid point
    !!output: value_interpolate: interpolated value
    !!left-hand coordinate system, (x|y)grid(2)>(x|y)grid(1)

    real(fp), intent(in) :: x, y
    real(fp), intent(in) :: xgrid(2), ygrid(2), value(2, 2)
    real(fp), intent(out) :: value_interpolate

    integer :: i, j

    value_interpolate = 0.0_dp
    if(xgrid(1) .ne. xgrid(2) .and. ygrid(1) .ne. ygrid(2)) then
      do j = 1, 2
        do i = 1, 2
          value_interpolate = value_interpolate + value(i, j) &
          &  * (1.0_fp - abs((x - xgrid(i))/(xgrid(2) - xgrid(1)))) &
          &  * (1.0_fp - abs((y - ygrid(j))/(ygrid(2) - ygrid(1))))
          !print *, x, y, xgrid(i), ygrid(j), value_interpolate
        enddo
      enddo
    else
      value_interpolate = value(1, 1)
    endif

    return
  end subroutine linear_interpolation_2d

!!2D block interpolation
  subroutine block_interpolation_2d(x, y, xgrid, ygrid, value, value_interpolate)
    !!input: (x, y): location where interpolated value required
    !!      xgrid(2), ygrid(2): locations of adjacent grid points
    !!      value(2, 2)              : values at grid points
    !!output: value_interpolate: interpolated value at (x, y)

    real(fp), intent(in) :: x, y
    real(fp), intent(in) :: xgrid(2), ygrid(2), value(2, 2)
    real(fp), intent(out) :: value_interpolate

    integer :: i, j

    value_interpolate = 0.0_dp
    if((x - xgrid(1)) .le. (xgrid(2) - x)) then
      i = 1
    else
      i = 2
    endif
    if((y - ygrid(1)) .le. (ygrid(2) - y)) then
      j = 1
    else
      j = 2
    endif
    value_interpolate = value(i, j)

    return
  end subroutine block_interpolation_2d
  
!!2D linear interpolation (derivative)
  subroutine linear_interpolation_dxdy_2d(x, y, xgrid, ygrid, value, dvdx, dvdy)
    real(fp), intent(in) :: x, y
    real(fp), intent(in) :: xgrid(2), ygrid(2), value(2, 2)
    real(fp), intent(out) :: dvdx, dvdy

    integer :: i

    dvdx = 0.0_fp
    dvdy = 0.0_fp
    do i = 1, 2
      dvdx = dvdx &
      &  + value(1, i) * (-1.0_fp) / (xgrid(2) - xgrid(1)) * (1.0_fp - abs(y - ygrid(i)) / (ygrid(2) - ygrid(1)))
      dvdy = dvdy &
      &  + value(i, 1) * (-1.0_fp) / (ygrid(2) - ygrid(1)) * (1.0_fp - abs(x - xgrid(i)) / (xgrid(2) - xgrid(1)))
    enddo
    do i = 1, 2
      dvdx = dvdx + value(2, i) / (xgrid(2) - xgrid(1)) * (1.0_fp - abs(y - ygrid(i)) / (ygrid(2) - ygrid(1)))
      dvdy = dvdy + value(i, 2) / (ygrid(2) - ygrid(1)) * (1.0_fp - abs(x - xgrid(i)) / (xgrid(2) - xgrid(1)))
    enddo

    return
  end subroutine linear_interpolation_dxdy_2d
  
  
!!3D linear interpolation
!!value_interpolate = sigma_x sigma_y sigma_z (value(xgrid(i), ygrid(j), zgrid(k)) * (1 - abs((x - xgrid(i))/(xgrid(2) - xgrid(1))))
!!                                                                                 * (1 - abs((y - ygrid(j))/(ygrid(2) - ygrid(1))))
!!                                                                                 * (1 - abs((z - zgrid(k))/(zgrid(2) - zgrid(1)))))
  subroutine linear_interpolation_3d(x, y, z, xgrid, ygrid, zgrid, value, value_interpolate)
    !!input (x, y, z): location where interpolated value required
    !!      xgrid(2), ygrid(2), zgrid(2): locations of adjacent grid points
    !!      value(2, 2, 2)              : values at grid points
    !!output: value_interpolate: interpolated value at (x, y, z)
    !!left-hand coordinate system, (x|y|z)grid(2)>(x|y|z)grid(1)
    real(fp), intent(in) :: x, y, z
    real(fp), intent(in) :: xgrid(2), ygrid(2), zgrid(2), value(2, 2, 2)
    real(fp), intent(out) :: value_interpolate

    integer :: i, j, k

    value_interpolate = 0.0_fp
    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          value_interpolate = value_interpolate + value(i, j, k) &
          &  * (1.0_fp - abs((x - xgrid(i))/(xgrid(2) - xgrid(1)))) &
          &  * (1.0_fp - abs((y - ygrid(j))/(ygrid(2) - ygrid(1)))) &
          &  * (1.0_fp - abs((z - zgrid(k))/(zgrid(2) - zgrid(1))))
          !print *, x, y, z, xgrid(i), ygrid(j), zgrid(k), value_interpolate
        enddo
      enddo
    enddo

    return
  end subroutine linear_interpolation_3d
  
!!3D block interpolation
  subroutine block_interpolation_3d(x, y, z, xgrid, ygrid, zgrid, value, value_interpolate)
    !!input: (x, y, z): location where interpolated value required
    !!      xgrid(2), ygrid(2), zgrid(2): locations of adjacent grid points
    !!      value(2, 2, 2)              : values at grid points
    !!output: value_interpolate: interpolated value at (x, y, z)
    !!left-hand coordinate system, (x|y|z)grid(2)>(x|y|z)grid(1)

    real(fp), intent(in) :: x, y, z
    real(fp), intent(in) :: xgrid(2), ygrid(2), zgrid(2), value(2, 2, 2)
    real(fp), intent(out) :: value_interpolate

    integer :: i, j, k

    value_interpolate = 0.0_fp
    if((x - xgrid(1)) .le. (xgrid(2) - x)) then
      i = 1
    else
      i = 2
    endif
    if((y - ygrid(1)) .le. (ygrid(2) - y)) then
      j = 1
    else
      j = 2
    endif
    if((z - zgrid(1)) .le. (zgrid(2) - z)) then
      k = 1
    else
      k = 2
    endif
    value_interpolate = value(i, j, k)

    return
  end subroutine block_interpolation_3d
  
!!3D linear interpolation (derivative)
  subroutine linear_interpolation_dxdydz_3d(x, y, z, xgrid, ygrid, zgrid, value, dvdx, dvdy, dvdz)
    !!spatial derivatives of interpolated value
    !!value_interpolate = sigma_x sigma_y sigma_z (value(xgrid(i), ygrid(j), zgrid(k)) * (1 - abs((x - xgrid(i))/(xgrid(2) - xgrid(1))))
    !!                                                                                 * (1 - abs((y - ygrid(j))/(ygrid(2) - ygrid(1))))
    !!                                                                                 * (1 - abs((z - zgrid(k))/(zgrid(2) - zgrid(1)))))
    !!calculate d(value_interpolate)/dx, d(value_interpolate)/dy, d(value_interpolate)/dz

    real(fp), intent(in) :: x, y, z
    real(fp), intent(in) :: xgrid(2), ygrid(2), zgrid(2), value(2, 2, 2)
    real(fp), intent(out) :: dvdx, dvdy, dvdz

    integer :: j, k

    dvdx = 0.0_dp
    dvdy = 0.0_dp
    dvdz = 0.0_dp
    do k = 1, 2
      do j = 1, 2
        dvdx = dvdx + value(1, j, k) &
        &  * (1.0_fp - abs((y - ygrid(j)) / (ygrid(2) - ygrid(1)))) &
        &  * (1.0_fp - abs((z - zgrid(k)) / (zgrid(2) - zgrid(1)))) &
        &  * (-1.0_fp / (xgrid(2) - xgrid(1)))
        dvdy = dvdy + value(j, 1, k) &
        &  * (1.0_fp - abs((x - xgrid(j)) / (xgrid(2) - xgrid(1)))) &
        &  * (1.0_fp - abs((z - zgrid(k)) / (zgrid(2) - zgrid(1)))) &
        &  * (-1.0_fp / (ygrid(2) - ygrid(1)))
        dvdz = dvdz + value(j, k, 1) &
        &  * (1.0_fp - abs((x - xgrid(j)) / (xgrid(2) - xgrid(1)))) &
        &  * (1.0_fp - abs((y - ygrid(k)) / (ygrid(2) - ygrid(1)))) &
        &  * (-1.0_fp / (zgrid(2) - zgrid(1)))
      enddo
    enddo
    do k = 1, 2
      do j = 1, 2
        dvdx = dvdx + value(2, j, k) &
        &  * (1.0_fp - abs((y - ygrid(j)) / (ygrid(2) - ygrid(1)))) &
        &  * (1.0_fp - abs((z - zgrid(k)) / (zgrid(2) - zgrid(1)))) &
        &  * (1.0_fp / (xgrid(2) - xgrid(1)))
        dvdy = dvdy + value(j, 2, k) &
        &  * (1.0_fp - abs((x - xgrid(j)) / (xgrid(2) - xgrid(1)))) &
        &  * (1.0_fp - abs((z - zgrid(k)) / (zgrid(2) - zgrid(1)))) &
        &  * (1.0_fp / (ygrid(2) - ygrid(1)))
        dvdz = dvdz + value(j, k, 2) &
        &  * (1.0_fp - abs((x - xgrid(j)) / (xgrid(2) - xgrid(1)))) &
        &  * (1.0_fp - abs((y - ygrid(k)) / (ygrid(2) - ygrid(1)))) &
        &  * (1.0_fp / (zgrid(2) - zgrid(1)))
      enddo
    enddo

    return
  end subroutine linear_interpolation_dxdydz_3d

end module linear_interpolation

