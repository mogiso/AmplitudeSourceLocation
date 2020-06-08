!!constant.F90
!!modules contains physical constants
!!Author: Masashi Ogiso (masashi.ogiso@gmail.com)
!!Copyright: (c) Masashi Ogiso 2020
!!License: MIT License https://opensource.org/licenses/MIT


module constants
  use nrtype, only : fp, dp
  implicit none
  private

  real(fp), public, parameter :: pi = 4.0_fp * atan(1.0_fp)
  real(fp), public, parameter :: deg2rad = pi / 180.0_fp, rad2deg = 180.0_fp / pi
  real(fp), public, parameter :: r_earth = 6371.0_fp                             !!averaged earth radius (km)
  real(dp), public, parameter :: r1 = 6378137.0_dp                               !!major axis length (m) of GRS80 ellipsoid
  real(dp), public, parameter :: ec = 1.0_dp / 298.257222101_dp                   !!eccentricity of GRS80 ellipsoid
  real(dp), public, parameter :: r2 = r1 * (1.0_dp - ec)                         !!minor axis length (m) of GRS80 ellipsoid
  real(fp), public, parameter :: r2r1 = real(r2 / r1, kind = fp)


end module constants

