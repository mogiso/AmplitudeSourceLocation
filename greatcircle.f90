module greatcircle
  private
  public :: greatcircle_dist, latctog, latgtoc

  contains

  subroutine greatcircle_dist(lat0, lon0, lat1, lon1, delta_out, distance, azimuth, backazimuth)
    use nrtype, only : fp
    use constants, only : pi, r_earth, r2r1, deg2rad, rad2deg
    implicit none

    !!Calculate distance, azimuth, backazimuth using Vincenty (1975, Survey review)
    !!Input: (lat0, lon0) (point A) , (lat1, lon1) (point B)
    !!output: delta(rad), dist(km), az(rad, point A -> point B), baz(rad, point B -> point A)
  
    real(kind = fp), intent(in) :: lat0, lon0, lat1, lon1
    real(kind = fp), intent(out), optional :: delta_out, distance, azimuth, backazimuth
  
    real(kind = fp) :: c, b, delta, lat_0, lat_1, sin_lat0, cos_lat0, sin_lat1, cos_lat1, lon_diff, &
    &                  cos_ldiff, sin_ldiff
  
    !!geographical latitude to geocentral latitude
    call latgtoc(lat0, lat_0)
    call latgtoc(lat1, lat_1)

    sin_lat0 = sin(lat_0)
    cos_lat0 = cos(lat_0)
    sin_lat1 = sin(lat_1)
    cos_lat1 = cos(lat_1)
    lon_diff = (lon1 - lon0) * deg2rad
    cos_ldiff = cos(lon_diff)
    sin_ldiff = sin(lon_diff)
    delta = atan2(sqrt(cos_lat1 * sin_ldiff * cos_lat1 * sin_ldiff &
    &                  + (cos_lat0 * sin_lat1 - sin_lat0 * cos_lat1 * cos_ldiff) &
    &                  * (cos_lat0 * sin_lat1 - sin_lat0 * cos_lat1 * cos_ldiff)), &
    &                  (sin_lat0 * sin_lat1 + cos_lat0 * cos_lat1 * cos_ldiff))

    c = atan2(cos_lat1 * sin_ldiff, cos_lat0 * sin_lat1 - sin_lat0 * cos_lat1 * cos_ldiff)
    if(c .lt. 0.0_fp) c = c + 2.0_fp * pi
    b = pi + atan2(cos_lat0 * sin_ldiff, -sin_lat0 * cos_lat1 + cos_lat0 * sin_lat1 * cos_ldiff)
    if(b .ge. 2.0_fp * pi) b = b - 2.0_fp * pi
  
    if(present(delta_out)) delta_out = delta
    if(present(distance)) distance = r_earth * delta
    if(present(azimuth)) azimuth = c
    if(present(backazimuth)) backazimuth = b
  
    return
  end subroutine greatcircle_dist

  subroutine latgtoc(lat_in, lat_out)
    use nrtype, only : fp
    use constants, only : r2r1, deg2rad
    real(kind = fp), intent(in) :: lat_in    !!degree
    real(kind = fp), intent(out) :: lat_out  !!radian

    lat_out = atan(r2r1 * r2r1 * tan(lat_in * deg2rad))
    return
  end subroutine latgtoc

  subroutine latctog(lat_in, lat_out)
    use nrtype, only : fp
    use constants, only : r2r1, deg2rad, rad2deg
    real(kind = fp), intent(in) :: lat_in    !!radian
    real(kind = fp), intent(out) :: lat_out  !!degree

    lat_out = atan(tan(lat_in) / (r2r1 * r2r1)) * rad2deg
    return
  end subroutine latctog


end module greatcircle
