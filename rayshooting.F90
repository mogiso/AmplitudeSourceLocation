!!provide ray shooting subroutine (2D/3D)
!!shooting equation: based on Koketsu (1991, Journal of Seismological society of Japan)
!!Last update: 2019-12-04 17:31:15 fix the case in which az_out = 2*pi

module rayshooting

  implicit none
  private
  public :: rayshooting2D, rayshooting3D, dist_to_latlon

contains

  subroutine rayshooting2D(lon_in, lat_in, az_in, dtime_step, velocity, dvdlon, dvdlat, &
  &                        lon_out, lat_out, az_out)
    use nrtype, only : fp
    use constants, only : pi, r_earth, deg2rad, rad2deg
    use greatcircle, only : latctog, latgtoc

    real(fp), intent(in) :: lon_in, lat_in, az_in        !!longitude (deg), latitude(deg), azimuth (rad)
    real(fp), intent(in) :: dtime_step                   !!time step (s), velocity at current location (km/s) and derivatives (lon, lat)
    real(fp), intent(in) :: velocity, dvdlon, dvdlat
    real(fp), intent(out) :: lon_out, lat_out, az_out    !!new longitude (deg), new latitude (deg), new azimuth (rad)

    real(fp) :: loc_lontmp, loc_lattmp
    real(fp) :: delta_lattmp, delta_lontmp, daz

    loc_lontmp = lon_in * deg2rad
    call latgtoc(lat_in, loc_lattmp)
    loc_lattmp = pi * 0.5_fp - loc_lattmp       !!from latitude to co-latitude

    !!Koketsu (1991) eq. (33) with incident ancle = pi/2, depth = earth surface
    delta_lattmp = -velocity / r_earth * (cos(az_in) * dtime_step)
    delta_lontmp = velocity / (r_earth * sin(loc_lattmp)) * sin(az_in) * dtime_step
    daz = (dvdlat / r_earth * sin(az_in) &
    &  - dvdlon / (r_earth * sin(loc_lattmp)) * cos(az_in)) &
    &  + velocity / r_earth * sin(az_in) / tan(loc_lattmp)
    az_out = az_in + daz * dtime_step
    if(az_out .ge. 2.0_fp * pi) az_out = az_out - 2.0_fp * pi
    if(az_out .lt. 0.0_fp) az_out = az_out + 2.0_fp * pi

    lon_out = (loc_lontmp + delta_lontmp) * rad2deg
    if(lon_out .gt. 180.0_fp) lon_out = lon_out - 360.0_fp
    if(lon_out .le. -180.0_fp) lon_out = lon_out + 360.0_fp
    loc_lattmp = 90.0_fp - (loc_lattmp + delta_lattmp) * rad2deg
    call latctog(loc_lattmp, lat_out)

    return
  end subroutine rayshooting2D

  subroutine rayshooting3D(lon_in, lat_in, z_in, az_in, inc_angle_in, dtime_step, velocity, dvdlon, dvdlat, dvdz, &
  &                        lon_out, lat_out, z_out, az_out, inc_angle_out)
    use nrtype, only : fp
    use constants, only : pi, r_earth, deg2rad, rad2deg
    use greatcircle, only : latctog, latgtoc

    real(fp), intent(in) :: lon_in, lat_in, z_in            !!longitude (deg), geographical latitude(deg), depth from surface (km)
    real(fp), intent(in) :: az_in, inc_angle_in             !!azimuth (rad), incident angle measured from depth direction (rad)
    real(fp), intent(in) :: dtime_step                      !!time step (s), velocity at current location (km/s) and derivatives (lon, lat, depth)
    real(fp), intent(in) :: velocity, dvdlon, dvdlat, dvdz
    real(fp), intent(out) :: lon_out, lat_out, z_out        !!new longitude (deg), new geographical latitude (deg), new depth from surface (km)
    real(fp), intent(out) :: az_out, inc_angle_out          !!new azimuth (rad), new incident angle(rad)

    real(fp) :: loc_lontmp, loc_lattmp, loc_ztmp, sin_inc, sin_az, inv_sin_lat, cos_inc, cos_az, inv_loc_ztmp
    real(fp) :: delta_lattmp, delta_lontmp, delta_ztmp, dinc_angle, daz

    loc_lontmp = lon_in * deg2rad
    call latgtoc(lat_in, loc_lattmp)
    loc_lattmp = pi * 0.5_fp - loc_lattmp   !!from latitude to co-latitude
    loc_ztmp = r_earth - z_in               !!depth from surface -> distance from earth center

    sin_inc = sin(inc_angle_in)
    if(abs(sin_inc) .le. 1.0e-5_fp) sin_inc = 1.0e-5_fp * sin_inc / abs(sin_inc)
    sin_az = sin(az_in)
    inv_sin_lat = 1.0_fp / sin(loc_lattmp)
    cos_inc = cos(inc_angle_in)
    cos_az = cos(az_in)
    inv_loc_ztmp =  1.0_fp / loc_ztmp

    !!Koketsu (1991) eq. (33)
    !delta_lattmp = -velocity / loc_ztmp * sin(inc_angle_in) * cos(az_in) * dtime_step
    !delta_lontmp = velocity / (loc_ztmp * sin(loc_lattmp)) * sin(inc_angle_in) * sin(az_in) * dtime_step
    !delta_ztmp = -velocity * dtime_step * cos(inc_angle_in)
    delta_lattmp = -velocity * inv_loc_ztmp * sin_inc     * cos_az  * dtime_step
    delta_lontmp =  velocity * inv_loc_ztmp * inv_sin_lat * sin_inc * sin_az * dtime_step
    delta_ztmp   = -velocity * dtime_step   * cos_inc

    !dinc_angle = (dvdz + velocity / loc_ztmp) * sin(inc_angle_in) &
    !&  - cos(inc_angle_in) * (dvdlat / loc_ztmp * cos(az_in) + dvdlon / (loc_ztmp * sin(loc_lattmp)) * sin(az_in))
    !daz = (dvdlat / loc_ztmp * sin(az_in) - dvdlon / (loc_ztmp * sin(loc_lattmp)) * cos(az_in)) / sin(inc_angle_in) &
    !&  + velocity / loc_ztmp * sin(inc_angle_in) * sin(az_in) / tan(loc_lattmp)
    dinc_angle = (dvdz + velocity * inv_loc_ztmp) * sin_inc &
    &  - cos_inc * (dvdlat * inv_loc_ztmp * cos_az + dvdlon * inv_loc_ztmp * inv_sin_lat * sin_az)
    daz = (dvdlat * inv_loc_ztmp * sin_az - dvdlon * inv_loc_ztmp * inv_sin_lat * cos_az) / sin_inc &
    &  + velocity * inv_loc_ztmp * sin_inc * sin_az * cos(loc_lattmp) * inv_sin_lat

    inc_angle_out = inc_angle_in + dinc_angle * dtime_step
    az_out = az_in + daz * dtime_step
    if(az_out .ge. 2.0_fp * pi) az_out = az_out - 2.0_fp * pi
    if(az_out .lt. 0.0_fp) az_out = az_out + 2.0_fp * pi

    lon_out = (loc_lontmp + delta_lontmp) * rad2deg
    if(lon_out .gt. 180.0_fp) lon_out = lon_out - 360.0_fp
    if(lon_out .le. -180.0_fp) lon_out = lon_out + 360.0_fp
    loc_lattmp = pi * 0.5_fp - (loc_lattmp + delta_lattmp)
    call latctog(loc_lattmp, lat_out)
    z_out = r_earth - (loc_ztmp + delta_ztmp)


    return
  end subroutine rayshooting3D

  subroutine dist_to_latlon(lon_in, lat_in, az_in, dtime_step, velocity, lon_out, lat_out, az_out)
    !!ray trace in homogeneous space, based on spherical trigonometry
    use nrtype, only : fp
    use constants, only : pi, r1, deg2rad, rad2deg
    use greatcircle, only : latctog, latgtoc
    implicit none
    real(kind = fp), intent(in) :: lon_in, lat_in, az_in, dtime_step
    real(kind = fp), intent(in) :: velocity
    real(kind = fp), intent(out) :: lon_out, lat_out, az_out

    real(kind = fp) :: lat_tmp, delta, sin_lat, cos_lat, sin_lon, cos_lon

    call latgtoc(lat_in, lat_tmp)
    delta = (velocity * dtime_step) / (real(r1, kind = fp) * 0.001_fp)
    sin_lat = sin(lat_tmp) * cos(delta) + cos(lat_tmp) * sin(delta) * cos(az_in)
    cos_lat = sqrt(1.0_fp - sin_lat * sin_lat)
    sin_lon = sin(delta) * sin(az_in) / cos_lat
    cos_lon = (sin_lat * cos(lat_tmp) - sin(delta) * cos(az_in)) / cos_lat / sin(lat_tmp)
    lon_out = atan2(sin_lon, cos_lon) * rad2deg + lon_in
    lat_tmp = atan2(sin_lat, cos_lat)
    call latctog(lat_tmp, lat_out)
    az_out = az_in

  end subroutine dist_to_latlon
    
end module rayshooting

