module def_gridpoint
  use nrtype, only : fp

  type gridpoint
    real(kind = fp) :: lon, lat, dep, r, theta, phi 
                    !!longitude, latitude, depth, distance from earth center, colatitude (rad), azimuth (rad)
  end type gridpoint

end module def_gridpoint
