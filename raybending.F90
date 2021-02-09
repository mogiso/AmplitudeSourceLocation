module raybending
  implicit none
  private

contains

  subroutine move_node(node1, node2, node3, velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, enhance_factor)
    !!move node2 based on Koketsu and Sekine (1998) eq. 19
    !!node1 and node3 are fixed
    use def_gridpoint
    use nrtype,               only : fp
    use constants,            only : r_earth, pi, rad2deg, deg2rad
    use greatcircle,          only : latctog
    use linear_interpolation, only : linear_interpolation_3d
    implicit none
    type(gridpoint), intent(in)    :: node1, node3
    type(gridpoint), intent(inout) :: node2
    real(kind = fp), intent(in)    :: velocity(:, :, :)       !!either velocity of  P- or S-waves
    real(kind = fp), intent(in)    :: lon_w, lat_s, dep_min, dlon, dlat, ddep, enhance_factor

    real(kind = fp) :: dist_l, normal_vector_mid, slowness_mid, const, dist_move
    real(kind = fp) :: tangentvector_mid(3), normal_vector_mid(3), v_mid(3), &
    &                  lon_grid(2), lat_grid(2), z_grid(2), val_3d(2, 2, 2)
    integer :: lon_index, lat_index, z_index
    type(gridpoint) :: node_tmp

    !!calculate location of midpoint between node1 and node3
    node_mid%r = 0.5_fp * (node1%r + node3%r)
    node_mid%theta = 0.5_fp * (node1%theta + node3%theta)
    node_mid%phi = 0.5_fp * (node1%phi + node3%phi)
    node_mid%dep = r_earth - node_mid%r
    node_mid%lon = node_mid%phi * rad2deg
    call latctog(pi * 0.5_fp - node_mid%theta, node_mid%lat)

    !!calculate distance between node1 and node3, divide by 2
    call hypodist(node1%lat, node1%lon, node1%dep, node3%lat, node3%lon, node3%dep, dist_l)
    dist_l = dist_l * 0.5_fp

    !!calculate tangent vector at node_mid
    tangentvector_mid(1) = 0.5_fp / dist_l * (node3%r - node1%r)
    tangentvector_mid(2) = 0.5_fp / dist_l * node_mid%r * (node3%theta - node1%theta)
    tangentvector_mid(3) = 0.5_fp / dist_l * node_mid%r * sin(node_mid%theta) * (node3%phi - node1%phi)

    !!calculate velocity at node_mid
    lon_index = int((node_mid%lon - lon_w) / dlon) + 1
    lat_index = int((node_mid%lat - lat_s) / dlat) + 1
    z_index   = int((node_mid%dep - dep_min) / ddep) + 1
    lon_grid(1 : 2) = lon_w + dlon * real(lon_index - 1, kind = fp)
    lat_grid(1 : 2) = lat_s + dlat * real(lat_index - 1, kind = fp)
    z_grid(1 : 2) = dep_min + ddep * real(z_index - 1, kind = fp)
    val_3d(1 : 2, 1 : 2, 1 : 2) = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)

    call linear_interpolation_3d(node_mid%lon, node_mid%lat, node_mid%dep, lon_grid, lat_grid, z_grid, val_3d, v_mid)

    !!calculate velocity gradient at node_mid
    call linear_interpolation_3d(node_mid%lon, node_mid%lat, z_grid(1), lon_grid, lat_grid, z_grid, val_3d, v1)
    call linear_interpolation_3d(node_mid%lon, node_mid%lat, z_grid(2), lon_grid, lat_grid, z_grid, val_3d, v2)
    grad_v_mid(1) = (v2 - v1) / ddep
    call linear_interpolation_3d(node_mid%lon, lat_grid(1), node_mid%dep, lon_grid, lat_grid, z_grid, val_3d, v_mid)
    call linear_interpolation_3d(node_mid%lon, lat_grid(2), node_mid%dep, lon_grid, lat_grid, z_grid, val_3d, v_mid)
    grad_v_mid(2) = -(v2 - v1) / (dlat * node_mid%r)
    call linear_interpolation_3d(lon_grid(1), node_mid%lat, node_mid%dep, lon_grid, lat_grid, z_grid, val_3d, v_mid)
    call linear_interpolation_3d(lon_grid(2), node_mid%lat, node_mid%dep, lon_grid, lat_grid, z_grid, val_3d, v_mid)
    grad_v_mid(3) = (v2 - v1) / (dlon * node_mid%r * sin(node_mid%theta))

    !!calculate normal vector of ray at node_mid
    normal_vector_mid = -1.0_fp / v_mid * (grad_v_mid - dot_product(grad_v_mid, tangentvector_mid) * tangentvector_mid)
    normal_vector_mid = -normal_vector_mid / abs(normal_vector_mid)

    !!calculate distance to move node2
    lon_index = int((node1%lon - lon_w) / dlon) + 1
    lat_index = int((node1%lat - lat_s) / dlat) + 1
    z_index   = int((node1%dep - dep_min) / ddep) + 1
    lon_grid(1 : 2) = lon_w + dlon * real(lon_index - 1, kind = fp)
    lat_grid(1 : 2) = lat_s + dlat * real(lat_index - 1, kind = fp)
    z_grid(1 : 2) = dep_min + ddep * real(z_index - 1, kind = fp)
    val_3d(1 : 2, 1 : 2, 1 : 2) = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
    call linear_interpolation_3d(node1%lon, node1%lat, node1%dep, lon_grid, lat_grid, z_grid, val_3d, v1)
    lon_index = int((node2%lon - lon_w) / dlon) + 1
    lat_index = int((node2%lat - lat_s) / dlat) + 1
    z_index   = int((node2%dep - dep_min) / ddep) + 1
    lon_grid(1 : 2) = lon_w + dlon * real(lon_index - 1, kind = fp)
    lat_grid(1 : 2) = lat_s + dlat * real(lat_index - 1, kind = fp)
    z_grid(1 : 2) = dep_min + ddep * real(z_index - 1, kind = fp)
    val_3d(1 : 2, 1 : 2, 1 : 2) = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
    call linear_interpolation_3d(node2%lon, node2%lat, node2%dep, lon_grid, lat_grid, z_grid, val_3d, v2)
    slowness_mid = (1.0_fp / v1 + 1.0_fp / v2) * 0.5_fp
    const = (slowness_mid * v_mid + 1.0_fp) / (4.0_fp * slowness_mid * dot_product(normal_vector_mid, grad_v_mid))
    dist_move = -const + sqrt(const * const + dist_l ** 2 * 0.5_fp / (slowness_mid * v_mid))

    !!calculate new location of node
    normal_vector_mid = -normal_vector_mid
    node_tmp%r = node_mid%r + dist_move * normal_vector_mid(1)
    node_tmp%theta = node_mid%theta + dist_move * normal_vector_mid(2) / node_mid%r
    node_tmp%phi = node_mid%phi + dist_move * normal_vector_mid(3) / (node_mid%r * sin(node_mid%theta))

    node2%r = node2%r + enhance_factor * (node_tmp%r - node2%r)
    node2%theta = node2%theta + enhance_factor * (node_tmp%theta - node2%theta)
    node2%phi = node2%phi + enhance_factor * (node_tmp%phi - node2%phi)


    node2%dep = r_earth - node2%r
    call latctog(pi * 0.5_fp - node2%theta, node2%lat)
    node2%lon = node2%phi * rad2deg

    return
  end subroutine move_node

  subroutine hypodist(lon0, lat0, dep0, lon1, lat1, dep1, hdist)
    use constants,   only : r_earth
    use greatcircle, only : greatcircle_dist
    implicit none
    real(kind = fp), intent(in)  :: lon0, lat0, dep0, lon1, lat1, dep1
    real(kind = fp), intent(out) :: hdist

    real(kind = fp) :: delta

    call greatcircle_dist(lat0, lon0, lat1, lon1, delta_out = delta)
    hdist = sqrt((r_earth - dep0) ** 2 + (r_earth - dep1) ** 2 &
    &             - 2.0_fp * (r_earth - dep0) * (r_earth - dep1) * cos(delta))

    return
  end subroutine hypodist

  subroutine calc_traveltime_element(node1, node2, velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, ttime_element)
    use linear_interpolation, only : linear_interpolation_3d
    implicit none
    type(gridpoint), intent(in)  :: node1, node2
    real(kind = fp), intent(in)  :: velocity(:, :, :)       !!either velocity of  P- or S-waves
    real(kind = fp), intent(in)  :: lon_w, lat_s, dep_min, dlon, dlat, ddep
    real(kind = fp), intent(out) :: ttime_element
  
    real(kind = fp) :: lon_grid(2), lat_grid(2), z_grid(2), val_3d(2, 2, 2)
    real(kind = fp) :: v1, v2, hdist
    integer :: lon_index, lat_index, z_index
    
    !!calculate mean slowness between two grid points
    lon_index = int((node1%lon - lon_w) / dlon) + 1
    lat_index = int((node1%lat - lat_s) / dlat) + 1
    z_index   = int((node1%dep - dep_min) / ddep) + 1
    lon_grid(1 : 2) = lon_w + dlon * real(lon_index - 1, kind = fp)
    lat_grid(1 : 2) = lat_s + dlat * real(lat_index - 1, kind = fp)
    z_grid(1 : 2) = dep_min + ddep * real(z_index - 1, kind = fp)
    val_3d(1 : 2, 1 : 2, 1 : 2) = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
    call linear_interpolation_3d(node1%lon, node1%lat, node1%dep, lon_grid, lat_grid, z_grid, val_3d, v1)
    lon_index = int((node2%lon - lon_w) / dlon) + 1
    lat_index = int((node2%lat - lat_s) / dlat) + 1
    z_index   = int((node2%dep - dep_min) / ddep) + 1
    lon_grid(1 : 2) = lon_w + dlon * real(lon_index - 1, kind = fp)
    lat_grid(1 : 2) = lat_s + dlat * real(lat_index - 1, kind = fp)
    z_grid(1 : 2) = dep_min + ddep * real(z_index - 1, kind = fp)
    val_3d(1 : 2, 1 : 2, 1 : 2) = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
    call linear_interpolation_3d(node2%lon, node2%lat, node2%dep, lon_grid, lat_grid, z_grid, val_3d, v2)
    slowness_mid = (1.0_fp / v1 + 1.0_fp / v2) * 0.5_fp

    !!calculate distance between two grid points
    call hypodist(node1%lon, node1%lat, node1%dep, node2%lon, node2%lat, node2%dep, hdist)

    !!calculate element of travel time
    ttime_element = hypodist * slowness_mid

    return
  end subroutine calc_traveltime_element


end module raybending

