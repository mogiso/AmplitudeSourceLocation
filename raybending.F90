module raybending
  implicit none
  private
  public :: pseudobending3D

contains

  subroutine pseudobending3D(lon_source, lat_source, dep_source, &
  &                          lon_receiver, lat_receiver, dep_receiver, &
  &                          velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, ndiv_raypath)
    use def_gridpoint
    use nrtype, only : fp
    use constants, only : r_earth, deg2rad, rad2deg, pi
    use greatcircle, only : latgtoc, latctog

    real(kind = fp), intent(in) :: lon_source, lat_source, dep_source, lon_receiver, lat_receiver, dep_receiver
    real(kind = fp), intent(in) :: velocity(:, :, :)
    real(kind = fp), intent(in) :: lon_w, lat_s, dep_min, dlon, dlat, ddep
    integer,         intent(in) :: ndiv_raypath

    real(kind = fp),  parameter :: traveltime_diff_threshold = 0.01_fp
    real(kind = fp),  parameter :: enhance_factor = 1.3_fp

    type(gridpoint)             :: raynode(2 + 2 ** ndiv_raypath - 1)
    real(kind = fp)             :: traveltime_ini, traveltime_new, traveltime_tmp, traveltime_raybend_new, traveltime_raybend_old
    real(kind = fp)             :: cos_psi, dist_tmp
    integer                     :: i, j, nlon, nlat, ndep, nraynode, raynode_index_mid

    nlon = ubound(velocity, 1)
    nlat = ubound(velocity, 2)
    ndep = ubound(velocity, 3)

    raynode(1)%lon = lon_source
    raynode(1)%lat = lat_source
    raynode(1)%dep = dep_source
    raynode(1)%r = r_earth - raynode(1)%dep
    raynode(1)%phi = raynode(1)%lon * deg2rad
    call latgtoc(raynode(1)%lat, raynode(1)%theta); raynode(1)%theta = pi * 0.5_fp - raynode(1)%theta

    raynode(2)%lon = lon_receiver
    raynode(2)%lat = lat_receiver
    raynode(2)%dep = dep_receiver
    raynode(2)%r = r_earth - raynode(2)%dep
    raynode(2)%phi = raynode(2)%lon * deg2rad
    call latgtoc(raynode(2)%lat, raynode(2)%theta); raynode(2)%theta = pi * 0.5_fp - raynode(2)%theta

    !!state initial situation
    nraynode = 2 
    call calc_traveltime_element(raynode(1), raynode(2), velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, traveltime_ini)
    print *, "traveltime_ini = ", traveltime_ini
    
    raypath_divide: do j = 1, ndiv_raypath   !!do loop for dividing raypath
      !!change array indices of exsisting nodes
      do i = nraynode, 2, -1
       raynode(i + 2 ** (j - 1) - (nraynode - i)) = raynode(i)
      enddo

      !!add new raynodes
      !!each new raynode locates at the midpoint of existing nodes
      nraynode = nraynode + 2 ** (j - 1)
      do i = 2, nraynode - 1, 2
        !raynode(i)%r = (raynode(i + 1)%r + raynode(i - 1)%r) * 0.5_fp
        call hypodist(raynode(i + 1)%lat, raynode(i + 1)%lon, raynode(i + 1)%dep, &
        &             raynode(i - 1)%lat, raynode(i - 1)%lon, raynode(i - 1)%dep, dist_tmp)
        cos_psi = (dist_tmp ** 2 + raynode(i - 1)%r ** 2 - raynode(i + 1)%r ** 2) / (2.0_fp * dist_tmp * raynode(i - 1)%r)
        raynode(i)%r = sqrt(raynode(i - 1)%r ** 2 + dist_tmp ** 2 * 0.25_fp - raynode(i - 1)%r * dist_tmp * cos_psi)
        raynode(i)%theta = (raynode(i + 1)%theta + raynode(i - 1)%theta) * 0.5_fp
        raynode(i)%phi = (raynode(i + 1)%phi + raynode(i - 1)%phi) * 0.5_fp
        raynode(i)%dep = r_earth - raynode(i)%r
        raynode(i)%lon = raynode(i)%phi * rad2deg
        raynode(i)%lat = pi * 0.5_fp - raynode(i)%theta; call latctog(raynode(i)%lat, raynode(i)%lat)
      enddo
      do i = 1, nraynode
         print '(a, i0, a, 3(e15.7, 1x))', "raynode(", i, ") lon, lat, dep = ", raynode(i)%lon, raynode(i)%lat, raynode(i)%dep
      enddo

      print *, "ndiv_raypath_count, j = ", j, "nraynode (aft. added) = ", nraynode

      !!calculate travel time
      traveltime_new = 0.0_fp
      do i = 1, nraynode - 1
        call calc_traveltime_element(raynode(i), raynode(i + 1), velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, &
        &                            traveltime_tmp)
        traveltime_new = traveltime_new + traveltime_tmp
      enddo
      print *, "traveltime_new = ", traveltime_new
      traveltime_ini = traveltime_new
 
      traveltime_raybend_old = traveltime_ini
      raybend_loop: do
        raynode_index_mid = nraynode / 2 + 1
        if(raynode_index_mid .gt. 2) then
          do i = 2, raynode_index_mid - 1
            print *, "raybend node = ", i, nraynode - i + 1
            call move_node(raynode(i - 1), raynode(i), raynode(i + 1), &
            &              velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, enhance_factor)
            call move_node(raynode(nraynode - i), raynode(nraynode - i + 1), raynode(nraynode - i + 2), &
            &              velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, enhance_factor)
          enddo
        endif
        !print *, raynode(raynode_index_mid)%lon, raynode(raynode_index_mid)%lat, raynode(raynode_index_mid)%dep
        print *, raynode(raynode_index_mid)%r, raynode(raynode_index_mid)%theta * rad2deg, &
        &        raynode(raynode_index_mid)%phi * rad2deg, "old"
        call move_node(raynode(raynode_index_mid - 1), raynode(raynode_index_mid), raynode(raynode_index_mid + 1), &
        &              velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, enhance_factor)
        !print *, raynode(raynode_index_mid)%lon, raynode(raynode_index_mid)%lat, raynode(raynode_index_mid)%dep
        print *, raynode(raynode_index_mid)%r, raynode(raynode_index_mid)%theta * rad2deg, &
        &        raynode(raynode_index_mid)%phi * rad2deg, "new"

        traveltime_raybend_new = 0.0_fp
        do i = 1, nraynode - 1
          call calc_traveltime_element(raynode(i), raynode(i + 1), velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, &
          &                            traveltime_tmp)
          traveltime_raybend_new = traveltime_raybend_new + traveltime_tmp
        enddo
        print *, "traveltime_raybend_new = ", traveltime_raybend_new

        if(abs(traveltime_raybend_old - traveltime_raybend_new) .lt. traveltime_diff_threshold) then
          exit raybend_loop
        endif
        traveltime_raybend_old = traveltime_raybend_new
      enddo raybend_loop
    enddo raypath_divide

    return
  end subroutine pseudobending3D

  subroutine move_node(node1, node2, node3, velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, enhance_factor)
    !!three-point perturbation scheme of raypath based on Koketsu and Sekine (1998) eq. 19
    !!node1 and node3 are fixed, move node2
    use def_gridpoint
    use nrtype,               only : fp
    use constants,            only : r_earth, pi, rad2deg, deg2rad
    use greatcircle,          only : latctog
    use linear_interpolation, only : linear_interpolation_3d
    implicit none
    type(gridpoint), intent(in)    :: node1, node3
    type(gridpoint), intent(inout) :: node2
    real(kind = fp), intent(in)    :: velocity(:, :, :)       !!either velocity of P- or S-waves
    real(kind = fp), intent(in)    :: lon_w, lat_s, dep_min, dlon, dlat, ddep, enhance_factor

    real(kind = fp) :: dist_l, slowness_mid, const, dist_move, v1, v2, vmid, cos_psi, vector_len
    real(kind = fp) :: tangentvector_mid(3), normalvector_mid(3), grad_vmid(3), &
    &                  lon_grid(2), lat_grid(2), z_grid(2), val_3d(2, 2, 2)
    integer :: lon_index, lat_index, z_index, i
    type(gridpoint) :: node_mid, node_tmp

    !!calculate distance between node1 and node3, then divide by 2
    call hypodist(node1%lat, node1%lon, node1%dep, node3%lat, node3%lon, node3%dep, dist_l)
    dist_l = dist_l * 0.5_fp

    !!calculate location of midpoint between node1 and node3
    node_mid%theta = 0.5_fp * (node1%theta + node3%theta)
    node_mid%phi = 0.5_fp * (node1%phi + node3%phi)
    !!psi is the angle between dist_l and node1%r
    cos_psi = (4.0_fp * dist_l ** 2 + node1%r ** 2 - node3%r ** 2) / (4.0_fp * dist_l * node1%r)
    node_mid%r = sqrt(node1%r ** 2 + dist_l ** 2 - 2.0_fp * node1%r * dist_l * cos_psi)
    node_mid%dep = r_earth - node_mid%r
    node_mid%lon = node_mid%phi * rad2deg
    call latctog(pi * 0.5_fp - node_mid%theta, node_mid%lat)

    !!calculate tangent vector at node_mid
    tangentvector_mid(1) = 0.5_fp / dist_l * (node3%r - node1%r)
    tangentvector_mid(2) = 0.5_fp / dist_l * node_mid%r * (node3%theta - node1%theta)
    tangentvector_mid(3) = 0.5_fp / dist_l * node_mid%r * sin(node_mid%theta) * (node3%phi - node1%phi)

    !!calculate velocity at node_mid
    lon_index = int((node_mid%lon - lon_w) / dlon) + 1
    lat_index = int((node_mid%lat - lat_s) / dlat) + 1
    z_index   = int((node_mid%dep - dep_min) / ddep) + 1
    lon_grid(1) = lon_w + dlon * real(lon_index - 1, kind = fp); lon_grid(2) = lon_grid(1) + dlon
    lat_grid(1) = lat_s + dlat * real(lat_index - 1, kind = fp); lat_grid(2) = lat_grid(1) + dlat
    z_grid(1) = dep_min + ddep * real(z_index - 1, kind = fp); z_grid(2) = z_grid(1) + ddep
    val_3d(1 : 2, 1 : 2, 1 : 2) = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)

    call linear_interpolation_3d(node_mid%lon, node_mid%lat, node_mid%dep, lon_grid, lat_grid, z_grid, val_3d, vmid)

    !!calculate velocity gradient at node_mid
    grad_vmid(1) = (velocity(lon_index, lat_index, z_index + 1) - velocity(lon_index, lat_index, z_index)) &
    &            + (velocity(lon_index + 1, lat_index, z_index + 1) - velocity(lon_index + 1, lat_index, z_index)) &
    &            + (velocity(lon_index, lat_index + 1, z_index + 1) - velocity(lon_index, lat_index + 1, z_index)) &
    &            + (velocity(lon_index + 1, lat_index + 1, z_index + 1) - velocity(lon_index + 1, lat_index + 1, z_index))
    grad_vmid(1) = grad_vmid(1) * (-0.25_fp) / ddep

    grad_vmid(2) = (velocity(lon_index, lat_index + 1, z_index) - velocity(lon_index, lat_index, z_index)) &
    &            + (velocity(lon_index + 1, lat_index + 1, z_index) - velocity(lon_index + 1, lat_index, z_index)) &
    &            + (velocity(lon_index, lat_index + 1, z_index + 1) - velocity(lon_index, lat_index, z_index + 1)) &
    &            + (velocity(lon_index + 1, lat_index + 1, z_index + 1) - velocity(lon_index + 1, lat_index, z_index + 1))
    grad_vmid(2) = grad_vmid(2) * (-0.25_fp) / (dlat * node_mid%r)

    grad_vmid(3) = (velocity(lon_index + 1, lat_index, z_index) - velocity(lon_index, lat_index, z_index)) &
    &            + (velocity(lon_index + 1, lat_index + 1, z_index) - velocity(lon_index, lat_index + 1, z_index)) &
    &            + (velocity(lon_index + 1, lat_index, z_index + 1) - velocity(lon_index, lat_index, z_index + 1)) &
    &            + (velocity(lon_index + 1, lat_index + 1, z_index + 1) - velocity(lon_index, lat_index + 1, z_index + 1))
    grad_vmid(3) = grad_vmid(3) * 0.25_fp / (dlon * node_mid%r * sin(node_mid%theta))

    !print *, "dot_product(grad_vmid, tangentvector) = ", dot_product(grad_vmid, tangentvector_mid)

    !!calculate normal vector of ray at node_mid
    normalvector_mid(1 : 3) = -1.0_fp / vmid &
    &                        * (grad_vmid(1 : 3) - dot_product(grad_vmid, tangentvector_mid) * tangentvector_mid(1 : 3))
    vector_len = sqrt(normalvector_mid(1) ** 2 + normalvector_mid(2) ** 2 + normalvector_mid(3) ** 2)
    normalvector_mid(1 : 3) = normalvector_mid(1 : 3) / vector_len

    print *, "tangentvector", (tangentvector_mid(i), i = 1, 3)
    print *, "normalvector", (normalvector_mid(i), i = 1, 3)
    print *, "grad_vmid", (grad_vmid(i), i = 1, 3)

    !check
    print *, dot_product(tangentvector_mid, normalvector_mid)

    !!calculate distance to move node2
    lon_index = int((node1%lon - lon_w) / dlon) + 1
    lat_index = int((node1%lat - lat_s) / dlat) + 1
    z_index   = int((node1%dep - dep_min) / ddep) + 1
    lon_grid(1) = lon_w + dlon * real(lon_index - 1, kind = fp); lon_grid(2) = lon_grid(1) + dlon
    lat_grid(1) = lat_s + dlat * real(lat_index - 1, kind = fp); lat_grid(2) = lat_grid(1) + dlat
    z_grid(1) = dep_min + ddep * real(z_index - 1, kind = fp); z_grid(2) = z_grid(1) + ddep
    val_3d(1 : 2, 1 : 2, 1 : 2) = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
    call linear_interpolation_3d(node1%lon, node1%lat, node1%dep, lon_grid, lat_grid, z_grid, val_3d, v1)

    lon_index = int((node2%lon - lon_w) / dlon) + 1
    lat_index = int((node2%lat - lat_s) / dlat) + 1
    z_index   = int((node2%dep - dep_min) / ddep) + 1
    lon_grid(1) = lon_w + dlon * real(lon_index - 1, kind = fp); lon_grid(2) = lon_grid(1) + dlon
    lat_grid(1) = lat_s + dlat * real(lat_index - 1, kind = fp); lat_grid(2) = lat_grid(1) + dlat
    z_grid(1) = dep_min + ddep * real(z_index - 1, kind = fp); z_grid(2) = z_grid(1) + ddep
    val_3d(1 : 2, 1 : 2, 1 : 2) = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
    call linear_interpolation_3d(node2%lon, node2%lat, node2%dep, lon_grid, lat_grid, z_grid, val_3d, v2)

    slowness_mid = (1.0_fp / v1 + 1.0_fp / v2) * 0.5_fp
    const = (slowness_mid * vmid + 1.0_fp) / (4.0_fp * slowness_mid * dot_product(normalvector_mid, grad_vmid))
    dist_move = -const + sqrt(const * const + dist_l ** 2 * 0.5_fp / (slowness_mid * vmid))

    print *, "dist_move = ", dist_move
  
    !!calculate new location of node
    normalvector_mid = -normalvector_mid
    node_tmp%r = node_mid%r + dist_move * normalvector_mid(1)
    node_tmp%theta = node_mid%theta + dist_move * normalvector_mid(2) / node_mid%r
    node_tmp%phi = node_mid%phi + dist_move * normalvector_mid(3) / (node_mid%r * sin(node_mid%theta))

    node2%r = node2%r + enhance_factor * (node_tmp%r - node2%r)
    node2%theta = node2%theta + enhance_factor * (node_tmp%theta - node2%theta)
    node2%phi = node2%phi + enhance_factor * (node_tmp%phi - node2%phi)


    node2%dep = r_earth - node2%r
    call latctog(pi * 0.5_fp - node2%theta, node2%lat)
    node2%lon = node2%phi * rad2deg

    return
  end subroutine move_node

  subroutine hypodist(lon0, lat0, dep0, lon1, lat1, dep1, hdist)
    use nrtype,      only : fp
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
    use nrtype,               only : fp
    use linear_interpolation, only : linear_interpolation_3d
    use def_gridpoint
    implicit none
    type(gridpoint), intent(in)  :: node1, node2
    real(kind = fp), intent(in)  :: velocity(:, :, :)       !!either velocity of  P- or S-waves
    real(kind = fp), intent(in)  :: lon_w, lat_s, dep_min, dlon, dlat, ddep
    real(kind = fp), intent(out) :: ttime_element
  
    real(kind = fp) :: lon_grid(2), lat_grid(2), z_grid(2), val_3d(2, 2, 2)
    real(kind = fp) :: v1, v2, hdist, slowness_mid
    integer :: lon_index, lat_index, z_index
    
    !!calculate mean slowness between two grid points
    lon_index = int((node1%lon - lon_w) / dlon) + 1
    lat_index = int((node1%lat - lat_s) / dlat) + 1
    z_index   = int((node1%dep - dep_min) / ddep) + 1
    lon_grid(1) = lon_w + dlon * real(lon_index - 1, kind = fp); lon_grid(2) = lon_grid(1) + 1
    lat_grid(1) = lat_s + dlat * real(lat_index - 1, kind = fp); lat_grid(2) = lat_grid(1) + 1
    z_grid(1) = dep_min + ddep * real(z_index - 1, kind = fp); z_grid(2) = z_grid(1) + 1
    val_3d(1 : 2, 1 : 2, 1 : 2) = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
    call linear_interpolation_3d(node1%lon, node1%lat, node1%dep, lon_grid, lat_grid, z_grid, val_3d, v1)
    lon_index = int((node2%lon - lon_w) / dlon) + 1
    lat_index = int((node2%lat - lat_s) / dlat) + 1
    z_index   = int((node2%dep - dep_min) / ddep) + 1
    lon_grid(1) = lon_w + dlon * real(lon_index - 1, kind = fp); lon_grid(2) = lon_grid(1) + 1
    lat_grid(1) = lat_s + dlat * real(lat_index - 1, kind = fp); lat_grid(2) = lat_grid(1) + 1
    z_grid(1) = dep_min + ddep * real(z_index - 1, kind = fp); z_grid(2) = z_grid(1) + 1
    val_3d(1 : 2, 1 : 2, 1 : 2) = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
    call linear_interpolation_3d(node2%lon, node2%lat, node2%dep, lon_grid, lat_grid, z_grid, val_3d, v2)
    slowness_mid = (1.0_fp / v1 + 1.0_fp / v2) * 0.5_fp

    !!calculate distance between two grid points
    call hypodist(node1%lon, node1%lat, node1%dep, node2%lon, node2%lat, node2%dep, hdist)

    !!calculate element of travel time
    ttime_element = hdist * slowness_mid

    return
  end subroutine calc_traveltime_element


end module raybending

