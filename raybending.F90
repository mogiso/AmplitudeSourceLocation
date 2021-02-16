module raybending
  implicit none
  private
  public :: pseudobending3D

  !!Ray tracing using pseudobending method in a 3-D medium (spherical coordinate)
  !!Reference: Um and Thurber (1987, BSSA); Koketsu and Sekine (1998, GJI)
  !!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2021
  !!License  : MIT License (https://opensource.org/licenses/MIT)

contains

  subroutine pseudobending3D(raypath_lon, raypath_lat, raypath_dep, nraypath, ndiv_raypath, &
  &                          velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, &
  &                          traveltime, ray_az, ray_incangle, &
  &                          qinv, lon_w_qinv, lat_s_qinv, dep_min_qinv, dlon_qinv, dlat_qinv, ddep_qinv, &
  &                          pulsewidth)
    use def_gridpoint
    use nrtype, only : fp
    use constants, only : r_earth, deg2rad, rad2deg, pi
    use greatcircle, only : latgtoc, latctog, greatcircle_dist

    !!array size of raypath_(lon|lat|dep) is (nraypath(initial val) - 1) * 2 ** ndiv_raypath + 1
    real(kind = fp), intent(inout) :: raypath_lon(:), raypath_lat(:), raypath_dep(:)
    integer,         intent(inout) :: nraypath
    integer,         intent(in)    :: ndiv_raypath
    real(kind = fp), intent(out)   :: traveltime                              !!second
    real(kind = fp), intent(in)    :: velocity(:, :, :)                       !!either velocity of P- or S-waves
    real(kind = fp), intent(in)    :: lon_w, lat_s, dep_min, dlon, dlat, ddep
    real(kind = fp), intent(out), optional :: ray_az, ray_incangle          !!radian
    real(kind = fp), intent(in),  optional :: qinv(:, :, :)
    real(kind = fp), intent(in),  optional :: lon_w_qinv, lat_s_qinv, dep_min_qinv, dlon_qinv, dlat_qinv, ddep_qinv
    real(kind = fp), intent(out), optional :: pulsewidth                    !!second (t-star)

    real(kind = fp),  parameter :: traveltime_diff_threshold = 0.01_fp
    real(kind = fp),  parameter :: enhance_factor = 1.5_fp

    type(gridpoint)             :: raynode(size(raypath_lon)), raynode_old(size(raypath_lon))
    real(kind = fp)             :: traveltime_ini, traveltime_double, traveltime_tmp, &
    &                              traveltime_raybend_new, traveltime_raybend_old, pulsewidth_tmp
    real(kind = fp)             :: cos_psi, dist_tmp
    integer                     :: i, j, raynode_index_mid, k

    do i = 1, nraypath
      raynode(i)%lon = raypath_lon(i)
      raynode(i)%lat = raypath_lat(i)
      raynode(i)%dep = raypath_dep(i)
      raynode(i)%r = r_earth - raynode(i)%dep
      raynode(i)%phi = raynode(i)%lon * deg2rad
      call latgtoc(raynode(i)%lat, raynode(i)%theta)
      raynode(i)%theta = pi * 0.5_fp - raynode(i)%theta
    enddo

    !!calculate traveltime of initial raypath
    traveltime_ini = 0.0_fp
    do i = 1, nraypath - 1
      call calc_traveltime_element(raynode(i), raynode(i + 1), velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, traveltime_tmp)
      traveltime_ini = traveltime_ini + traveltime_tmp
    enddo
    !write(0, '(a, e15.7)') "traveltime of initial raypath = ", traveltime_ini
    
    !!do psuedobending
    raypath_divide: do j = 1, ndiv_raypath   !!do loop for dividing raypath
      !write(0, '(a, i0)') "count of ndiv_raypath = ", j

      !!add new raynodes
      !!change array indices of existing nodes
      do i = nraypath, 2, -1
       !write(0, *) j, "move ", i, " to ", 2 * i - 1
       raynode(2 * i - 1) = raynode(i)
      enddo

      !!each new raynode is located at the midpoint of existing nodes
      nraypath = nraypath * 2 - 1
      do i = 2, nraypath - 1, 2
        call hypodist(raynode(i + 1)%lon, raynode(i + 1)%lat, raynode(i + 1)%dep, &
        &             raynode(i - 1)%lon, raynode(i - 1)%lat, raynode(i - 1)%dep, dist_tmp)
        cos_psi = (dist_tmp ** 2 + raynode(i - 1)%r ** 2 - raynode(i + 1)%r ** 2) / (2.0_fp * dist_tmp * raynode(i - 1)%r)
        raynode(i)%r = sqrt(raynode(i - 1)%r ** 2 + dist_tmp ** 2 * 0.25_fp - raynode(i - 1)%r * dist_tmp * cos_psi)
        raynode(i)%theta = (raynode(i + 1)%theta + raynode(i - 1)%theta) * 0.5_fp
        raynode(i)%phi = (raynode(i + 1)%phi + raynode(i - 1)%phi) * 0.5_fp
        raynode(i)%dep = r_earth - raynode(i)%r
        raynode(i)%lon = raynode(i)%phi * rad2deg
        raynode(i)%lat = pi * 0.5_fp - raynode(i)%theta; call latctog(raynode(i)%lat, raynode(i)%lat)
      enddo

      !!calculate travel time
      traveltime_double = 0.0_fp
      do i = 1, nraypath - 1
        call calc_traveltime_element(raynode(i), raynode(i + 1), velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, &
        &                            traveltime_tmp)
        traveltime_double = traveltime_double + traveltime_tmp
      enddo

      if(abs(traveltime_ini - traveltime_double) .lt. traveltime_diff_threshold .and. j .gt. 1) then
        traveltime = traveltime_double
        exit raypath_divide
      endif

      !!perturb raynode
      traveltime_raybend_old = traveltime_double
      raynode_old(1 : nraypath) = raynode(1 : nraypath)
      raybend_loop: do
        raynode_index_mid = nraypath / 2 + 1
        if(raynode_index_mid .gt. 2) then
          do i = 2, raynode_index_mid - 1
            !write(0, '(a, 2(i0, 1x))') "raybend node = ", i, nraypath - i + 1
            call move_node(raynode(i - 1), raynode(i), raynode(i + 1), &
            &              velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, enhance_factor)
            call move_node(raynode(nraypath - i), raynode(nraypath - i + 1), raynode(nraypath - i + 2), &
            &              velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, enhance_factor)
            !do k = 1, nraypath
            !  if(raynode(k)%lon .lt. 135.0 .or. raynode(k)%lon .gt. 138.0 .or. &
            !  &  raynode(k)%lat .lt. 32.0 .or. raynode(k)%lat .gt. 34.5 .or. &
            !  &  raynode(k)%dep .lt. 0.0 .or. raynode(k)%dep .gt. 50.0) then
            !    write(0, '(2(i0, 1x))') k, i, nraypath - i + 1
            !    write(0, '(a)') "i"
            !    write(0, '(3(e15.7, 1x))') raynode_old(i - 1)%lon, raynode_old(i - 1)%lat, raynode_old(i - 1)%dep
            !    write(0, '(3(e15.7, 1x))') raynode_old(i)%lon, raynode_old(i)%lat, raynode_old(i)%dep
            !    write(0, '(3(e15.7, 1x))') raynode_old(i + 1)%lon, raynode_old(i + 1)%lat, raynode_old(i + 1)%dep
            !    write(0, '(3(e15.7, 1x))') raynode(i - 1)%lon, raynode(i - 1)%lat, raynode(i - 1)%dep
            !    write(0, '(3(e15.7, 1x))') raynode(i)%lon, raynode(i)%lat, raynode(i)%dep
            !    write(0, '(3(e15.7, 1x))') raynode(i + 1)%lon, raynode(i + 1)%lat, raynode(i + 1)%dep
            !    write(0, '(a)') "nraypath - i + 1"
            !    write(0, '(3(e15.7, 1x))') raynode_old(nraypath - i)%lon, raynode_old(nraypath - i)%lat, raynode_old(nraypath - i)%dep
            !    write(0, '(3(e15.7, 1x))') raynode_old(nraypath - i + 1)%lon, raynode_old(nraypath - i + 1)%lat, raynode_old(nraypath - i + 1)%dep
            !    write(0, '(3(e15.7, 1x))') raynode_old(nraypath - i + 2)%lon, raynode_old(nraypath - i + 2)%lat, raynode_old(nraypath - i + 2)%dep
            !    write(0, '(3(e15.7, 1x))') raynode(nraypath - i)%lon, raynode(nraypath - i)%lat, raynode(nraypath - i)%dep
            !    write(0, '(3(e15.7, 1x))') raynode(nraypath - i + 1)%lon, raynode(nraypath - i + 1)%lat, raynode(nraypath - i + 1)%dep
            !    write(0, '(3(e15.7, 1x))') raynode(nraypath - i + 2)%lon, raynode(nraypath - i + 2)%lat, raynode(nraypath - i + 2)%dep
            !    write(0, *) j, traveltime_ini, traveltime_double
            !    write(0, *) traveltime_raybend_new, traveltime_raybend_old
            !  endif
            !enddo
          enddo
        endif
        if(mod(nraypath, 2) .ne. 0) then
          call move_node(raynode(raynode_index_mid - 1), raynode(raynode_index_mid), raynode(raynode_index_mid + 1), &
          &              velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, enhance_factor)
        endif

            

        traveltime_raybend_new = 0.0_fp
        do i = 1, nraypath - 1
          call calc_traveltime_element(raynode(i), raynode(i + 1), velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, &
          &                            traveltime_tmp)
          traveltime_raybend_new = traveltime_raybend_new + traveltime_tmp
        enddo

        if(abs(traveltime_raybend_new - traveltime_raybend_old) .lt. traveltime_diff_threshold) then
          if(traveltime_raybend_new .lt. traveltime_raybend_old) then
            traveltime_ini = traveltime_raybend_new
            exit raybend_loop
          else
            raynode(1 : nraypath) = raynode_old(1 : nraypath)
            traveltime_ini = traveltime_raybend_old
          endif
        endif
        traveltime_raybend_old = traveltime_raybend_new
      enddo raybend_loop
    enddo raypath_divide
    !write(0, '(a, e15.7)') "traveltime after bending = ", traveltime

    do i = 1, nraypath
      raypath_lon(i) = raynode(i)%lon
      raypath_lat(i) = raynode(i)%lat
      raypath_dep(i) = raynode(i)%dep
    enddo

    if(present(ray_az)) then
      call greatcircle_dist(raynode(1)%lat, raynode(1)%lon, raynode(2)%lat, raynode(2)%lon, azimuth=ray_az)
    endif
    if(present(ray_incangle)) then
      call greatcircle_dist(raynode(1)%lat, raynode(1)%lon, raynode(2)%lat, raynode(2)%lon, delta_out=ray_incangle)
      ray_incangle = raynode(1)%r * ray_incangle
      ray_incangle = atan2(ray_incangle, raynode(2)%dep - raynode(1)%dep)
    endif

    if(present(pulsewidth)) then
      pulsewidth = 0.0_fp
      do i = 1, nraypath - 1
        call calc_traveltime_element(raynode(i), raynode(i + 1), velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, &
        &                            traveltime_tmp, &
        &                            qinv, lon_w_qinv, lat_s_qinv, dep_min_qinv, dlon_qinv, dlat_qinv, ddep_qinv, &
        &                            pulsewidth_tmp)
        pulsewidth = pulsewidth + pulsewidth_tmp
      enddo
    endif

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
    integer :: lon_index, lat_index, z_index, i, nlon, nlat, ndep
    type(gridpoint) :: node_mid, node_tmp


    nlon = ubound(velocity, 1)
    nlat = ubound(velocity, 2)
    ndep = ubound(velocity, 3)

    !!calculate distance between node1 and node3, then divide by 2
    call hypodist(node1%lon, node1%lat, node1%dep, node3%lon, node3%lat, node3%dep, dist_l)
    dist_l = dist_l * 0.5_fp
    !print *, "dist node1-mid or mid-3", dist_l

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
    z_grid(1)   = dep_min + ddep * real(z_index - 1, kind = fp); z_grid(2)   = z_grid(1)   + ddep
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

    !!calculate normal vector of ray at node_mid
    if(grad_vmid(1) ** 2 + grad_vmid(2) ** 2 + grad_vmid(3) ** 2 .eq. 0.0_fp) then
      normalvector_mid(1 : 3) = 0.0_fp
    else
      normalvector_mid(1 : 3) = -1.0_fp / vmid &
      &                        * (grad_vmid(1 : 3) - dot_product(grad_vmid, tangentvector_mid) * tangentvector_mid(1 : 3))
      vector_len = sqrt(normalvector_mid(1) ** 2 + normalvector_mid(2) ** 2 + normalvector_mid(3) ** 2)
      normalvector_mid(1 : 3) = -normalvector_mid(1 : 3) / vector_len
    endif

    !print *, "tangentvector", (tangentvector_mid(i), i = 1, 3)
    !print *, "length of tangentvector", sqrt(tangentvector_mid(1) ** 2 + tangentvector_mid(2) ** 2 + tangentvector_mid(3) ** 2)
    !print *, "normalvector", (normalvector_mid(i), i = 1, 3)
    !print *, "grad_vmid", (grad_vmid(i), i = 1, 3)
    !check
    !print *, dot_product(tangentvector_mid, normalvector_mid)

    !!calculate distance to move node2
    lon_index = int((node1%lon - lon_w) / dlon) + 1
    lat_index = int((node1%lat - lat_s) / dlat) + 1
    z_index   = int((node1%dep - dep_min) / ddep) + 1
    lon_grid(1) = lon_w + dlon * real(lon_index - 1, kind = fp); lon_grid(2) = lon_grid(1) + dlon
    lat_grid(1) = lat_s + dlat * real(lat_index - 1, kind = fp); lat_grid(2) = lat_grid(1) + dlat
    z_grid(1)   = dep_min + ddep * real(z_index - 1, kind = fp); z_grid(2)   = z_grid(1)   + ddep
    val_3d(1 : 2, 1 : 2, 1 : 2) = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
    call linear_interpolation_3d(node1%lon, node1%lat, node1%dep, lon_grid, lat_grid, z_grid, val_3d, v1)

    lon_index = int((node3%lon - lon_w) / dlon) + 1
    lat_index = int((node3%lat - lat_s) / dlat) + 1
    z_index   = int((node3%dep - dep_min) / ddep) + 1
    lon_grid(1) = lon_w + dlon * real(lon_index - 1, kind = fp); lon_grid(2) = lon_grid(1) + dlon
    lat_grid(1) = lat_s + dlat * real(lat_index - 1, kind = fp); lat_grid(2) = lat_grid(1) + dlat
    z_grid(1)   = dep_min + ddep * real(z_index - 1, kind = fp); z_grid(2)   = z_grid(1)   + ddep
    val_3d(1 : 2, 1 : 2, 1 : 2) = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
    call linear_interpolation_3d(node3%lon, node3%lat, node3%dep, lon_grid, lat_grid, z_grid, val_3d, v2)

    if(normalvector_mid(1) ** 2 + normalvector_mid(2) ** 2 + normalvector_mid(3) ** 2 .eq. 0.0_fp) then
      dist_move = 0.0_fp
    else
      slowness_mid = (1.0_fp / v1 + 1.0_fp / v2) * 0.5_fp
      const = (slowness_mid * vmid + 1.0_fp) / (4.0_fp * slowness_mid * dot_product(normalvector_mid, grad_vmid))
      dist_move = -const + sqrt(const * const + dist_l ** 2 * 0.5_fp / (slowness_mid * vmid))
    endif

    !!calculate new location of node
    node_tmp%r = node_mid%r + dist_move * normalvector_mid(1)
    node_tmp%theta = node_mid%theta + dist_move * normalvector_mid(2) / node_mid%r
    node_tmp%phi = node_mid%phi + dist_move * normalvector_mid(3) / (node_mid%r * sin(node_mid%theta))

    node2%r = node2%r + enhance_factor * (node_tmp%r - node2%r)
    node2%theta = node2%theta + enhance_factor * (node_tmp%theta - node2%theta)
    node2%phi = node2%phi + enhance_factor * (node_tmp%phi - node2%phi)
    node2%dep = r_earth - node2%r
    call latctog(pi * 0.5_fp - node2%theta, node2%lat)
    node2%lon = node2%phi * rad2deg

    !!ad-hoc modification when the ray intersects boundary of velocity array
    if(node2%lon .lt. lon_w .or. node2%lon .gt. lon_w + dlon * real(nlon - 1, kind = fp) .or. &
    &  node2%lat .lt. lat_s .or. node2%lat .gt. lat_s + dlat * real(nlat - 1, kind = fp) .or. &
    &  node2%dep .lt. dep_min .or. node2%dep .gt. dep_min + ddep * real(ndep - 1, kind = fp)) then

      !write(0, '(a)') "Ad-hoc modification"
      node2%r = node_mid%r + dist_move * 0.1_fp * normalvector_mid(1)
      node2%theta = node_mid%theta + dist_move * 0.1_fp * normalvector_mid(2) / node_mid%r
      node2%phi = node_mid%phi + dist_move * 0.1_fp * normalvector_mid(3) / (node_mid%r * sin(node_mid%theta))
      node2%dep = r_earth - node2%r
      call latctog(pi * 0.5_fp - node2%theta, node2%lat)
      node2%lon = node2%phi * rad2deg
    endif

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

  subroutine calc_traveltime_element(node1, node2, velocity, lon_w, lat_s, dep_min, dlon, dlat, ddep, ttime_element, &
  &                                  qinv, lon_w_qinv, lat_s_qinv, dep_min_qinv, dlon_qinv, dlat_qinv, ddep_qinv, &
  &                                  pulsewidth_element)
    use nrtype,               only : fp
    use linear_interpolation, only : linear_interpolation_3d, block_interpolation_3d
    use def_gridpoint
    implicit none
    type(gridpoint), intent(in)  :: node1, node2
    real(kind = fp), intent(in)  :: velocity(:, :, :)       !!either velocity of  P- or S-waves
    real(kind = fp), intent(in)  :: lon_w, lat_s, dep_min, dlon, dlat, ddep
    real(kind = fp), intent(out) :: ttime_element
    real(kind = fp), intent(in),  optional :: qinv(:, :, :)       !!either Qinv of  P- or S-waves
    real(kind = fp), intent(in),  optional :: lon_w_qinv, lat_s_qinv, dep_min_qinv, dlon_qinv, dlat_qinv, ddep_qinv
    real(kind = fp), intent(out), optional :: pulsewidth_element
  
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
    if(lon_index .ge. ubound(velocity, 1)) then
      print *, lon_index, node2%lon, node2%lat, node2%dep
    endif
    val_3d(1 : 2, 1 : 2, 1 : 2) = velocity(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
    call linear_interpolation_3d(node2%lon, node2%lat, node2%dep, lon_grid, lat_grid, z_grid, val_3d, v2)
    slowness_mid = (1.0_fp / v1 + 1.0_fp / v2) * 0.5_fp

    !!calculate distance between two grid points
    call hypodist(node1%lon, node1%lat, node1%dep, node2%lon, node2%lat, node2%dep, hdist)

    !!calculate element of travel time
    ttime_element = hdist * slowness_mid

    if(present(pulsewidth_element)) then
      lon_index = int((node1%lon - lon_w_qinv) / dlon_qinv) + 1
      lat_index = int((node1%lat - lat_s_qinv) / dlat_qinv) + 1
      z_index   = int((node1%dep - dep_min_qinv) / ddep_qinv) + 1
      lon_grid(1) = lon_w_qinv + dlon_qinv * real(lon_index - 1, kind = fp); lon_grid(2) = lon_grid(1) + 1
      lat_grid(1) = lat_s_qinv + dlat_qinv * real(lat_index - 1, kind = fp); lat_grid(2) = lat_grid(1) + 1
      z_grid(1) = dep_min_qinv + ddep_qinv * real(z_index - 1, kind = fp); z_grid(2) = z_grid(1) + 1
      val_3d(1 : 2, 1 : 2, 1 : 2) = qinv(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
      call block_interpolation_3d(node1%lon, node1%lat, node1%dep, lon_grid, lat_grid, z_grid, val_3d, v1)
      lon_index = int((node2%lon - lon_w_qinv) / dlon_qinv) + 1
      lat_index = int((node2%lat - lat_s_qinv) / dlat_qinv) + 1
      z_index   = int((node2%dep - dep_min_qinv) / ddep_qinv) + 1
      lon_grid(1) = lon_w_qinv + dlon_qinv * real(lon_index - 1, kind = fp); lon_grid(2) = lon_grid(1) + 1
      lat_grid(1) = lat_s_qinv + dlat_qinv * real(lat_index - 1, kind = fp); lat_grid(2) = lat_grid(1) + 1
      z_grid(1) = dep_min_qinv + ddep_qinv * real(z_index - 1, kind = fp); z_grid(2) = z_grid(1) + 1
      val_3d(1 : 2, 1 : 2, 1 : 2) = qinv(lon_index : lon_index + 1, lat_index : lat_index + 1, z_index : z_index + 1)
      call block_interpolation_3d(node2%lon, node2%lat, node2%dep, lon_grid, lat_grid, z_grid, val_3d, v2)

      pulsewidth_element = ttime_element * (v1 + v2) * 0.5_fp
    endif

    return
  end subroutine calc_traveltime_element


end module raybending

