module grdfile_io
  !!module file for I/O of netcdf-format 2d grdfile
  !!Author: Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2020
  !!License  : MIT License https://opensource.org/licenses/MIT

  private
  public :: read_grdfile_2d, write_grdfile_2d, write_grdfile_fp_2d

  contains

  subroutine read_grdfile_2d(netcdf_file, xval, yval, zval)
    use netcdf
    use nrtype, only : fp
    implicit none
    character(len = *),           intent(in)  :: netcdf_file
    real(kind = fp), allocatable, intent(out) :: xval(:), yval(:)
    real(kind = fp), allocatable, intent(out) :: zval(:, :)

    integer :: ncstatus, ncid, varid_x, varid_y, varid_z, varlen_x, varlen_y

    !!open netcdf-format grd file
    write(0, '(2a)') "reading netcdf-format grdfile ", trim(netcdf_file)
    ncstatus = nf90_open(netcdf_file, nf90_nowrite, ncid)
    
    call get_varinfo(ncid, 'x', varid = varid_x, varlen = varlen_x)
    call get_varinfo(ncid, 'y', varid = varid_y, varlen = varlen_y)
    call get_varinfo(ncid, 'z', varid = varid_z)

    allocate(xval(varlen_x), yval(varlen_y), zval(varlen_x, varlen_y))

    !!read variables
    ncstatus = nf90_get_var(ncid, varid_x, xval)
    ncstatus = nf90_get_var(ncid, varid_y, yval)
    ncstatus = nf90_get_var(ncid, varid_z, zval)

    !!close file
    ncstatus = nf90_close(ncid)

    return
  end subroutine read_grdfile_2d

  subroutine write_grdfile_2d(xmin, ymin, dx, dy, nx, ny, zval, netcdf_file, nanval)
    use netcdf
    use nrtype, only : sp, dp
    implicit none

    real(kind = dp),    intent(in)    :: xmin, ymin, dx, dy
    integer,            intent(in)    :: nx, ny
    real(kind = sp),    intent(in)    :: zval(1 : nx, 1 : ny)
    character(len = *), intent(inout) :: netcdf_file
    real(kind = sp),    intent(in), optional :: nanval

    integer :: i, j, ncstatus, ncid, dimid_x, dimid_y, varid_x, varid_y, varid_z
    real(kind = dp), allocatable :: tmp_array(:)
    real(kind = sp), allocatable :: tmp_array2d(:, :)

    !!open file
    ncstatus = nf90_create(trim(netcdf_file), NF90_NETCDF4, ncid)

    !!dimensions
    ncstatus = nf90_def_dim(ncid, 'x', nx, dimid_x)
    ncstatus = nf90_def_dim(ncid, 'y', ny, dimid_y)

    !!variables
    ncstatus = nf90_def_var(ncid, 'x', NF90_DOUBLE, [dimid_x], varid_x)
    ncstatus = nf90_def_var(ncid, 'y', NF90_DOUBLE, [dimid_y], varid_y)
    ncstatus = nf90_def_var(ncid, 'z', NF90_FLOAT,  [dimid_x, dimid_y], varid_z)

    allocate(tmp_array(1 : nx))
    do i = 1, nx
      tmp_array(i) = xmin + dx * real(i - 1, kind = dp)
    enddo
    ncstatus = nf90_put_var(ncid, varid_x, tmp_array)
    deallocate(tmp_array)

    allocate(tmp_array(1 : ny))
    do i = 1, ny
      tmp_array(i) = ymin + dy * real(i - 1, kind = dp)
    enddo
    ncstatus = nf90_put_var(ncid, varid_y, tmp_array)
    deallocate(tmp_array)
    
    allocate(tmp_array2d(1 : nx, 1 : ny))
    do j = 1, ny
      do i = 1, nx
        tmp_array2d(i, j) = zval(i, j)
        if(present(nanval) .and. zval(i, j) .eq. nanval) tmp_array2d(i, j) = 0.0_sp / 0.0_sp
      enddo
    enddo
    ncstatus = nf90_put_var(ncid, varid_z, tmp_array2d)
    deallocate(tmp_array2d)
    
    !!close file
    ncstatus = nf90_close(ncid)
    return
  end subroutine write_grdfile_2d

  subroutine write_grdfile_fp_2d(xmin, ymin, dx, dy, nx, ny, zval, netcdf_file, nanval)
    use netcdf
    use nrtype, only : fp
    implicit none

    real(kind = fp),    intent(in)    :: xmin, ymin, dx, dy
    integer,            intent(in)    :: nx, ny
    real(kind = fp),    intent(in)    :: zval(1 : nx, 1 : ny)
    character(len = *), intent(inout) :: netcdf_file
    real(kind = fp),    intent(in), optional :: nanval

    integer :: i, j, ncstatus, ncid, dimid_x, dimid_y, varid_x, varid_y, varid_z
    real(kind = fp), allocatable :: tmp_array(:)
    real(kind = fp), allocatable :: tmp_array2d(:, :)

    !!open file
    ncstatus = nf90_create(trim(netcdf_file), NF90_NETCDF4, ncid)

    !!dimensions
    ncstatus = nf90_def_dim(ncid, 'x', nx, dimid_x)
    ncstatus = nf90_def_dim(ncid, 'y', ny, dimid_y)

    !!variables
    ncstatus = nf90_def_var(ncid, 'x', NF90_DOUBLE, [dimid_x], varid_x)
    ncstatus = nf90_def_var(ncid, 'y', NF90_DOUBLE, [dimid_y], varid_y)
    ncstatus = nf90_def_var(ncid, 'z', NF90_FLOAT,  [dimid_x, dimid_y], varid_z)

    allocate(tmp_array(1 : nx))
    do i = 1, nx
      tmp_array(i) = xmin + dx * real(i - 1, kind = fp)
    enddo
    ncstatus = nf90_put_var(ncid, varid_x, tmp_array)
    deallocate(tmp_array)

    allocate(tmp_array(1 : ny))
    do i = 1, ny
      tmp_array(i) = ymin + dy * real(i - 1, kind = fp)
    enddo
    ncstatus = nf90_put_var(ncid, varid_y, tmp_array)
    deallocate(tmp_array)
    
    allocate(tmp_array2d(1 : nx, 1 : ny))
    do j = 1, ny
      do i = 1, nx
        tmp_array2d(i, j) = zval(i, j)
        if(present(nanval) .and. zval(i, j) .eq. nanval) tmp_array2d(i, j) = 0.0_fp / 0.0_fp
      enddo
    enddo
    ncstatus = nf90_put_var(ncid, varid_z, tmp_array2d)
    deallocate(tmp_array2d)
    
    !!close file
    ncstatus = nf90_close(ncid)
    return
  end subroutine write_grdfile_fp_2d


  subroutine get_varinfo(fileid, varname, varid, varlen)
    use netcdf 
    implicit none
    integer,            intent(in) :: fileid
    character(len = *), intent(in) :: varname
    integer, intent(out), optional :: varid
    integer, intent(out), optional :: varlen

    integer :: ncstatus, id_tmp

    if(present(varid)) then
      ncstatus = nf90_inq_varid(fileid, varname, varid)
    endif 
    if(present(varlen)) then
      ncstatus = nf90_inq_dimid(fileid, varname, id_tmp)
      ncstatus = nf90_inquire_dimension(fileid, id_tmp, len = varlen)
    endif

    return
  end subroutine get_varinfo

end module grdfile_io

