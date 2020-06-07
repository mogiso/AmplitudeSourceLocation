module read_grdfile
  private
  public :: read_grdfile_2d

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


  subroutine get_varinfo(fileid, varname, varid, varlen)
    use netcdf 
    implicit none
    integer,            intent(in)  :: fileid
    character(len = *), intent(in)  :: varname
    integer, intent(out), optional  :: varid
    integer, intent(out), optional  :: varlen

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

end module read_grdfile
