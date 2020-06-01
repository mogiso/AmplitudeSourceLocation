!
! grd handling module
!
!           AUTHOR: KATO Masaya <mkato@pastel.ocn.ne.jp>
!    FIRST RELEASE: 2002-09-04
!   LATEST RELEASE: 2002-09-05
!
module GMT
  implicit none
  include 'netcdf.inc'
  private

  ! NaN value
  real    :: GMT_fnan
  real(8) :: GMT_dnan

  data GMT_fnan/Z'7fffffff'/
  data GMT_dnan/Z'7fffffffffffffff'/

  ! common parameter
  integer, parameter :: GRD_NOERROR  =  1
  integer, parameter :: GRD_NOCREATE = -1
  integer, parameter :: GRD_NOPUTVAR = -2

  integer, parameter :: OneByteInt    = selected_int_kind(2)
  integer, parameter :: TwoByteInt    = selected_int_kind(4)
  integer, parameter :: FourByteInt   = selected_int_kind(9)
  integer, parameter :: FourByteReal  = selected_real_kind( P =  6, R =  37 )
  integer, parameter :: EightByteReal = selected_real_kind( P = 13, R = 307 )
  
  public :: grd_create

  interface grd_create
     module procedure grd_create_float
  end interface
contains
  integer function grd_create_float( grdfile, var, &
       & x_range, y_range, spacing, node_offset, &
       & overwrite, NaN, iscan, jscan )
    character(len=*)                       , intent(IN) :: grdfile
    real(kind=FourByteReal), dimension(:,:), intent(IN) :: var
    real(kind=EightByteReal)               , intent(IN) :: x_range(2)
    real(kind=EightByteReal)               , intent(IN) :: y_range(2)
    real(kind=EightByteReal)               , intent(IN) :: spacing(2)
    integer(kind=FourByteInt)              , optional   :: node_offset
    integer(kind=FourByteInt)              , optional   :: iscan, jscan
    real(kind=FourByteReal)                , optional   :: NaN
    logical                                , optional   :: overwrite

    ! description of optional arguments
    !
    !      * node_offset   see Appendix B in users guide
    !                         0: grid line registration
    !                         1: pixel registration
    !      * NaN           NaN value. If var(i,j) = NaN, var(i,j) set to 'NaN'.
    !      * iscan         Scaning method in x direction
    !                         0: left  (west)  -> right (east) [ default ]
    !                         1: right (east)  -> left  (west)
    !      * jscan         Scaning method in y direction
    !                         0: up   (north) -> down (south)  [ default ]
    !                         1: down (south) -> up   (north)

    integer(kind=FourByteInt) :: ix, jx, id, status, ivar, i, j, ip, jp, iv
    integer(kind=FourByteInt) :: over
    real(kind=EightByteReal)  :: z_range(2), scale_factor, add_offset
    real(kind=FourByteReal), allocatable :: out(:)

    integer(kind=FourByteInt) :: node = 0
    integer(kind=FourByteInt) :: isep = 0
    integer(kind=FourByteInt) :: jsep = 0
    logical                   :: owrite = .false.

    if( present( node_offset ) ) node = node_offset
    if( present( iscan       ) ) isep = iscan
    if( present( jscan       ) ) jsep = jscan
    if( present( overwrite   ) ) owrite = overwrite
    if( node /= 0 .and. node /= 1 ) node = 0
    scale_factor = 1.
    add_offset   = 0.
    over         = NF_NOCLOBBER
    if( owrite ) over = NF_CLOBBER

    ix = size( var, dim=1 )
    jx = size( var, dim=2 )

    z_range(1) = minval( var )
    z_range(2) = maxval( var )

    status = write_grd_header( grdfile, id, ix, jx, NF_FLOAT, &
         & scale_factor, add_offset, node, &
         & x_range, y_range, z_range, spacing, ivar, over )
    if( status /= GRD_NOERROR ) then
       grd_create_float = GRD_NOCREATE

    else
       allocate( out(ix*jx) )
       do j = 1, jx
          jp = j
          if( jsep == 1 ) jp = jx - j + 1
          do i = 1, ix
             ip = i
             if( isep == 1 ) ip = ix - i + 1
             iv = ( j - 1 ) * ix + i
             out(iv) = var(ip,jp)
          end do
       end do
       if( present( NaN ) ) then
          where( out == NaN ) out = GMT_fnan
       end if
       status = nf_put_var_real( id, ivar, out )
       deallocate( out )
       if( status /= NF_NOERR ) then
          grd_create_float = GRD_NOPUTVAR
       else
          grd_create_float = GRD_NOERROR
       end if
    end if

    status = nf_close( id )
    return
  end function grd_create_float
  
  integer function write_grd_header( grdfile, id, ix, jx, xtype, &
       & scale_factor, add_offset, node_offset, &
       & x_range, y_range, z_range, spacing, iz, over )
    character(len=*), intent(IN)             :: grdfile
    integer(kind=FourByteInt), intent(IN)    :: ix, jx, xtype
    integer(kind=FourByteInt), intent(IN)    :: node_offset, over
    integer(kind=FourByteInt), intent(INOUT) :: id, iz
    real(kind=EightByteReal) , intent(IN)    :: scale_factor, add_offset
    real(kind=EightByteReal)               , intent(IN) :: x_range(2)
    real(kind=EightByteReal)               , intent(IN) :: y_range(2)
    real(kind=EightByteReal)               , intent(IN) :: z_range(2)
    real(kind=EightByteReal)               , intent(IN) :: spacing(2)

    integer :: status
    integer :: iside, ixysize
    integer :: ix_range, iy_range, iz_range, ispacing, idimension

    ! create grd file
    status = nf_create( grdfile, over, id )
    if( status /= NF_NOERR ) then
       write_grd_header = GRD_NOCREATE
       write(6,*) nf_strerror( status )
    else
       write_grd_header = GRD_NOERROR
    end if

    ! define of dimensions
    status = nf_def_dim( id, 'side'  , 2    , iside   )
    status = nf_def_dim( id, 'xysize', ix*jx, ixysize )

    ! define of variables
    status = nf_def_var( id, 'x_range'  , NF_DOUBLE, 1, iside  , ix_range   )
    status = nf_def_var( id, 'y_range'  , NF_DOUBLE, 1, iside  , iy_range   )
    status = nf_def_var( id, 'z_range'  , NF_DOUBLE, 1, iside  , iz_range   )
    status = nf_def_var( id, 'spacing'  , NF_DOUBLE, 1, iside  , ispacing   )
    status = nf_def_var( id, 'dimension', NF_INT   , 1, iside  , idimension )
    status = nf_def_var( id, 'z'        , xtype    , 1, ixysize, iz         )

    ! attributes of each variables
    status = nf_put_att_text( id, ix_range, 'units', 11, 'user_x_unit' )
    status = nf_put_att_text( id, iy_range, 'units', 11, 'user_y_unit' )
    status = nf_put_att_text( id, iz_range, 'units', 11, 'user_z_unit' )
    status = nf_put_att_double( id, iz, 'scale_factor', &
         & NF_DOUBLE, 1, scale_factor )
    status = nf_put_att_double( id, iz, 'add_offset', &
         & NF_DOUBLE, 1, add_offset   )
    status = nf_put_att_int( id, iz, 'node_offset', NF_INT, 1, node_offset )

    ! global attributes
    status = nf_put_att_text( id, NF_GLOBAL, 'title' , 1, ' ' )
    status = nf_put_att_text( id, NF_GLOBAL, 'source', 21, &
         & 'created by GMT module' )

    ! end define mode
    status = nf_enddef( id )

    ! put range informations
    status = nf_put_var_double( id, ix_range  , x_range    )
    status = nf_put_var_double( id, iy_range  , y_range    )
    status = nf_put_var_double( id, iz_range  , z_range    )
    status = nf_put_var_double( id, ispacing  , spacing    )
    status = nf_put_var_int( id, idimension, (/ix, jx/) )

  end function write_grd_header
end module GMT
