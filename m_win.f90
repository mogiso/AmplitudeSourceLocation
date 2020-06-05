!-------------------------------------------------------------------------------------------------!
!> Read Win/Win32-formatted seismograph data
!!
!! @copyright
!! Copyright (c) 2019 Takuto Maeda. All rights reserved. 
!!
!! @license 
!! This software is released under the MIT license. See LICENSE for details. 
!--
module m_win
  
  use, intrinsic :: iso_fortran_env
  use m_util
  implicit none
  private

  !! main subroutines for win/win32 files
  public :: win__read_files
  public :: win__read_file

  !! supporting public procedures & header type
  public :: win__start_file_buf
  public :: win__finish_file_buf
  public :: win__decode_buf
  public :: win__scan_buf
  public :: win__get_chid
  public :: win__ach2ich
  public :: win__ich2ach
  public :: win__init
  public :: win__free
  public :: win__hdr
  
  !! win header
  type win__hdr    
    integer              :: nb                       !< number of second blocks
    integer, allocatable :: pb(:)                    !< pointer to the second-blocks
    integer, allocatable :: yr(:), mo(:), dy(:)      !< date of the second-blocks
    integer, allocatable :: hr(:), mi(:), sc(:)      !< time of the second-blocks
    integer              :: fmt                      !< WIN1(1) or WIN32(32)
    integer              :: sz                       !< file size in bytes
    logical              :: chkbuf = .false. 
  end type win__hdr

  !! private parameters
  integer, private, parameter :: NS_MAX = 200 ! maximum samples per second
  integer, private, parameter :: WIN1  =  1
  integer, private, parameter :: WIN32 = 32
  integer, private, parameter :: NB_BUF_INIT = 3600 !< initial size of the buffer memory (3600 s)
  integer, private, parameter :: NO_CH = -99999
  
  logical,        save, private :: initialized = .false.
  integer(int16), save, private :: twelve_bit
  character(1),   save, private :: azero           !< binary image of zero in single-bit integer


  integer, parameter :: YR_MIN = 1913 !! avoid misreading 201912 as 1912
  integer, parameter :: YR_MAX = 2100

contains

  !----------------------------------------------------------------------------------------------!
  !> read single win data file for given channel ids chid. 
  !! if specified channel is not available, returns with npts(ich) = 0
  !--
  subroutine win__read_file(fn_win, chid, sfreq, nsec, tim, dat, npts)

    character(*), intent(in)  :: fn_win
    character(4), intent(in)  :: chid(:)
    integer,      intent(out) :: sfreq(:)
    integer,      intent(out) :: nsec
    integer,      intent(out) :: tim
    integer,      intent(out), allocatable :: dat(:,:)
    integer,      intent(out), allocatable :: npts(:,:)
    !--
    integer :: io_win
    integer :: j
    integer, allocatable :: dat0(:,:)
    integer :: nch
    character, allocatable :: buf(:)
    type(win__hdr) :: wh
    !----

    nch = size(chid(:))

    if( .not. initialized ) call win__init()
    
    write(error_unit,'(A)') '[win__read_file]: ' // trim(fn_win)

    call win__start_file_buf(fn_win, wh, buf, io_win )
    call win__finish_file_buf(io_win)    
    call win__scan_buf(wh, buf)
    if( .not. allocated(npts) ) allocate(npts(wh%nb, nch))
    npts(:,:) = 0
    call win__decode_buf(wh, buf, nch, chid, dat0, npts, sfreq)
    nsec = wh%nb
    allocate(dat(maxval(sfreq)*nsec,nch))

    deallocate(buf)

    call util__timelocal(wh%yr(1), wh%mo(1), wh%dy(1), wh%hr(1), wh%mi(1), wh%sc(1), tim)

    do j=1, nch
      dat(1:sfreq(j)*nsec,j) = dat0(1:sfreq(j)*nsec,j)
    end do
    deallocate(dat0)

    call win__free(wh)

  end subroutine win__read_file

  !----------------------------------------------------------------------------------------------!
  !> read multiple win data files having same temporal lengths
  !--
  subroutine win__read_files(fn_win, chid, sfreq, nsec, tim, dat, npts)
    
    character(*), intent(in)  :: fn_win(:)
    character(4), intent(in)  :: chid(:)
    integer,      intent(out) :: sfreq(:)
    integer,      intent(out) :: nsec        !< number of seconds
    integer,      intent(out) :: tim         !< Unix (POSIX) time of the head of the first file
    integer,      intent(out), allocatable :: dat(:,:)
    integer,      intent(out), allocatable :: npts(:,:)
    !--
    integer :: nw
    type(win__hdr) :: wh1, wh2
    integer :: io_win1, io_win2
    integer :: i, j
    integer, allocatable :: dat0(:,:), npts0(:,:)
    integer :: nb
    integer :: nch
    character, allocatable :: buf1(:), buf2(:)
    logical :: first_touch
    integer :: sfreq_max
    !----

    nw = size(fn_win(:))
    nch = size(chid(:))
    if( .not. initialized ) call win__init()
    write(error_unit,'(A)') '[win__read_files]: ' // trim(fn_win(1))

    first_touch = .true. 

    call win__start_file_buf(fn_win(1), wh1, buf1, io_win1 )
    call win__finish_file_buf(io_win1)


    do i=2, nw-1, 2
      write(error_unit,'(A)') '[win__read_files]: ' // trim(fn_win(i))
      call win__start_file_buf(fn_win(i), wh2, buf2, io_win2 )
      call win__scan_buf(wh1, buf1)  
      if(.not. allocated(npts0)) allocate(npts0(wh1%nb,nch))
      call win__decode_buf(wh1, buf1, nch, chid, dat0, npts0, sfreq)

      if( first_touch ) then
        nb = wh1%nb
        nsec = nb * nw
        sfreq_max = maxval(sfreq)
        allocate(dat(sfreq_max*nsec,nch))
        dat(:,:) = 0
        call util__timelocal(wh1%yr(1), wh1%mo(1), wh1%dy(1), wh1%hr(1), wh1%mi(1), wh1%sc(1), tim)
        allocate(npts(nsec,nch))
        npts(:,:) = 0
        first_touch = .false. 
      else 
        if( maxval(sfreq) > sfreq_max ) then
          call expand_dat(dat, nch, sfreq_max*nsec, maxval(sfreq)*nb*nw)
          sfreq_max = maxval(sfreq)
        end if
      end if

      npts((i-2)*nb+1:(i-1)*nb,:) = npts0(1:nb,:)

      do j=1, nch
        dat((i-2)*sfreq(j)*nb+1:(i-1)*sfreq(j)*nb, j) = dat0(1:sfreq(j)*nb,j)
      end do
      call win__free(wh1)  
              
      call win__finish_file_buf(io_win2)

      write(error_unit,'(A)') '[win__read_files]: ' // trim(fn_win(i+1))
      call win__start_file_buf(fn_win(i+1), wh1, buf1, io_win1 )
      
      call win__scan_buf(wh2, buf2)
      call win__decode_buf(wh2, buf2, nch, chid, dat0, npts0, sfreq)
     
      deallocate(buf2)

      if( maxval(sfreq) > sfreq_max ) then
        call expand_dat(dat, nch, sfreq_max*nsec, maxval(sfreq)*nb*nw)
        sfreq_max = maxval(sfreq)
      end if

      npts((i-1)*nb+1:i*nb,:) = npts0(1:nb,:)
      do j=1, nch
        dat((i-1)*sfreq(j)*nb+1:i*sfreq(j)*nb, j) = dat0(1:sfreq(j)*nb,j)
      end do
      call win__free(wh2)  

      call win__finish_file_buf(io_win1)

    end do

    if( mod(nw, 2) == 0 ) then
        write(error_unit,'(A)') '[win__read_files]: ' // trim(fn_win(nw))
        call win__start_file_buf(fn_win(nw), wh2, buf2, io_win2 )
      i = nw 
    else
      i = nw + 1 
    end if

    if( allocated(dat0) ) deallocate(dat0)

    call win__scan_buf(wh1, buf1)
    if(.not. allocated(npts0)) allocate(npts0(wh1%nb,nch))
    call win__decode_buf(wh1, buf1, nch, chid, dat0, npts0, sfreq)
    deallocate(buf1)
    if( first_touch ) then
      nb = wh1%nb
      nsec = nb * nw
      sfreq_max = maxval(sfreq)
      allocate(dat(sfreq_max*nsec,nch))
      dat(:,:) = 0
      call util__timelocal(wh1%yr(1), wh1%mo(1), wh1%dy(1), wh1%hr(1), wh1%mi(1), wh1%sc(1), tim)
      allocate(npts(nsec,nch))
      npts(:,:) = 0
      first_touch = .false. 
    else 
      if( maxval(sfreq) > sfreq_max ) then
        call expand_dat(dat, nch, sfreq_max*nsec, maxval(sfreq)*nb*nw)
        sfreq_max = maxval(sfreq)
      end if
    end if
    npts((i-2)*nb+1:(i-1)*nb,:) = npts0(1:nb,:)
    do j=1, nch
      dat((i-2)*sfreq(j)*nb+1:(i-1)*sfreq(j)*nb, j) = dat0(1:sfreq(j)*nb,j)
    end do
    call win__free(wh1)  
    if( mod(nw, 2) == 0 ) then
      call win__finish_file_buf(io_win2)
      call win__scan_buf(wh2, buf2)
      call win__decode_buf(wh2, buf2, nch, chid, dat0, npts0, sfreq)
      deallocate(buf2)

      if( maxval(sfreq) > sfreq_max ) then
        call expand_dat(dat, nch, sfreq_max*nsec, maxval(sfreq)*nb*nw)
        sfreq_max = maxval(sfreq)
      end if

      npts((i-1)*nb+1:i*nb,:) = npts0(1:nb,:)
      do j=1, nch
        dat((i-1)*sfreq(j)*nb+1:i*sfreq(j)*nb, j) = dat0(1:sfreq(j)*nb,j)
      end do
      call win__free(wh2)  
    end if

    deallocate(dat0, npts0)

  end subroutine win__read_files
  !----------------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------------!
  !> initialize module
  !--
  subroutine win__init()

    integer       :: iptr
    integer(int8) :: izero
    !----

    !! this variable will be used for reading 12-bit-width integer data from buffer
    twelve_bit = 0

    do iptr=0, 11
      twelve_bit = ibset(twelve_bit,iptr)
    end do

    izero = 0
    azero = transfer( izero, azero )
    
    initialized = .true.

  end subroutine win__init
  !----------------------------------------------------------------------------------------------! 
  
  !----------------------------------------------------------------------------------------------!
  !> Open WIN1/WIN32 file and start reading all data into buffer buf(:)
  !! The buffer memory will be allocated.
  !! If buf(:) are already being used, these data will be destoroied
  !--
  subroutine win__start_file_buf(fn, wh, buf, io)
    
    character(*),           intent(in)    :: fn     !< input filename
    type(win__hdr),         intent(inout) :: wh     !< header data
    character, allocatable, intent(inout) :: buf(:)
    integer, intent(out) :: io
    !--
    integer :: ierr
    !----

    !! initialize module
    if( .not. initialized ) call win__init()
    
    open(newunit=io, file=fn, access='stream', iostat=ierr, action='read', asynchronous='yes')
    if(ierr == 0) then
      inquire(io, size=wh%sz)
      if(allocated(buf)) then
        deallocate(buf)
      end if
      allocate(buf(wh%sz))
      read(io, asynchronous='yes') buf(:)
    else
      write(error_unit,*) 'file open error: ' // trim(fn)
      stop
    end if
    
  end subroutine win__start_file_buf
  !----------------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------------!
  !> Finalize reading process (win__start_file_buf) and close the file
  !--
  subroutine win__finish_file_buf(io)
    
    integer, intent(in) :: io
    logical :: is_open
    !----

    inquire(io, opened=is_open)
    if( is_open ) close(io)

  end subroutine win__finish_file_buf
  !----------------------------------------------------------------------------------------------!
  
  !----------------------------------------------------------------------------------------------!  !>
  !> Check the data buffer buf to
  !!   - investigate file type (WIN1 or WIN32)
  !!   - allocate necessary memory space based on the number of second-blocks
  !!   - create channel ID tables contained in the data buffer
  !! This routine must be called before start reading data
  !--
  subroutine win__scan_buf(wh, buf)
    
    type(win__hdr), intent(inout) :: wh     !< win header
    character,      intent(in)    :: buf(:) !< buffer
    !--
    integer(int16) :: idb
    integer :: sz_block
    integer, allocatable :: pb(:)  !! address of second blocks
    integer :: nb_buf
    integer :: n
    integer :: block_header_size
    integer :: dt1, dt2 !< location of date&time record
    integer :: bs1, bs2 !< location of blocksize record
    integer, allocatable :: yr(:), mo(:), dy(:), hr(:), mi(:), sc(:)
    integer :: ichk
    integer :: yr_32, mo_32, dy_32, hr_32, mi_32, sc_32
    !----
    
    !! initialize module
    if( .not. initialized ) call win__init()

    !! this routine is called only once per file
    if( wh%chkbuf ) return
    
    !! initial size of buffer memory. Expanded if necessary 
    nb_buf = NB_BUF_INIT
    allocate( pb(nb_buf), yr(nb_buf), mo(nb_buf), dy(nb_buf), hr(nb_buf), mi(nb_buf), sc(nb_buf) )
    
    !! file type check: assume win32-type and evaluate if date & time are valid
    ichk = int( transfer(buf(1:2), idb) )
    call decode_datetime(32, buf(5:12), yr_32, mo_32, dy_32, hr_32, mi_32, sc_32)
    if( YR_MIN <= yr_32 .and. yr_32 <= YR_MAX .and. &
             1 <= mo_32 .and. mo_32 <= 12     .and. &
             1 <= dy_32 .and. dy_32 <= 31     .and. &
             0 <= hr_32 .and. hr_32 <= 23     .and. &
             0 <= mi_32 .and. mi_32 <= 59     .and. ichk == 0 ) then
      wh%fmt = WIN32
      pb(1) = 5
    else ! if not, it must be WIN1 file
      wh%fmt = WIN1
      pb(1) = 1
    end if

    if( wh%fmt == WIN1 ) then
      block_header_size = 0
      dt1 = 4
      dt2 = 9
      bs1 = 0
      bs2 = 3
    else
      block_header_size = 16
      dt1 = 0
      dt2 = 7
      bs1 = 12
      bs2 = 15
    end if

    n = 0
    !! obtain block locations & datetime information
    do
      
      n = n + 1

      call decode_datetime(wh%fmt, buf(pb(n)+dt1:pb(n)+dt2), &
                           yr(n), mo(n), dy(n), hr(n), mi(n), sc(n))

      sz_block = transfer(buf(pb(n)+bs2:pb(n)+bs1:-1), sz_block )

      !! location of the next block
      pb(n+1) = pb(n) + block_header_size + sz_block

      !! EOF detection
      if( pb(n+1)  > wh% sz ) then
        wh%nb = n
        exit
      end if
      
      !! expand buffer memory if necessary
      if( n == nb_buf - 1) then
        call expand_i( pb, nb_buf )
      end if
      
    end do
    
    allocate(wh%pb(wh%nb+1) )
    allocate(wh%yr(wh%nb), wh%mo(wh%nb), wh%dy(wh%nb), wh%hr(wh%nb), wh%mi(wh%nb), wh%sc(wh%nb))
    
    wh%pb(1:wh%nb+1) = pb(1:wh%nb+1)
    wh%yr(1:wh%nb)   = yr(1:wh%nb)
    wh%mo(1:wh%nb)   = mo(1:wh%nb)
    wh%dy(1:wh%nb)   = dy(1:wh%nb)
    wh%hr(1:wh%nb)   = hr(1:wh%nb)
    wh%mi(1:wh%nb)   = mi(1:wh%nb)
    wh%sc(1:wh%nb)   = sc(1:wh%nb)
    
    deallocate(pb, yr, mo, dy, hr, mi, sc)

    wh%chkbuf = .true. 
    
  end subroutine win__scan_buf
  !----------------------------------------------------------------------------------------------!  

  !----------------------------------------------------------------------------------------------!
  !> Return a list of channel IDs in the ib-th block of the file buffer
  !--
  subroutine win__get_chid( wh, buf, ib, nch, chid )
    
    type(win__hdr), intent(inout) :: wh
    character,      intent(in)    :: buf(:)
    integer,        intent(in)    :: ib     !< block number
    integer,        intent(out)   :: nch
    !--
    character(4), allocatable :: chid(:)
    character(4), allocatable :: chid_tmp(:)
    integer :: p
    integer :: nmaxch
    integer(int8) :: isb
    integer(int16) :: idb
    integer :: ss, ns
    integer :: p1, p2
    !----

    !! initialize module
    if( .not. initialized ) call win__init()

    if( .not. wh%chkbuf ) call win__scan_buf(wh, buf)
    
    if( wh%fmt == WIN32 ) then
      p1 = 16 !< block header size
      p2 = 2  !< institution ID & network ID 
    else
      p1 = 10
      p2 = 0  !< no institution&network IDs in WIN1
    end if
    
    if( ib > wh%nb ) then
      write(error_unit,*) "ERROR [win__chid]: No such block ", ib
      return
    end if
    
    nmaxch = 300
    allocate(chid_tmp(nmaxch))
    if( allocated( chid ) ) deallocate(chid)
    nch = 0
    p = wh%pb(ib) + p1

    do while ( p < wh%pb(ib+1) )

      !! temp memory expand
      if( nch == nmaxch ) call expand_a( 4, chid_tmp, nmaxch )

      ! read channel ID
      p = p + p2
      write(chid_tmp(nch+1),'(Z4)') transfer( buf(p+1:p:-1), idb )
      isb = transfer( buf(p+2), isb )
      ss = int( ishft(isb, -4) )
      ns = int( iand( transfer( buf(p+3:p+2:-1), idb), twelve_bit ))
      
      nch = nch + 1
      if( ss == 0 ) then
        p = p + 8 + ns/2
      else
        p = p + 8 + (ns-1) * ss
      end if
    end do

    allocate( chid(nch) )
    chid(1:nch) = chid_tmp(1:nch)
    deallocate(chid_tmp)

    !! zero padding
    do p1=1, nch
      do p2=1, 4
        if( chid(p1)(p2:p2) == ' ' ) chid(p1)(p2:p2) = '0'
      end do
    end do
    
  end subroutine win__get_chid
  !----------------------------------------------------------------------------------------------! 

  !----------------------------------------------------------------------------------------------!
  !> Decode win data buffer for specified channel ids
  !--
  subroutine win__decode_buf(wh, buf, nch, chid, dat, npts, sfreq)

    type(win__hdr), intent(inout)              :: wh
    character,        intent(in)               :: buf(:)            !< win data buffer
    integer,          intent(in)               :: nch               !< #channels to decode
    character(4),     intent(in)               :: chid(nch)         !< channel IDs 
    integer,          intent(out), allocatable :: dat(:,:)          !< (npts,nch)
    integer,          intent(out)              :: npts(wh%nb, nch)  !< #samples. 0 for not-found
    integer,          intent(out)              :: sfreq(nch)        !< sampling frequency
    !--
    integer :: ib
    integer :: ns(nch)
    integer :: dbuf(NS_MAX,nch)
    integer :: ich
    integer(int16) :: ichid(nch)
    integer :: chid_tbl(-32768:32767)
    integer :: i
    integer :: sfreq_max
    !----

    !! initialize module
    if( .not. initialized ) call win__init()
    if( .not. wh%chkbuf ) call win__scan_buf(wh, buf)
    do ich=1, nch
      ichid(ich) = win__ach2ich( chid(ich) )
    end do
    
    ns(:) = -1
    npts(:,:) = 0
    
    dbuf = 0
    sfreq(:) = 0
    !! set channel ID inverse table
    chid_tbl(:) = NO_CH
    do i=1, size(ichid)
      chid_tbl(ichid(i)) = i
    end do

    do ib=1, wh%nb
      call decode_1sec( wh, buf, ib, ichid, chid_tbl, ns, dbuf(:,:))
      if( ib == 1 ) then
        sfreq_max = max(maxval(ns(:)),1)
        if(allocated(dat)) deallocate(dat)
        allocate(dat(sfreq_max*wh%nb, nch))
        dat(:,:) = 0
      else
        if( maxval(ns(:)) > sfreq_max ) then
          !! expand memory size
          call expand_dat(dat, nch, sfreq_max*wh%nb, maxval(ns(:))*wh%nb)
          sfreq_max = maxval(ns(:))
        end if
      end if

      do ich=1, nch
        npts(ib,ich) = ns(ich)
        if( ns(ich) > 0 ) then
          dat((ib-1)*ns(ich)+1:ib*ns(ich), ich) = dbuf(1:ns(ich),ich)
          sfreq(ich) = ns(ich)
        end if
      end do
    end do
    
  end subroutine win__decode_buf
  !----------------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------------!
  !> Convert ascii channel ID to integer
  !--
  function win__ach2ich( achid ) result(chid)
    
    character(4), intent(in) :: achid
    !--
    integer(int16) :: chid
    !----

    read(achid,'(Z4)') chid

  end function win__ach2ich
  !----------------------------------------------------------------------------------------------!
  
  !----------------------------------------------------------------------------------------------!
  !> Convert integer channel ID to ascii
  !--
  function win__ich2ach( chid ) result(achid)
    
    integer(int16), intent(in) :: chid
    !--
    character(4) :: achid
    !----

    write(achid,'(Z4.4)') chid

  end function win__ich2ach
  !----------------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------------!
  !> Free memory
  !--
  subroutine win__free( wh )
    
    type(win__hdr), intent(inout) :: wh
    !----

    deallocate(wh%pb)
    deallocate(wh%yr, wh%mo, wh%dy, wh%hr, wh%mi, wh%sc)

    wh%chkbuf = .false.
    wh%nb = 0
    wh%sz = 0

  end subroutine win__free
  !----------------------------------------------------------------------------------------------!
  
  !----------------------------------------------------------------------------------------------!
  !> Double the size of the input allocable array x(1:n) to x(1:2*n)
  !! Array size n will be overwritten to 2*n
  !--
  subroutine expand_i(x, n)

    integer, intent(inout), allocatable :: x(:)
    integer, intent(inout) :: n
    !--
    integer, allocatable :: tmp(:)
    !----

    if( size(x) /= n ) then
      write(error_unit,*) "Error [expand]: array size mismatch"
      return
    end if
    
    allocate( tmp(n) )
    tmp(:) = x(:)
    deallocate(x)
    allocate(x(2*n))
    x(1:n) = tmp(:)
    deallocate(tmp)
    n = 2 * n
    
  end subroutine expand_i
  !----------------------------------------------------------------------------------------------!    
  !----------------------------------------------------------------------------------------------!
  !> Double the size of the input allocable array x(1:n) to x(1:2*n)
  !! Array size n will be overwritten to 2*n
  !--
  subroutine expand_a(m, x, n)

    integer,      intent(in) :: m
    character(m), intent(inout), allocatable :: x(:)
    integer, intent(inout) :: n
    !--
    character(m), allocatable :: tmp(:)
    !----

    if( size(x) /= n ) then
      write(error_unit,*) "Error [expand]: array size mismatch"
      return
    end if
    
    allocate( tmp(n) )
    tmp(:) = x(:)
    deallocate(x)
    allocate(x(2*n))
    x(1:n) = tmp(:)
    deallocate(tmp)
    n = 2 * n
    
  end subroutine expand_a
  !----------------------------------------------------------------------------------------------!    

  !----------------------------------------------------------------------------------------------!    
  !> memory size change dat(npts0,nch) -> dat(npts1,nch)
  !--
  subroutine expand_dat(dat, nch, npts0, npts1)
    
    integer, intent(inout), allocatable :: dat(:,:)
    integer, intent(in)                 :: nch, npts0, npts1
    !--
    integer, allocatable :: dtmp(:,:)
    !----

    if( size(dat, dim=1) /= npts0 ) then
      write(error_unit,*) '[expand_dat]: size mismatch'
      return
    end if
    allocate(dtmp(npts0,nch))
    dtmp(:,:) = dat(:,:)
    deallocate(dat)
    allocate(dat(npts1,nch))
    dat(1:npts0,:) = dtmp(:,:)
    dat(npts0+1:npts1,:) = 0
    deallocate(dtmp)

  end subroutine expand_dat
  !----------------------------------------------------------------------------------------------!    
  !----------------------------------------------------------------------------------------------!    
  !> memory size change dat(npts0,nch) -> dat(npts1,nch)
  !--
  subroutine expand_rdat(dat, nch, npts0, npts1)
    
    real, intent(inout), allocatable :: dat(:,:)
    integer, intent(in)              :: nch, npts0, npts1
    !--
    real, allocatable :: dtmp(:,:)
    !----

    if( size(dat, dim=1) /= npts0 ) then
      write(error_unit,*) '[expand_dat]: size mismatch'
      return
    end if
    allocate(dtmp(npts0,nch))
    dtmp(:,:) = dat(:,:)
    deallocate(dat)
    allocate(dat(npts1,nch))
    dat(1:npts0,:) = dtmp(:,:)
    dat(npts0+1:npts1,:) = 0
    deallocate(dtmp)

  end subroutine expand_rdat
  !----------------------------------------------------------------------------------------------!      
  !----------------------------------------------------------------------------------------------!    
  !> Read win/win32-formatted date and time in the second block header
  !--
  subroutine decode_datetime( fmt, tbuf, yr, mo, dy, hr, mi, sc )

    integer, intent(in) :: fmt
    character, intent(in) :: tbuf(:)
    integer, intent(out) :: yr, mo, dy, hr, mi, sc
    !--
    integer(int8) :: hb1, hb2
    integer     :: y1, y2
    !----

    if( fmt == WIN1 ) then
      call bcd_int(tbuf(1),hb1,hb2); y1 = hb1*10+hb2
      call bcd_int(tbuf(2),hb1,hb2); mo = hb1*10+hb2
      call bcd_int(tbuf(3),hb1,hb2); dy = hb1*10+hb2
      call bcd_int(tbuf(4),hb1,hb2); hr = hb1*10+hb2
      call bcd_int(tbuf(5),hb1,hb2); mi = hb1*10+hb2
      call bcd_int(tbuf(6),hb1,hb2); sc = hb1*10+hb2

      !! assume observation period 1981 -- 2080
      if( y1 > 80 ) then
        yr = y1 + 1900
      else
        yr = y1 + 2000
      end if
      
    else if ( fmt == WIN32 ) then
      call bcd_int(tbuf(1),hb1,hb2); y1 = hb1*10+hb2
      call bcd_int(tbuf(2),hb1,hb2); y2 = hb1*10+hb2
      call bcd_int(tbuf(3),hb1,hb2); mo = hb1*10+hb2
      call bcd_int(tbuf(4),hb1,hb2); dy = hb1*10+hb2
      call bcd_int(tbuf(5),hb1,hb2); hr = hb1*10+hb2
      call bcd_int(tbuf(6),hb1,hb2); mi = hb1*10+hb2
      call bcd_int(tbuf(7),hb1,hb2); sc = hb1*10+hb2
      yr = y1 * 100 + y2
    end if

  end subroutine decode_datetime
  !----------------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------------!
  !> Decode bcd-formatted buffer into two integers
  !--
  subroutine bcd_int( b, j, k )

    character,   intent(in)  :: b
    integer(int8), intent(out) :: j
    integer(int8), intent(out) :: k
    !--
    integer(int8) :: i
    integer :: l
    !----
    
    i = transfer( b, i )
    j=ishft(i,-4)
    k=i
    do l=4,7
      k=ibclr(k,l)
    end do

  end subroutine bcd_int
  !----------------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------------!
  subroutine hbyte_int( i, j, k )

    integer(int8), intent(in)  :: i
    integer(int8), intent(out) :: j  ! i(1:4)
    integer(int8), intent(out) :: k  ! i(5:8)
    !--
    integer :: l
    !----

    !! Decompose
    j=ishft(i,-4)
    k=i
    do l=4,7
      k=ibclr(k,l)
    end do

    !! Positive/Negative
    if( btest(j,3) ) then
      j=ibset(j,4)
      j=ibset(j,5)
      j=ibset(j,6)
      j=ibset(j,7)
    end if

    if( btest(k,3) ) then
      k=ibset(k,4)
      k=ibset(k,5)
      k=ibset(k,6)
      k=ibset(k,7)
    end if
    
  end subroutine hbyte_int
  !----------------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------------!
  !> Read specified second-block and decode it for given channels 
  !--
  subroutine decode_1sec( wh, buf, ib, chid, chid_tbl, ns, dat )
    
    type(win__hdr), intent(in)  :: wh
    character,      intent(in)  :: buf(:)
    integer,        intent(in)  :: ib !< block number
    integer(int16), intent(in)  :: chid(:)
    integer,        intent(in)  :: chid_tbl(-32768:32767)
    integer,        intent(out) :: ns(:)
    integer,        intent(out) :: dat(:,:)
    !--
    integer(int8) :: isb
    integer(int8) :: isb1, isb2
    integer(int16) :: idb
    integer(int16) :: chid_buf
    integer :: p, q
    integer :: ss, ns0
    integer :: p1, p2
    integer :: ich
    integer :: i
    integer :: nch
    !----

    nch = size(chid)
    
    if( wh%fmt == WIN32 ) then
      p1 = 16 !< block header size
      p2 = 2  !< institution ID & network ID 
    else
      p1 = 10
      p2 = 0  !< no institution&network IDs in WIN1
    end if
    
    if( ib > wh%nb ) then
      write(error_unit,*) "ERROR [decode_1sec]: No such block ", ib
      return
    end if
    
    p = wh%pb(ib) + p1

    !! initialize
    ns(:) = -1

    do while ( p < wh%pb(ib+1) )

      !! skip institution & network ID
      p = p + p2

      !! read channel ID
      chid_buf = transfer( buf(p+1:p:-1), idb )

      !! decompose sample size ss (0.5 byte) & sample length ns (1.5 byte)
      isb = transfer( buf(p+2), isb )
      ss = int( ishft(isb, -4) )
      ns0 = int( iand( transfer( buf(p+3:p+2:-1), idb), twelve_bit ))

      ich = chid_tbl(chid_buf)

      if( ich /= NO_CH ) then

        ns(ich) = ns0

        !! first data
        dat(1,ich) = transfer( buf(p+7:p+4:-1), ich )

        select case( ss ) 
        case( 0 ) !! 0.5 byte

          do i=2, ns0-2, 2
            q = 8 + p + (i-2) / 2
            call hbyte_int( transfer(buf(q:q), isb), isb1, isb2)
            dat(i,  ich) = dat(i-1,ich) + int(isb1)
            dat(i+1,ich) = dat(i-1,ich) + int(isb1) + int(isb2)
          end do
          if( ns0 /= 1 ) then
            q = 8 + p + (ns0-2) / 2
            call hbyte_int( transfer(buf(q:q), isb), isb1, isb2)
            dat(ns0,  ich) = dat(ns0-1,ich) + int(isb1)
          end if
        case( 1 )
          do i=2, ns0
            q = 8 + p + (i-2) 
            dat(i,ich) = dat(i-1,ich) + int( transfer( buf(q:q), isb) )
          end do
        case( 2 )
          do i=2, ns0
            q = 8 + p + (i-2) * 2
            dat(i,ich) = dat(i-1,ich) + int( transfer( buf(q+1:q:-1), idb) )
          end do
        case( 3 )
          do i=2, ns0
            q = 8 + p + (i-2) * 3
            dat(i,ich) = dat(i-1,ich) + int( transfer( azero // buf(q+2)//buf(q+1)//buf(q), i ) ) / 256
          end do
        case( 4 )
          do i=2, ns0
            q = 8 + p + (i-2) * 4
            dat(i,ich) = dat(i-1,ich) + transfer(buf(q+3:q:-1), i)
          end do
        case( 5 )
          do i=2, ns0
            q = 8 + p + (i-2) * 4
            dat(i,ich) = transfer(buf(q+3:q:-1), i)
          end do
        end select
      end if

      if( ns0 /= 1 ) then
        if ( ss == 0 ) then
          p = p + 8 + ns0/2
        else
          p = p + 8 + (ns0-1) * ss
        end if
      else
        p = p + 8
      end if

    end do

  end subroutine decode_1sec
  !----------------------------------------------------------------------------------------------! 

end module m_win
!------------------------------------------------------------------------------------------------!
