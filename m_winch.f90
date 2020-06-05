!-------------------------------------------------------------------------------------------------!
!> WIN/WIN32 Channel Table
!!
!! @copyright
!! Copyright (c) 2019 Takuto Maeda. All rights reserved. 
!!
!! @license 
!! This software is released under the MIT license. See LICENSE for details. 
!--
module m_winch

  use, intrinsic :: iso_fortran_env
  use m_win
  implicit none

  public :: winch__hdr

  type winch__hdr

    character(4)    :: achid     !< channel id in 4-digit character
    integer(int16)  :: ichid     !< channel id in integer form
    character(16)   :: stnm      !< station name
    character(16)   :: cmpnm     !< component name
    character(16)   :: unit      !< measurement unit
    real(real64)    :: period    !< natural period of the sensor
    real(real64)    :: damp      !< damping constant
    real(real64)    :: lon       !< station longitude (degrees-east)
    real(real64)    :: lat       !< station latitude (degrees-north)
    real(real64)    :: elev      !< station elevation (m)
    real(real64)    :: sens      !< sensor sensitivity
    integer         :: ampl      !< amplification factor
    real(real64)    :: step      !< 1-digit step in voltage
    real(real64)    :: conv      !< conversion coef. from DC to physical unit

  end type winch__hdr

contains

  !----------------------------------------------------------------------------------------------!
  !> Initialize win channel header
  !--
  subroutine winch__init(wch)

    type(winch__hdr), intent(inout) :: wch
    !----

    wch%achid    = '-123'         
    wch%ichid    = 0              
    wch%stnm     = '-12345.0'     
    wch%cmpnm    = '-12345.0'     
    wch%unit     = '-12345.0'     
    wch%ampl     = -12345         
    wch%period   = -12345.0_real64
    wch%damp     = -12345.0_real64
    wch%lon      = -12345.0_real64
    wch%lat      = -12345.0_real64
    wch%elev     = -12345.0_real64
    wch%sens     = -12345.0_real64
    wch%step     = -12345.0_real64
    wch%conv     = -12345.0_real64

  end subroutine winch__init
  !----------------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------------!
  !> Read the spacified channel table file
  !--
  subroutine winch__read_tbl(fn_chtbl, wch)

    character(*), intent(in) :: fn_chtbl
    type(winch__hdr), intent(out), allocatable :: wch(:)
    !! --
    integer        :: io
    integer(int16) :: ich
    integer        :: ierr
    character(512) :: line
    type(winch__hdr) :: wch0
    !! ----

    allocate(wch(0))
    
    open(newunit=io, file=trim(fn_chtbl), iostat=ierr, status='old', action='read')
    if( ierr /= 0 ) then
       write(error_unit,'(A)') 'error [winch__readchtbl]: file not found.'
       return
    end if
        
    do 
       read(io, '(A)', iostat=ierr) line
       if( ierr/= 0 ) exit

       line = adjustl(line)

       ! skip comment or blank lines
       if( line(1:1) .eq. "#" ) cycle
       if( len_trim(line) == 0 ) cycle

       call winch__init(wch0)

       ich = win__ach2ich(line(1:4))

       ! store the channel information
       wch0 % ichid = ich
       wch0 % achid = line(1:4)

       ! the space (separation) should be searched via index function rather than automated read
       ! command, as the channel table info may contain other separator characters such as '/'

       call readline_skip()
       call readline_skip()
       call readline_skip()

       call readline_ascii  (wch0%stnm)                     !  station name
       call readline_ascii  (wch0%cmpnm)                    !  station name
 
       call readline_skip()
       call readline_skip()

       call readline_real64 (wch0%sens, 1.0_real64)         ! sensor sensitivity
       call readline_ascii  (wch0%unit)                     ! sensor unit
       call readline_real64 (wch0%period, -12345.0_real64)  ! natural period
       call readline_real64 (wch0%damp,   -12345.0_real64)  ! damping constant
       call readline_integer(wch0%ampl,   -12345)           ! amplification factor 
       call readline_real64 (wch0%step,   1.0_real64)       ! 1 digit step in voltage
       call readline_real64 (wch0%lat,    -12345.0_real64)  ! station latitude
       call readline_real64 (wch0%lon,    -12345.0_real64)  ! station longitude
       call readline_real64 (wch0%elev,   -12345.0_real64)  ! station elevation

       ! conversion coefficient
       wch0%conv = wch0%step / (wch0%sens * 10**(dble(wch0%ampl)/20.0_real64))

       !! expand array
       wch = [wch, wch0]

    end do

  contains

    subroutine readline_skip()
      integer :: idx
      !----
      idx=index(line,' ')
      line=adjustl(line(idx:len_trim(line)))
      
    end subroutine readline_skip
 
    subroutine readline_ascii(var)
      character(*), intent(out) :: var
      integer :: idx
      !----
      idx=index(line,' ')
      read(line(1:idx),'(A)') var
      line=adjustl(line(idx:len_trim(line)))
      
    end subroutine readline_ascii

    subroutine readline_real64(var, default_var)
      real(real64), intent(out) :: var
      real(real64), intent(in)  :: default_var
      integer :: idx

      idx=index(line,' ')
      if( trim(line(1:idx)) == "*" ) then
        var= default_var
      else
        read(line(1:idx),*) var
      end if
      line=adjustl(line(idx:len_trim(line)))

    end subroutine readline_real64
  
    subroutine readline_integer(var, default_var)
      integer, intent(out) :: var
      integer, intent(in)  :: default_var
      integer :: idx

      idx=index(line,' ')
      if( trim(line(1:idx)) == "*" ) then
         var = default_var
      else
         read(line(1:idx),*) var
      end if
      line=adjustl(line(idx:len_trim(line)))
    end subroutine readline_integer

  end subroutine winch__read_tbl
  !----------------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------------!
  !> Return all channel IDs in the list
  !--
  subroutine winch__get_all_chid(wch, chid)

    type(winch__hdr), intent(in) :: wch(:)
    character(4), intent(out), allocatable :: chid(:)

    allocate(chid(size(wch)))
    chid(:) = wch(:)%achid

  end subroutine winch__get_all_chid

  !----------------------------------------------------------------------------------------------!
  !> Return all station names in the channel list without duplications
  !--
  subroutine winch__get_all_stnm(wch, stnm)
    
    type(winch__hdr), intent(in) :: wch(:)
    character(*), intent(out), allocatable :: stnm(:)
    !--
    integer :: i, j
    character(len(stnm)) :: stnm0
    logical :: found
    integer :: nch, nst
    !----

    nch = size(wch)

    if( nch <= 0 ) then
      nst = 0
      return
    end if

    nst = 0
    do i=1, nch
      if( .not. allocated(stnm) ) then
        nst = 1
        allocate(stnm(nst))
        stnm(nst) = trim(adjustl(wch(i)%stnm))
      end if

      stnm0 = trim(adjustl(wch(i)%stnm))
      found = .false.

      do j=1, nst
        found = trim(adjustl(stnm0)) == trim(adjustl(stnm(j)))
        if( found ) exit
      end do

      if( .not. found ) then
        nst = nst + 1
        stnm = [stnm, stnm0] !! expand the array
      end if

    end do
    
  end subroutine winch__get_all_stnm
  !----------------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------------!
  !> Return all station names in the channel list without duplications
  !--
  subroutine winch__get_all_cmpnm(wch, cmpnm)
    
    type(winch__hdr), intent(in) :: wch(:)
    character(*), intent(out), allocatable :: cmpnm(:)
    !--
    integer :: i, j
    character(len(cmpnm)) :: cmpnm0
    logical :: found
    integer :: nch, ncmp
    !----

    nch = size(wch)

    if( nch <= 0 ) then
      ncmp = 0
      return
    end if

    ncmp = 0
    do i=1, nch
      if( .not. allocated(cmpnm) ) then
        ncmp = 1
        allocate(cmpnm(ncmp))
        cmpnm(ncmp) = trim(adjustl(wch(i)%cmpnm))
      end if

      cmpnm0 = trim(adjustl(wch(i)%cmpnm))
      found = .false.

      do j=1, ncmp
        found = trim(adjustl(cmpnm0)) == trim(adjustl(cmpnm(j)))
        if( found ) exit
      end do

      if( .not. found ) then
        ncmp = ncmp + 1
        cmpnm = [cmpnm, cmpnm0] !! expand the array
      end if

    end do
    
  end subroutine winch__get_all_cmpnm
  !----------------------------------------------------------------------------------------------!


  !----------------------------------------------------------------------------------------------!
  !> Return channel ID corresponding to given station and component name
  !--
  subroutine winch__st2chid(wch, stnm, cmpnm, chid, ikey)

    type(winch__hdr), intent(in)  :: wch(:)
    character(*), intent(in)  :: stnm  !< station name
    character(*), intent(in)  :: cmpnm !< component name
    character(4), intent(out) :: chid  !< channel id
    integer,      intent(out) :: ikey  !< location in the wch(:), -1 for failure
    !--
    integer :: i
    !----

    do i=1, size(wch)
      if( trim(adjustl( wch(i)%stnm  )) == trim(adjustl(stnm)) .and. &
          trim(adjustl( wch(i)%cmpnm )) == trim(adjustl(cmpnm)) ) then
        chid = wch(i)%achid
        ikey = i
        return
       end if
       
    end do

    ! not found
    ikey = -1
    chid = 'XXXX'
    return

  end subroutine winch__st2chid
  !----------------------------------------------------------------------------------------------!  

 
end module m_winch
!-------------------------------------------------------------------------------------------------!