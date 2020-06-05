!-------------------------------------------------------------------------------------------------!
!> Utility routines
!!
!! @copyright
!! Copyright (c) 2019 Takuto Maeda. All rights reserved. 
!!
!! @license 
!! This software is released under the MIT license. See LICENSE for details. 
!<
!--
module m_util

  use iso_fortran_env
  implicit none
  public

contains

  !-----------------------------------------------------------------------------------------------!
  !> Read command line buffer
  !! given buffer (command-line option) may be one of the following:
  !!  - a filename contains data
  !!  - 'all' or 'ALL'
  !!  - comma-separated data
  !--
  subroutine util__read_arglst( buf, n, is_all, lst )

    character(*),              intent(in)  :: buf
    integer,                   intent(out) :: n
    logical,                   intent(out) :: is_all
    character(*), allocatable, intent(out) :: lst(:)

    !--
    logical :: is_exist
    integer :: io
    integer :: i
    !----
    
    is_all = .false.

    !! first assume buf is a filename
    inquire( file=buf, exist=is_exist )

    if( is_exist ) then

      !! if it is filename, read list from it     
      open(newunit=io, file=buf, action='read', status='old')

      call util__countline(io, n)

      allocate( lst(n) )
      do i=1, n
        read(io,'(A)') lst(i)
      end do
      close(io)
      
      is_all = .false.

    else if( trim(buf) == 'ALL' .or. trim(buf) == 'all' ) then
      
      is_all = .true.
      n = 0
      
    else
      
      !! buf should be a comma-separated channel number list 
      call util__split( buf, ',', n, lst )
      is_all = .false.
      
    end if

  end subroutine util__read_arglst
  !----------------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------------!
  subroutine util__split(buf, sep, n, nm)
    
    character(*), intent(in)  :: buf  !< input buffer
    character(*), intent(in)  :: sep  !< separator
    integer,      intent(out) :: n    !< number of separated fields
    !--
    character(*), allocatable, intent(out) :: nm(:)  !< separated variables
    integer :: i, j, k
    !----
    
    ! count separator
    n = 0
    i = 1
    j = 1
    do
      i = index(buf(j:), SEP)
      if( i == 0 ) exit
      j = j + i
      n = n + 1
    end do
    n = n + 1  !! #words = #separator + 1
    allocate(nm(n))

    i = 1
    j = 1
    do k=1, n-1
      i = index(buf(j:), SEP)
      nm(k) = trim(adjustl(buf(j:i+j-2)))
      j = j + i
    end do
    nm(n) = trim(adjustl(buf(j:)))

  end subroutine util__split
  !----------------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------------!
  !> Count the line number except for the comment line
  !--
  subroutine util__countline(io, n, comment)

    integer,   intent(in)  :: io
    integer,   intent(out) :: n
    character, intent(in), optional :: comment
    !--
    integer :: stat
    character (256) :: line
    !----

    n = 0
    rewind(io)
    do
       read( io, '(A256)', iostat=stat) line
       if( stat /= 0 ) exit

       if( present( comment ) ) then
          if( line(1:1) == comment ) cycle
       end if

       if( trim(line) /= '' ) then
          n = n + 1
       end if

    end do
    rewind( io )

  end subroutine util__countline
  !----------------------------------------------------------------------------------------------!
  
  !-----------------------------------------------------------------------------------------------!
  subroutine util__readlst(fn, n, lst)

    character(*), intent(in) :: fn  !< filename
    integer, intent(out)     :: n   !< size
    character(*), intent(out), allocatable :: lst(:) !< lst
    !--
    character(256) :: abuf
    integer :: io, ierr
    integer :: i
    !----

    open(newunit=io, file=fn, action='read', status='old', iostat=ierr)
    !if( ierr /= 0 ) error stop "file "//trim(fn) //" not found"
    if( ierr /= 0 ) then
      write(0, '(3a)') "file "//trim(fn) //" not found"
      error stop
    endif
    call util__countline(io, n)

    allocate(lst(n))
    do i=1, n
      read(io,'(A)') abuf
      lst(i) = trim(adjustl(abuf))
    end do

  end subroutine util__readlst
  !-----------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------!
  !> Number of days since 0001-01-01 by Fairfield's formula
  !--
  integer function fairfield(y, m, d) result(nd)
    integer, intent(in) :: y, m, d
    integer :: y0, m0

    if( m==1 .or. m==2 ) then
      y0 = y - 1
      m0 = m + 12
    else  
      y0 = y
      m0 = m
    end if

    nd = 365*y0 + floor(y0/4.) - floor(y/100.) + floor(y/400.) + floor(306*(m0+1)/10.) + d - 428

  end function fairfield
  !-----------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------!
  !> calculate UNIX(POSIX) time from date & time
  !--
  subroutine util__timelocal(yr, mo, dy, hr, mi, sc, tim)

    integer, intent(in) :: yr, mo, dy, hr, mi, sc
    integer, intent(out) :: tim
    integer, parameter :: tim0 = 719163

    tim = (((fairfield(yr, mo, dy) - tim0) * 24 + hr )*60 + mi)*60 + sc

  end subroutine util__timelocal
  !-----------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------!
  !> calculate date & time & number of days in the year from UNIX(POSIX) time
  !--
  subroutine util__localtime(tim, yr, mo, dy, hr, mi, sc, jday)

    integer, intent(in) :: tim
    integer, intent(out) :: yr, mo, dy, hr, mi, sc, jday
    integer :: dy0, sc0
    dy0 = floor( tim / (24*60*60.) ) + 719163
    sc0 = tim - (dy0-719163) * 24*60*60

    do yr=1970, 2099
      if( fairfield(yr+1,1,1) > dy0 ) exit
    end do
    do mo=1, 11
      if( fairfield(yr, mo+1,1) > dy0 ) exit
    end do
    dy = dy0 - fairfield(yr, mo, 1) + 1
    hr = sc0 / 3600
    mi = (sc0 - hr*3600) / 60
    sc = sc0 - hr*3600 - mi * 60
    jday = fairfield(yr, mo, dy) - fairfield(yr, 1, 1) + 1

  end subroutine util__localtime
  !-----------------------------------------------------------------------------------------------!

end module m_util
!-------------------------------------------------------------------------------------------------!
