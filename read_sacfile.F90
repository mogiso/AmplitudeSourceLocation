module read_sacfile
  private
  public :: read_sachdr, read_sacdata

  contains

  subroutine read_sachdr(infile, &
  &  header, delta, begin, origin, ptime, t0, t1, t2, t3, stlat, stlon, stdp, eqlat, eqlon, eqdep, npts, &
  &  user7, user8, user9, stname)
    use nrtype, only : fp, sp
    implicit none
    character(len = *), intent(in) :: infile
    character(len = 4), intent(out), optional :: header(158)
    real(kind = fp), intent(out), optional :: delta
    real(kind = fp), intent(out), optional :: begin
    real(kind = fp), intent(out), optional :: origin
    real(kind = fp), intent(out), optional :: ptime
    real(kind = fp), intent(out), optional :: t0
    real(kind = fp), intent(out), optional :: t1
    real(kind = fp), intent(out), optional :: t2
    real(kind = fp), intent(out), optional :: t3
    real(kind = fp), intent(out), optional :: stlat
    real(kind = fp), intent(out), optional :: stlon
    real(kind = fp), intent(out), optional :: stdp
    real(kind = fp), intent(out), optional :: eqlat
    real(kind = fp), intent(out), optional :: eqlon
    real(kind = fp), intent(out), optional :: eqdep
    integer,         intent(out), optional :: npts
    real(kind = fp), intent(out), optional :: user7
    real(kind = fp), intent(out), optional :: user8
    real(kind = fp), intent(out), optional :: user9
    character(len = 8), intent(out), optional :: stname
  
    real(kind = sp) :: buf
    integer :: i
    character(len = 4) :: stname1, stname2
  
    open(unit = 10, file = infile, form = "unformatted", access = "direct", recl = 4)
    print '(2a)', "reading sacfile, ", trim(infile)
    if(present(header)) then
      do i = 1, 158
        read(10, rec = i) header(i)
      enddo
    endif
    if(present(delta)) then
      read(10, rec = 1) buf; delta = real(buf, kind = fp)
    endif
    if(present(begin)) then
       read(10, rec = 6) buf; begin = real(buf, kind = fp)
    endif
    if(present(origin)) then
      read(10, rec = 8) buf; origin = real(buf, kind = fp)
    endif
    if(present(ptime)) then
       read(10, rec = 9) buf; ptime = real(buf, kind = fp)
    endif
    if(present(t0)) then
      read(10, rec = 11) buf; t0 = real(buf, kind = fp)
    endif
    if(present(t1)) then
      read(10, rec = 12) buf; t1 = real(buf, kind = fp)
    endif
    if(present(t2)) then
      read(10, rec = 13) buf; t2 = real(buf, kind = fp)
    endif
    if(present(t3)) then
      read(10, rec = 14) buf; t3 = real(buf, kind = fp)
    endif
    if(present(stlat)) then
      read(10, rec = 32) buf; stlat = real(buf, kind = fp)
    endif
    if(present(stlon)) then
      read(10, rec = 33) buf; stlon = real(buf, kind = fp)
    endif
    if(present(stdp)) then
      read(10, rec = 35) buf; stdp = real(buf, kind = fp)
    endif
    if(present(eqlat)) then
      read(10, rec = 36) buf; eqlat = real(buf, kind = fp)
    endif
    if(present(eqlon)) then
      read(10, rec = 37) buf; eqlon = real(buf, kind = fp)
    endif
    if(present(eqdep)) then
      read(10, rec = 39) buf; eqdep = real(buf, kind = fp)
    endif
    if(present(user7)) then
      read(10, rec = 48) buf; user7 = real(buf, kind = fp)
    endif
    if(present(user8)) then
      read(10, rec = 49) buf; user8 = real(buf, kind = fp)
    endif
    if(present(user9)) then
      read(10, rec = 50) buf; user9 = real(buf, kind = fp)
    endif
    if(present(npts)) then
      read(10, rec = 80) npts
    endif
    if(present(stname)) then
      read(10, rec = 111) stname1
      read(10, rec = 112) stname2
      stname = stname1 // stname2
    endif
    close(10)
    return
  end subroutine read_sachdr
  
  subroutine read_sacdata(infile, npts, wavedata)
    use nrtype, only : sp, dp
    implicit none
    character(len = *), intent(in)  :: infile
    integer,              intent(in)  :: npts
    real(kind = dp),      intent(out) :: wavedata(npts)
    integer         :: i, ios
    real(kind = sp) :: buf

    wavedata(1 : npts) = 0.0_dp
    open(unit = 10, file = infile, form = "unformatted", access = "direct", recl = 4)
    do i = 1, npts
      read(10, rec = 158 + i, iostat = ios) buf
      if(ios .ne. 0) exit
      wavedata(i) = real(buf, kind = dp)
    enddo
    close(10)
    return
  end subroutine read_sacdata
end module read_sacfile 
  
