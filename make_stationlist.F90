program make_stationlist
  use nrtype, only : sp
  implicit none

  integer :: i, nfile
  character(len = 8), allocatable :: stname(:)
  character(len = 129) :: outfile
  character(len = 129), allocatable :: sacfile(:)
  real(kind = sp), allocatable :: stlon(:), stlat(:), stdp(:), ttime_cor(:, :), siteamp(:, :), noiseamp(:)
  logical, allocatable :: use_flag(:)
  

  nfile = iargc()

  allocate(sacfile(1 : nfile))
  do i = 1, nfile
    call getarg(i, sacfile(i))
  enddo

  allocate(stname(1 : nfile), stlon(1 : nfile), stlat(1 : nfile), stdp(1 : nfile), &
  &        use_flag(1 : nfile), ttime_cor(1 : 2, 1 : nfile), siteamp(1 : 2, 1 : nfile), noiseamp(1 : nfile))
  use_flag(1 : nfile) = .true.
  ttime_cor(1 : 2, 1 : nfile) = 0.0_sp
  do i = 1, nfile
    siteamp(1 : 2, i) = [1.0_sp, 0.0_sp]
  enddo
  noiseamp(1 : nfile) = 0.0_sp

  do i = 1, nfile
    open(unit = 10, file = trim(sacfile(i)), form = "unformatted", access = "direct", recl = 4)
    read(10, rec = 32) stlat(i)
    read(10, rec = 33) stlon(i)
    read(10, rec = 35) stdp(i)
    if(stdp(i) .eq. -12345.0_sp) stdp(i) = 0.0_sp
    read(10, rec = 111) stname(i)(1 : 4)
    read(10, rec = 112) stname(i)(5 : 8)
    close(10)
  enddo

  do i = 1, nfile
    print '(2(f0.4, 1x), f5.2, 1x, a8 , 1x, l1, 1x, 5(f3.1, 1x))', &
    &     stlon(i), stlat(i), stdp(i), stname(i), &
    &     use_flag(i), ttime_cor(1 : 2, i), siteamp(1 : 2, i), noiseamp(i)
  enddo

  stop
end program make_stationlist
