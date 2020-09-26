program calc_3comp_envelope
  use nrtype, only : fp, sp
  implicit none

  !!make amplitude list from sac event files for AmplitudeSourceLocation_masterevent.F90
  !!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2020
  !!License  : MIT License (https://opensource.org/licenses/MIT)


  integer, parameter :: npts_offset = 500                  !!number of data used for removal of offset measured from
                                                           !!initial of the waveform
  real(kind = fp), parameter :: order = 1.0_fp
  real(kind = fp), parameter :: amp_timewindow = 10.0_fp
  !integer, parameter :: nsta = 6
  !character(len = 6), parameter :: stname(nsta) = ["V.MEAB", "V.MEAA", "V.NSYM", "V.MNDK", &
  !&                                                "V.KNGM", "V.PNMM"]
  integer, parameter :: nsta = 5
  character(len = 6), parameter :: stname(nsta) = ["V.MEAB", "V.MEAA", "V.PMNS", "V.NSYM", "V.MNDK"]
  !integer, parameter :: nsta = 1
  !character(len = 6), parameter :: stname(nsta) = ["V.MEAB"]

  integer :: iarg, ndate, npts, i, j, k, ios, icount_p, icount_s, ptime_index, stime_index, buf_i
  real(kind = fp) :: sample, avg_ns, avg_ew, avg_ud, mean_amp_p(nsta), mean_amp_s(nsta)
  real(kind = fp), allocatable :: data_ns(:), data_ew(:), data_ud(:)
  real(kind = sp) :: buf, ptime, begin, end, stime
  real(kind = fp) :: amp_p
  character(len = 128) :: infile_ns, infile_ew, infile_ud, outfile, outfile_P, outfile_S
  character(len = 128), allocatable :: infile(:) 
  character(len = 4) :: header(158)
  character(len = 17) :: cfmt 

  !!filter variables
  character(len = 8) :: fl_t, fh_t, fs_t
  real*8 :: fl, fh, fs, gn, c
  real*8, allocatable :: h(:)
  integer :: m, n

  
  iarg = iargc()
  call getarg(1, fl_t); read(fl_t, *) fl
  call getarg(2, fh_t); read(fh_t, *) fh
  call getarg(3, fs_t); read(fs_t, *) fs
  call getarg(4, outfile)
  ndate = iarg - 4
  allocate(infile(ndate))
  do i = 1, ndate
    call getarg(i + 4, infile(i))
  enddo

  outfile_P = trim(outfile) // "_P.txt"
  outfile_S = trim(outfile) // "_S.txt"
  open(30, file = outfile_P)
  open(31, file = outfile_S)
  cfmt = "(i(a, 1x))"
  write(cfmt(2 : 2), '(i1)') nsta + 1
  write(30, trim(cfmt)) "# ", (trim(stname(i)), i = 1, nsta)
  write(31, trim(cfmt)) "# ", (trim(stname(i)), i = 1, nsta)

  cfmt = "(i(e15.7, 1x), a)"
  write(cfmt(2 : 2), '(i1)') nsta

  date_loop: do j = 1, ndate
    station_loop: do i = 1, nsta
      !infile_ns = trim(infile(j)) // "." // trim(stname(i)) // ".N.sac"
      !infile_ew = trim(infile(j)) // "." // trim(stname(i)) // ".E.sac"
      infile_ud = trim(infile(j)) // "." // trim(stname(i)) // ".U.sac"

      !write(0, *) "infile_ns = ", trim(infile_ns)
      !write(0, *) "infile_ew = ", trim(infile_ew)
      write(0, *) "infile_ud = ", trim(infile_ud)

      !open(unit = 10, file = infile_ns, form = "unformatted", access = "direct", recl = 4, iostat = ios)
      !if(ios .ne. 0) then
      !  close(10)
      !  cycle date_loop
      !endif
      !open(unit = 11, file = infile_ew, form = "unformatted", access = "direct", recl = 4, iostat = ios)
      !if(ios .ne. 0) then
      !  close(10)
      !  close(11)
      !  cycle date_loop
      !endif
      open(unit = 12, file = infile_ud, form = "unformatted", access = "direct", recl = 4, iostat = ios)
      if(ios .ne. 0) then
        !close(10)
        !close(11)
        close(12)
        cycle date_loop
      endif

      !!read header
      !!sample
      read(12, rec = 1) buf; sample = dble(buf)
      !!begin, end
      read(12, rec = 6) begin
      if(abs(buf - 0.01) .gt. 0.0001) then
        write(0, *) "sample time error", buf
        !close(10)
        !close(11)
        close(12)
        cycle date_loop
      endif
      !!npts
      read(12, rec = 80) npts
      write(0, *) "npts = ", npts
      if(npts .lt. 0) then
        !close(10)
        !close(11)
        close(12)
        cycle date_loop
      endif

      !!P-time
      !read(10, rec = 6) begin
      !read(10, rec = 9, iostat = ios) ptime
      !if(ptime .lt. 0.0 .or. ios .ne. 0) then
        !close(10)
        !close(11)
        !close(12)
        !cycle date_loop
      !endif
      !read(11, rec = 9, iostat = ios) ptime
      !if(ptime .lt. 0.0 .or. ios .ne. 0) then
        !close(10)
        !close(11)
        !close(12)
        !cycle date_loop
      !endif
      read(12, rec = 9, iostat = ios) ptime
      if(ptime .lt. 0.0 .or. ios .ne. 0) then
        !close(10)
        !close(11)
        close(12)
        cycle date_loop
      endif
      read(12, rec = 11) stime
      if(ptime .eq. -12345.0 .or. stime .eq. -12345.0) then
        !close(10)
        !close(11)
        close(12)
        cycle date_loop
      endif
      ptime_index = int((ptime - begin) / real(sample) + 0.5) + 1
      stime_index = int((stime - begin) / real(sample) + 0.5) + 1
      write(0, '(a, f0.2)') "begin = ", begin
      write(0, '(a, 2(i0, 1x))') "(P|S)time index = ", ptime_index, stime_index

      allocate(data_ns(npts), data_ew(npts), data_ud(npts))
      do k = 1, npts
        !read(10, rec = 158 + k) buf; data_ns(k) = dble(buf) / order
        !read(11, rec = 158 + k) buf; data_ew(k) = dble(buf) / order
        read(12, rec = 158 + k) buf; data_ud(k) = dble(buf) / order
        !read(10, rec = 158 + k) buf_i
        !if(buf_i .le. -2147483647) then
          !close(10)
          !close(11)
          !close(12)
          !deallocate(data_ns, data_ew, data_ud)
          !cycle date_loop
        !endif
        !read(11, rec = 158 + k) buf_i
        !if(buf_i .le. -2147483647) then
          !close(10)
          !close(11)
          !close(12)
          !deallocate(data_ns, data_ew, data_ud)
          !cycle date_loop
        !endif
        read(12, rec = 158 + k) buf_i
        if(buf_i .le. -2147483647) then
          !close(10)
          !close(11)
          close(12)
          deallocate(data_ns, data_ew, data_ud)
          cycle date_loop
        endif
      enddo
      !close(10)
      !close(11)
      close(12)

      !!remove offset
      !avg_ns = 0.0d0
      !avg_ew = 0.0d0
      avg_ud = 0.0d0
      do k = 1, npts_offset
        !avg_ns = avg_ns + data_ns(k)
        !avg_ew = avg_ew + data_ew(k)
        avg_ud = avg_ud + data_ud(k)
      enddo
      !data_ns(1 : npts) = data_ns(1 : npts) - avg_ns / dble(npts_offset)
      !data_ew(1 : npts) = data_ew(1 : npts) - avg_ew / dble(npts_offset)
      data_ud(1 : npts) = data_ud(1 : npts) - avg_ud / dble(npts_offset)

      !!フィルタパラメータ
      call calc_bpf_order(fl, fh, fs, 0.5d0, 5.0d0, sample, m, n, c)
      allocate(h(4 * m))
      call calc_bpf_coef(fl, fh, sample, m, n, h, c, gn)
      !!フィルタ
      !call tandem1(data_ns, data_ns, npts, h, m, 1)
      !call tandem1(data_ew, data_ew, npts, h, m, 1)
      call tandem1(data_ud, data_ud, npts, h, m, 1)
      !data_ns(1 : npts) = data_ns(1 : npts) * gn
      !data_ew(1 : npts) = data_ew(1 : npts) * gn
      data_ud(1 : npts) = data_ud(1 : npts) * gn
      deallocate(h)

      mean_amp_p(i) = 0.0_fp
      mean_amp_s(i) = 0.0_fp
      icount_p = 0
      icount_s = 0
      do k = 1, int(amp_timewindow / sample + 0.5)
        !mean_amp_p(i) = mean_amp_p(i) + sqrt((data_ns(ptime_index + k - 1) ** 2  &
        !&                                   + data_ew(ptime_index + k - 1) ** 2  &
        !&                                   + data_ud(ptime_index + k - 1) ** 2) &
        !&                                   / 3.0_fp)
        mean_amp_p(i) = mean_amp_p(i) + sqrt(data_ud(ptime_index + k - 1) ** 2)
        icount_p = icount_p + 1
        !mean_amp_s(i) = mean_amp_s(i) + sqrt((data_ns(stime_index + k - 1) ** 2  &
        !&                                   + data_ew(stime_index + k - 1) ** 2  &
        !&                                   + data_ud(stime_index + k - 1) ** 2) &
        !&                                   / 3.0_fp)
        mean_amp_s(i) = mean_amp_s(i) + sqrt(data_ud(stime_index + k - 1) ** 2)
        icount_s = icount_s + 1
      enddo
      mean_amp_p(i) = mean_amp_p(i) / real(icount_p, kind = fp)
      mean_amp_s(i) = mean_amp_s(i) / real(icount_s, kind = fp)

      deallocate(data_ns, data_ew, data_ud)
    enddo station_loop
    write(30, trim(cfmt)) (mean_amp_p(i), i = 1, nsta), trim(infile(j))
    write(31, trim(cfmt)) (mean_amp_s(i), i = 1, nsta), trim(infile(j))

  enddo date_loop
  close(30)
  close(31)

  stop
end program calc_3comp_envelope
    
  
