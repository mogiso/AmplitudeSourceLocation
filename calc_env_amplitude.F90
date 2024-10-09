program calc_env_amplitude
  use nrtype, only : fp, sp
  implicit none

  !!make amplitude list from sac event files for AmplitudeSourceLocation_masterevent.F90
  !!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2020
  !!License  : MIT License (https://opensource.org/licenses/MIT)


  integer, parameter :: npts_offset = 500                  !!number of data used for removal of offset measured from
                                                           !!initial of the waveform
  real(kind = fp), parameter :: order = 1.0_fp

  integer :: iarg, nsta, ndate, npts, i, j, k, ios, icount_p, icount_s, ptime_index, stime_index, buf_i
  real(kind = fp) :: sample, avg_offset, stlon, stlat, stdp
  real(kind = fp), allocatable :: wavedata(:), mean_amp_p(:), mean_amp_s(:)
  real(kind = sp) :: buf, ptime, begin, stime
  real(kind = fp) :: amp_p, amp_timewindow
  character(len = 128) :: stationlist, infile_sac, outfile_p, outfile_s
  character(len = 128), allocatable :: infile(:) 
  character(len = 6),   allocatable :: stname(:)
  character(len = 4) :: header(158), nsta_c
  character(len = 20) :: cfmt, amp_timewindow_t, cmpnm

  !!filter variables
  character(len = 8) :: fl_t, fh_t, fs_t
  real*8 :: fl, fh, fs, gn, c
  real*8, allocatable :: h(:)
  integer :: m, n

  
  iarg = iargc()
  if(iarg .lt. 9) then
    write(0, '(a)', advance = "no") "usage: ./calc_env_amplitude (station_param_file) (component_name) (fl) (fh) (fs) "
    write(0, '(a)')              "(rms_time_window_length) (outfile_txt_P) (outfile_txt_S) (sacfile_index1) (sacfile_index2) ..."
    error stop
  endif
  call getarg(1, stationlist)
  call getarg(2, cmpnm)
  call getarg(3, fl_t); read(fl_t, *) fl
  call getarg(4, fh_t); read(fh_t, *) fh
  call getarg(5, fs_t); read(fs_t, *) fs
  call getarg(6, amp_timewindow_t); read(amp_timewindow_t, *) amp_timewindow
  call getarg(7, outfile_p)
  call getarg(8, outfile_s)
  ndate = iarg - 8
  allocate(infile(ndate))
  do i = 1, ndate
    call getarg(i + 8, infile(i))
  enddo

  open(unit = 10, file = stationlist)
  nsta = 0
  do
    read(10, *, iostat = ios)
    if(ios .ne. 0) exit
    nsta = nsta + 1
  enddo
  allocate(stname(nsta), mean_amp_p(nsta), mean_amp_s(nsta))
  rewind(10)
  do i = 1, nsta
    read(10, *) stlon, stlat, stdp, stname(i)
  enddo
  close(10)
  print *, "nsta = ", nsta

  write(nsta_c, '(i0)') nsta
  cfmt = "(a, " // trim(nsta_c) // "(a, 1x))"

  open(30, file = outfile_p)
  open(31, file = outfile_s)
  write(30, trim(cfmt)) "# ", (trim(stname(i)), i = 1, nsta)
  write(31, trim(cfmt)) "# ", (trim(stname(i)), i = 1, nsta)

  cfmt = "(" // trim(nsta_c) // "(e15.7, 1x), a)"

  date_loop: do j = 1, ndate
    station_loop: do i = 1, nsta
      print *, i
#ifdef TESTDATA
      infile_sac = trim(infile(j)) // "." // trim(stname(i)) // ".env_cal.sac"
#else
      infile_sac = trim(infile(j)) // "." // trim(stname(i)) // "." // trim(cmpnm) // ".sac"
#endif

      write(0, *) "infile_sac = ", trim(infile_sac)

      open(unit = 12, file = infile_sac, form = "unformatted", access = "direct", recl = 4, iostat = ios)
      if(ios .ne. 0) then
        close(12)
        cycle date_loop
      endif

      !!read header
      !!sample
      read(12, rec = 1) buf; sample = dble(buf)
      !!begin, end
      read(12, rec = 6) begin
      !!npts
      read(12, rec = 80) npts
      write(0, *) "npts = ", npts
      if(npts .lt. 0) then
        close(12)
        cycle date_loop
      endif

      read(12, rec = 11, iostat = ios) stime
      if(stime .eq. -12345.0) then
        close(12)
        cycle date_loop
      endif
      read(12, rec = 9, iostat = ios) ptime
      if(ptime .eq. -12345.0) then
        ptime = stime
      endif
      ptime_index = int((ptime - begin) / real(sample) + 0.5) + 1
      stime_index = int((stime - begin) / real(sample) + 0.5) + 1
      write(0, '(a, f0.2)') "begin = ", begin
      write(0, '(a, 2(i0, 1x))') "(P|S)time index = ", ptime_index, stime_index

      allocate(wavedata(npts))
      do k = 1, npts
        read(12, rec = 158 + k) buf; wavedata(k) = dble(buf) / order
        read(12, rec = 158 + k) buf_i
        if(buf_i .le. -2147483647) then
          close(12)
          deallocate(wavedata)
          cycle date_loop
        endif
      enddo
      close(12)

#ifndef TESTDATA
      avg_offset = 0.0d0
      do k = 1, npts_offset
        avg_offset = avg_offset + wavedata(k)
      enddo
      wavedata(1 : npts) = wavedata(1 : npts) - avg_offset / dble(npts_offset)

      !!calculate coefficients of band-pass filter
      call calc_bpf_order(fl, fh, fs, 0.5d0, 5.0d0, sample, m, n, c)
      allocate(h(4 * m))
      call calc_bpf_coef(fl, fh, sample, m, n, h, c, gn)
      call tandem1(wavedata, wavedata, npts, h, m, 1)
      wavedata(1 : npts) = wavedata(1 : npts) * gn
      deallocate(h)
#endif

      mean_amp_p(i) = 0.0_fp
      mean_amp_s(i) = 0.0_fp
      icount_p = 0
      icount_s = 0
      do k = 1, int(amp_timewindow / sample + 0.5_fp)
        mean_amp_p(i) = mean_amp_p(i) + wavedata(ptime_index + k - 1) ** 2
        icount_p = icount_p + 1
        mean_amp_s(i) = mean_amp_s(i) + wavedata(stime_index + k - 1) ** 2
        icount_s = icount_s + 1
      enddo
      mean_amp_p(i) = sqrt(mean_amp_p(i) / real(icount_p, kind = fp))
      mean_amp_s(i) = sqrt(mean_amp_s(i) / real(icount_s, kind = fp))

      deallocate(wavedata)
    enddo station_loop
    write(30, trim(cfmt)) (mean_amp_p(i), i = 1, nsta), trim(infile(j))
    write(31, trim(cfmt)) (mean_amp_s(i), i = 1, nsta), trim(infile(j))

  enddo date_loop
  close(30)
  close(31)

  stop
end program calc_env_amplitude
    
  
