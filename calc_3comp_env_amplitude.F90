program calc_3comp_envelope
  implicit none

  !!SACの3成分のファイルを読んで、ヘッダはそのまま、(ns ** 2 + ew ** 2 + ud ** 2) / 3のファイルを作る

  integer, parameter :: npts_offset = 500                  !!オフセット取りに使うデータ数。波形先頭からの個数。

  real*8, parameter :: order = 1.0d-6
  real, parameter :: coda_time_begin = 60.0, coda_time_end = 70.0  !!コーダの経過時間の最初と最後

  integer :: iarg, nfile, npts, i, j, count, ptime_index, ios, stime_index, tref_index
  integer :: coda_index_begin, coda_index_end
  integer*4 :: buf_i
  real*8 :: sample, avg_ns, avg_ew, avg_ud
  real*8, allocatable :: data_ns(:), data_ew(:), data_ud(:)
  real :: buf, ptime, begin, end, stime
  real*8 :: amp_p, amp_tref
  character(len = 128) :: infile_ns, infile_ew, infile_ud, outfile, infile_x, infile_y
  character(len = 128), allocatable :: infile(:) 
  character(len = 4) :: header(158)

  !!フィルター関係
  character(len = 8) :: fl_t, fh_t, fs_t
  real*8 :: fl, fh, fs, gn, c
  real*8, allocatable :: h(:)
  integer :: m, n

  !!特性補正
  real*8, allocatable :: deconv_acc(:), deconv_vel(:), deconv_disp(:), wave_past(:), wave_past2(:)
  real*8 :: acc_past
  real*8, parameter :: damp = 0.7d0, natural_freq = 1.0d0
   
  
  iarg = iargc()
  call getarg(1, fl_t); read(fl_t, *) fl
  call getarg(2, fh_t); read(fh_t, *) fh
  call getarg(3, fs_t); read(fs_t, *) fs
  nfile = iarg - 3
  allocate(infile(nfile))
  do i = 1, nfile
    call getarg(i + 3, infile(i))
  enddo

  do i = 1, nfile
    infile_ns = trim(infile(i)) // ".N.sac"
    infile_ew = trim(infile(i)) // ".E.sac"
    infile_ud = trim(infile(i)) // ".U.sac"
    infile_x = trim(infile(i)) // ".X.sac"
    infile_y = trim(infile(i)) // ".Y.sac"

    write(0, *) "infile_ns = ", trim(infile_ns)
    write(0, *) "infile_ew = ", trim(infile_ew)
    write(0, *) "infile_ud = ", trim(infile_ud)
    ios = access(infile_ns, " ")
    if(ios .ne. 0) then
      infile_ns = infile_y
      ios = access(infile_ns, " ")
      if(ios .ne. 0) cycle
    endif
    ios = access(infile_ew, " ")
    if(ios .ne. 0) then
      infile_ew = infile_x
      ios = access(infile_ew, " ")
      if(ios .ne. 0) cycle
    endif

    open(unit = 10, file = infile_ns, form = "unformatted", access = "direct", recl = 4, iostat = ios)
    if(ios .ne. 0) then
      close(10)
      cycle
    endif
    open(unit = 11, file = infile_ew, form = "unformatted", access = "direct", recl = 4, iostat = ios)
    if(ios .ne. 0) then
      close(11)
      cycle
    endif
    open(unit = 12, file = infile_ud, form = "unformatted", access = "direct", recl = 4, iostat = ios)
    if(ios .ne. 0) then
      close(10)
      close(11)
      close(12)
      cycle
    endif

    !!ヘッダ部
    do j = 1, 158
      read(10, rec = j) header(j)
    enddo
    !if(header(123) .ne. "Sman") then
    !  write(0, *) "infile not manual picked"
    !  close(10)
    !  close(11)
    !  close(12)
    !  cycle
    !endif
    !!sample
    read(10, rec = 1) buf; sample = dble(buf)
    !!begin, end
    read(10, rec = 6) begin
    if(abs(buf - 0.01) .gt. 0.0001) then
      write(0, *) "sample time error", buf
      close(10)
      close(11)
      close(12)
      cycle
    endif
    !!npts
    read(10, rec = 80) npts
    write(0, *) "npts = ", npts
    if(npts .lt. 0) then
      close(10)
      close(11)
      close(12)
      cycle
    endif

    end = begin + real(npts) * real(sample)
    if(end .lt. 70.0) then
      write(0, *) "end time error", end
      close(10)
      close(11)
      close(12)
      cycle
    endif

    !!P-time
    !read(10, rec = 6) begin
    read(10, rec = 9, iostat = ios) ptime
    if(ptime .lt. 0.0 .or. ios .ne. 0) then
      close(10)
      close(11)
      close(12)
      cycle
    endif
    read(11, rec = 9, iostat = ios) ptime
    if(ptime .lt. 0.0 .or. ios .ne. 0) then
      close(10)
      close(11)
      close(12)
      cycle
    endif
    read(12, rec = 9, iostat = ios) ptime
    if(ptime .lt. 0.0 .or. ios .ne. 0) then
      close(10)
      close(11)
      close(12)
      cycle
    endif
    read(12, rec = 11) stime
    ptime_index = int((ptime - begin - 5.0) / real(sample) + 0.5)
    stime_index = int((stime - begin - 0.5) / real(sample) + 0.5)
    !tref_index = int((65.0 - begin) / real(sample) + 0.5)
    coda_index_begin = int((coda_time_begin - begin) / real(sample) + 0.5)
    coda_index_end = int((coda_time_end - begin) / real(sample) + 0.5)

    allocate(data_ns(npts), data_ew(npts), data_ud(npts))
    do j = 1, npts
      read(10, rec = 158 + j) buf; data_ns(j) = dble(buf) / order
      read(11, rec = 158 + j) buf; data_ew(j) = dble(buf) / order
      read(12, rec = 158 + j) buf; data_ud(j) = dble(buf) / order
      read(10, rec = 158 + j) buf_i
      if(buf_i .le. -2147483647) exit
      read(11, rec = 158 + j) buf_i
      if(buf_i .le. -2147483647) exit
      read(12, rec = 158 + j) buf_i
      if(buf_i .le. -2147483647) exit
    enddo
    close(10)
    close(11)
    close(12)

    if(buf_i .le. -2147483647) then
      deallocate(data_ns, data_ew, data_ud)
      write(0, *) "buf_i = 0x80000000"
      cycle
    endif

    !!オフセット取り
    avg_ns = 0.0d0
    avg_ew = 0.0d0
    avg_ud = 0.0d0
    do j = 1, npts_offset
      avg_ns = avg_ns + data_ns(j)
      avg_ew = avg_ew + data_ew(j)
      avg_ud = avg_ud + data_ud(j)
    enddo
    data_ns(1 : npts) = data_ns(1 : npts) - avg_ns / dble(npts_offset)
    data_ew(1 : npts) = data_ew(1 : npts) - avg_ew / dble(npts_offset)
    data_ud(1 : npts) = data_ud(1 : npts) - avg_ud / dble(npts_offset)

    !!S/Nチェック用に、P部分の平均を計算
    avg_ns = 0.0d0
    avg_ew = 0.0d0
    avg_ud = 0.0d0
    count = 0
    do j = ptime_index, ptime_index + npts_offset
      avg_ns = avg_ns + data_ns(j) ** 2
      avg_ew = avg_ew + data_ew(j) ** 2
      avg_ud = avg_ud + data_ud(j) ** 2
    enddo
    avg_ns = sqrt(avg_ns)
    avg_ew = sqrt(avg_ew)
    avg_ud = sqrt(avg_ud)
    write(0, '(3(e15.7, 1x))') (avg_ns), (avg_ew), (avg_ud)
    !!波形先頭のmax ampl ratioが一ケタ以上違っていたらエンベロープを計算しない
    if(abs(log10(avg_ns) - log10(avg_ew)) .ge. 1.0d0 .or. abs(log10(avg_ns) - log10(avg_ud)) .ge. 1.0d0 .or. &
    &  abs(log10(avg_ew) - log10(avg_ud)) .ge. 1.0d0) then
      deallocate(data_ns, data_ew, data_ud)
      write(0, *) "amplitude ratio error"
      write(0, '(3(e15.7, 1x))') log10(avg_ns), log10(avg_ew), log10(avg_ud)
      cycle
    endif

#ifdef DECONVOLV
    allocate(deconv_acc(npts), deconv_vel(npts), deconv_disp(npts), wave_past(2), wave_past2(2))
    acc_past = 0.0d0
    wave_past(1 : 2) = 0.0d0
    call deconvolution2(data_ns, npts, damp, natural_freq, sample, deconv_acc, acc_past, wave_past)
    acc_past = 0.0d0
    wave_past(1 : 2) = 0.0d0
    wave_past2(1 : 2) = 0.0d0
    call integral(deconv_acc, npts, sample, acc_past, wave_past, wave_past2, deconv_vel, deconv_disp)
    data_ns(1 : npts) = deconv_vel(1 : npts)
    acc_past = 0.0d0
    wave_past(1 : 2) = 0.0d0
    call deconvolution2(data_ew, npts, damp, natural_freq, sample, deconv_acc, acc_past, wave_past)
    acc_past = 0.0d0
    wave_past(1 : 2) = 0.0d0
    wave_past2(1 : 2) = 0.0d0
    call integral(deconv_acc, npts, sample, acc_past, wave_past, wave_past2, deconv_vel, deconv_disp)
    data_ew(1 : npts) = deconv_vel(1 : npts)
    acc_past = 0.0d0
    wave_past(1 : 2) = 0.0d0
    call deconvolution2(data_ud, npts, damp, natural_freq, sample, deconv_acc, acc_past, wave_past)
    acc_past = 0.0d0
    wave_past(1 : 2) = 0.0d0
    wave_past2(1 : 2) = 0.0d0
    call integral(deconv_acc, npts, sample, acc_past, wave_past, wave_past2, deconv_vel, deconv_disp)
    data_ud(1 : npts) = deconv_vel(1 : npts)
    deallocate(deconv_acc, deconv_vel, deconv_disp, wave_past, wave_past2)
#endif
    


    !!S/Nチェック用に、再度波形先頭部分の平均を計算
    !avg_ns = 0.0d0
    !avg_ew = 0.0d0
    !avg_ud = 0.0d0
    !count = 0
    !do j = 1, npts_offset * 2
    !  if(j .gt. ptime_index) exit
    !  if(abs(data_ns(j)) .gt. avg_ns) avg_ns = data_ns(j)
    !  if(abs(data_ew(j)) .gt. avg_ew) avg_ew = data_ew(j)
    !  if(abs(data_ud(j)) .gt. avg_ud) avg_ud = data_ud(j)
    !enddo
    !write(0, '(3(e15.7, 1x))') (avg_ns), (avg_ew), (avg_ud)
    !!波形先頭のmax ampl ratioが一ケタ以上違っていたらエンベロープを計算しない
    !if(abs(log10(avg_ns) - log10(avg_ew)) .ge. 1.0d0 .or. abs(log10(avg_ns) - log10(avg_ud)) .ge. 1.0d0 .or. &
    !&  abs(log10(avg_ew) - log10(avg_ud)) .ge. 1.0d0) then
    !  deallocate(data_ns, data_ew, data_ud)
    !  write(0, *) "amplitude ratio error"
    !  write(0, '(3(e15.7, 1x))') log10(avg_ns), log10(avg_ew), log10(avg_ud)
    !  cycle
    !endif

    !!S/Nチェック用に、P部分の平均を計算
    !avg_ns = 0.0d0
    !avg_ew = 0.0d0
    !avg_ud = 0.0d0
    !count = 0
    !do j = ptime_index, ptime_index + npts_offset
    !  avg_ns = avg_ns + data_ns(j) ** 2
    !  avg_ew = avg_ew + data_ew(j) ** 2
    !  avg_ud = avg_ud + data_ud(j) ** 2
    !enddo
    !avg_ns = sqrt(avg_ns)
    !avg_ew = sqrt(avg_ew)
    !avg_ud = sqrt(avg_ud)
    !write(0, '(3(e15.7, 1x))') (avg_ns), (avg_ew), (avg_ud)
    !!波形先頭のmax ampl ratioが一ケタ以上違っていたらエンベロープを計算しない
    !if(abs(log10(avg_ns) - log10(avg_ew)) .ge. 1.0d0 .or. abs(log10(avg_ns) - log10(avg_ud)) .ge. 1.0d0 .or. &
    !&  abs(log10(avg_ew) - log10(avg_ud)) .ge. 1.0d0) then
    !  deallocate(data_ns, data_ew, data_ud)
    !  write(0, *) "amplitude ratio error"
    !  write(0, '(3(e15.7, 1x))') log10(avg_ns), log10(avg_ew), log10(avg_ud)
    !  cycle
    !endif

    !!Noise amplitude calculation (before P)
    !amp_p = 0.0d0
    !count = 0
    !do j = ptime_index, ptime_index - 1000, -1
    !  if(j .ge. 1) then
    !    amp_p = amp_p + log(data_ns(j) ** 2 + data_ew(j) ** 2 + data_ud(j) ** 2)
    !    count = count + 1
    !  endif
    !enddo
    !amp_p = sqrt(exp(amp_p / dble(count)))
    !!Tref amplitude calculation
    !amp_tref = 0.0d0
    !count = 0
    !do j = coda_time_begin, coda_time_end
    !  if(j .ge. 1 .and. j .le. npts) then
    !    amp_tref = amp_tref + log(data_ns(j) ** 2 + data_ew(j) ** 2 + data_ud(j) ** 2)
    !    count = count + 1
    !  endif
    !enddo
    !amp_tref = sqrt(exp(amp_tref / dble(count)))

    !if(amp_tref / amp_p .lt. 2.0d0) then
    !  write(0, *) "tref is not larger than amp_p", amp_tref, amp_p
    !  deallocate(data_ns, data_ew, data_ud)
    !  cycle
    !endif

    !data_ns(1 : stime_index - 1) = 0.0d0
    !data_ew(1 : stime_index - 1) = 0.0d0
    !data_ud(1 : stime_index - 1) = 0.0d0

    !!フィルタパラメータ
    call calc_bpf_order(fl, fh, fs, 0.5d0, 5.0d0, sample, m, n, c)
    allocate(h(4 * m))
    call calc_bpf_coef(fl, fh, sample, m, n, h, c, gn)
    !!フィルタ
    call tandem1(data_ns, data_ns, npts, h, m, 1)
    call tandem1(data_ew, data_ew, npts, h, m, 1)
    call tandem1(data_ud, data_ud, npts, h, m, 1)
    data_ns(1 : npts) = data_ns(1 : npts) * gn
    data_ew(1 : npts) = data_ew(1 : npts) * gn
    data_ud(1 : npts) = data_ud(1 : npts) * gn
    !call tandem1(data_ns, data_ns, npts, h, m, -1)
    !call tandem1(data_ew, data_ew, npts, h, m, -1)
    !call tandem1(data_ud, data_ud, npts, h, m, -1)
    !data_ns(1 : npts) = data_ns(1 : npts) * gn
    !data_ew(1 : npts) = data_ew(1 : npts) * gn
    !data_ud(1 : npts) = data_ud(1 : npts) * gn
    deallocate(h)
    !open(unit = 20, file = "log")
    !do j = 1, npts
    !  write(20, '(3(e15.7, 1x))') data_ns(j), data_ew(j), data_ud(j)
    !enddo
    !close(20)



    outfile = trim(infile(i)) // "_" // trim(fl_t) // "-" // trim(fh_t) // "_env.sac"
    header(151) = "ENV "
    write(0, *) "output file = ", trim(outfile)
    open(unit = 13, file = outfile, form = "unformatted", access = "direct", recl = 4)
    do j = 1, 158
      write(13, rec = j) header(j)
    enddo
    do j = 1, npts
      !buf = real((data_ns(j) / order) ** 2 + (data_ew(j) / order) ** 2 + (data_ud(j) / order) ** 2)
      buf = real(((data_ns(j)) ** 2 + (data_ew(j)) ** 2 + (data_ud(j)) ** 2) / 3.0d0)
      !if(j .lt. stime_index) buf = 0.0
      write(13, rec = 158 + j) buf
    enddo
    close(13)

    deallocate(data_ns, data_ew, data_ud)

  enddo

  stop
end program calc_3comp_envelope
    
  
