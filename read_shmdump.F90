module read_shmdump
  !!read win-text-formatted 1-s waveform (stdout of shmdump)
  !!Author: Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2021
  !!License: MIT License https://opensource.org/licenses/MIT

  private
  public :: read_shmdump_win

  contains

  subroutine read_shmdump_win(chid, conv, yr, mo, dy, hh, mm, ss, waveform, ndata_sec, filtered)
    use nrtype, only : fp
    implicit none

    character(len = 4), intent(in)    :: chid(:)
    real(kind = fp),    intent(in)    :: conv(:)
    integer,            intent(inout) :: yr(:), mo(:), dy(:), hh(:), mm(:), ss(:), ndata_sec(:, :)
    real(kind = fp),    intent(inout) :: waveform(:, :)
    logical,            intent(inout) :: filtered(:, :)
   
    integer,                parameter :: nsample_int_max = 205
    character(len = 4)                :: chid_tmp
    integer                           :: i, j, k, waveform_int(nsample_int_max), sample_int, nsta, nsec_buf, ndata_buf, &
    &                                    yr_tmp, mo_tmp, dy_tmp, hh_tmp, mm_tmp, ss_tmp, nch, ios

    nsta = ubound(chid, 1)
    nsec_buf = ubound(ndata_sec, 1)

    read(*, *, iostat = ios) yr_tmp, mo_tmp, dy_tmp, hh_tmp, mm_tmp, ss_tmp, nch
    if(ios .ne. 0) then
      write(0, *) "cannot read waveform from stdin"
      stop
    endif
    yr(1 : nsec_buf - 1) = yr(2 : nsec_buf); yr(nsec_buf) = yr_tmp
    mo(1 : nsec_buf - 1) = mo(2 : nsec_buf); mo(nsec_buf) = mo_tmp
    dy(1 : nsec_buf - 1) = dy(2 : nsec_buf); dy(nsec_buf) = dy_tmp
    hh(1 : nsec_buf - 1) = hh(2 : nsec_buf); hh(nsec_buf) = hh_tmp
    mm(1 : nsec_buf - 1) = mm(2 : nsec_buf); mm(nsec_buf) = mm_tmp
    ss(1 : nsec_buf - 1) = ss(2 : nsec_buf); ss(nsec_buf) = ss_tmp
  
    ch_loop: do k = 1, nch
      read(*, '(a4, 1x, i3)', advance='no') chid_tmp, sample_int
      read(*, *) (waveform_int(j), j = 1, sample_int)

      station_loop: do j = 1, nsta
        if(trim(chid_tmp) .eq. trim(chid(j))) then
          ndata_buf = 0
          do i = 1, nsec_buf
            if(ndata_sec(i, j) .eq. 0) then
              waveform(ndata_buf + 1 : ndata_buf + sample_int, j) = real(waveform_int(1 : sample_int), kind = fp) * conv(j)
              ndata_sec(i, j) = sample_int
              filtered(i, j) = .false.
              cycle station_loop
            endif
            ndata_buf = ndata_buf + ndata_sec(i, j)
          enddo
          if(ndata_sec(nsec_buf, j) .ne. 0) then
            waveform(1 : ndata_buf - sample_int, j) = waveform(1 + sample_int : ndata_buf, j)
            waveform(ndata_buf - sample_int + 1 : ndata_buf, j) = real(waveform_int(1 : sample_int), kind = fp) * conv(j)
            ndata_sec(1 : nsec_buf - 1, j) = ndata_sec(2 : nsec_buf, j); ndata_sec(nsec_buf, j) = sample_int
            filtered(nsec_buf, j) = .false.
          endif
        endif
      enddo station_loop

    enddo ch_loop

    return
  end subroutine read_shmdump_win

end module read_shmdump
