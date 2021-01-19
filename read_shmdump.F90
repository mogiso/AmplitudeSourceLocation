module read_shmdump
  !!read win-text-formatted 1-s waveform (stdout of shmdump)
  !!Author: Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2021
  !!License: MIT License https://opensource.org/licenses/MIT

  private
  public :: read_shmdump_win

  contains

  subroutine read_shmdump_win(chid_order, conv, yr, mo, dy, hh, mm, ss, waveform, ndata_sec)
    use nrtype, only : fp
    implicit none

    integer           , intent(in)    :: chid_order(:)
    real(kind = fp),    intent(in)    :: conv(:)
    integer,            intent(inout) :: yr(:), mo(:), dy(:), hh(:), mm(:), ss(:), ndata_sec(:, :)
    real(kind = fp),    intent(inout) :: waveform(:, :)
   
    integer,                parameter :: nsample_int_max = 205
    character(len = 4)                :: chid_tmp
    integer                           :: i, j, waveform_int(nsample_int_max), sample_int, nsta, nbuf_sec, ndata_tmp, &
    &                                    yr_tmp, mo_tmp, dy_tmp, hh_tmp, mm_tmp, ss_tmp, nch, ios

    nsta = ubound(chid_order, 1)
    nbuf_sec = ubound(ndata_sec, 1)

    read(6, *, iostat = ios) yr_tmp, mo_tmp, dy_tmp, hh_tmp, mm_tmp, ss_tmp, nch
    if(ios .ne. 0) then
      write(0, '(a)') "cannot read waveforms from stdin"
      error stop
    endif
    yr(1 : nbuf_sec - 1) = yr(2 : nbuf_sec); yr(nbuf_sec) = yr_tmp
    mo(1 : nbuf_sec - 1) = mo(2 : nbuf_sec); yr(nbuf_sec) = mo_tmp
    dy(1 : nbuf_sec - 1) = dy(2 : nbuf_sec); yr(nbuf_sec) = dy_tmp
    hh(1 : nbuf_sec - 1) = hh(2 : nbuf_sec); yr(nbuf_sec) = hh_tmp
    mm(1 : nbuf_sec - 1) = mm(2 : nbuf_sec); yr(nbuf_sec) = mm_tmp
    ss(1 : nbuf_sec - 1) = ss(2 : nbuf_sec); yr(nbuf_sec) = ss_tmp
  
    do j = 1, nch
      read(6, '(a, 1x, i0)', advance="no") chid_tmp, sample_int
      read(6, *) (waveform_int(i), i = 1, sample_int)
      ndata_tmp = sum(ndata_sec(1 : nbuf_sec, chid_order(j)))
      waveform(1 : ndata_tmp - sample_int, chid_order(j)) = waveform(1 + sample_int : ndata_tmp, chid_order(j))
      waveform(ndata_tmp + 1 : ndata_tmp + sample_int, chid_order(j)) = waveform_int(1 : sample_int) * conv(chid_order(j))
    enddo

    return
  end subroutine read_shmdump_win

end module read_shmdump
