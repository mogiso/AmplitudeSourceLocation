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
    integer                           :: waveform_int(nsample_int_max), sample_int, nsta 

    nsta = ubound(chid_order, 1)
    nbufsec = ubound(ndata_sec, 1)

    read(6, *) yr_tmp, mo_tmp, dy_tmp, hh_tmp, mm_tmp, ss_tmp, nch
    yr(1 : nbufsec - 1) = yr(2 : nbufsec); yr(nbufsec) = yr_tmp
    mo(1 : nbufsec - 1) = mo(2 : nbufsec); yr(nbufsec) = yr_tmp
    dy(1 : nbufsec - 1) = dy(2 : nbufsec); yr(nbufsec) = yr_tmp
    hh(1 : nbufsec - 1) = hh(2 : nbufsec); yr(nbufsec) = yr_tmp
    mm(1 : nbufsec - 1) = mm(2 : nbufsec); yr(nbufsec) = yr_tmp
    ss(1 : nbufsec - 1) = ss(2 : nbufsec); yr(nbufsec) = yr_tmp
  
    do j = 1, nch
      read(6, *, advance = "no") chid_tmp, sample_int
      read(6, *) (waveform_int(i), i = 1, sample_int)
      ndata_tmp = sum(ndata_sec(1 : nbuf_sec, chid_order(j)))
      waveform(1 : ndata_tmp - sample_int, chid_order(j)) = waveform(1 + sample_int : ndata_tmp, chid_order(j))
      waveform(ndata_tmp + 1 : ndata_tmp + sample_int, chid_order(j)) = waveform_int(1 : sample_int) * conv(chid_order(j))
    enddo

    return
  end subroutine read_shmdump_win

end module read_shmdump
