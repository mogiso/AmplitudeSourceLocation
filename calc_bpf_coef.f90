subroutine calc_bpf_coef(fl, fh, sample, m, n, h, c, gn)
  use constants, only : pi
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Butterworth band-pass filter by Saito(1978)
  !!! Input
  !!! fl, fh: low- and high- passband frequency (not normalized)
  !!! sample: sampling sequence (sec)
  !!! m, n: filter order (calculated by calc_bpf_order.f90)
  !!! h(4 * m) : filter coefficient
  !!! gn : gain factor
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, intent(IN) :: m, n
  real*8, intent(IN) :: c, fl, fh, sample
  real*8, intent(OUT) :: h(4 * m)
  real*8, intent(OUT) :: gn

  real*8 :: g, fj, a1, dp, sigma_fl, sigma_fh, re, ri, wpc, wmc, clh
  complex*16 :: oj, r(2), cq
  integer :: i, j, k

  dp = pi / 2.0d0 / dble(n)
  sigma_fl = pi * fl * sample
  sigma_fh = pi * fh * sample

  clh = 1.0d0 / (cos(sigma_fl) * cos(sigma_fh))
  g = 1.0d0
  fj = 1.0d0
  k = 0

  do j = 1, int(n / 2)
    oj = dcmplx(cos(dp * fj), sin(dp * fj)) * 0.5d0
    fj = fj + 2.0d0
    cq = sqrt(oj ** 2 + c ** 2 * tan(sigma_fl) * tan(sigma_fh))
    r(1) = oj + cq
    r(2) = oj - cq
    g = g * c ** 2
    do i = 1, 2
      re = dble(r(i)) ** 2
      ri = imag(r(i))
      a1 = 1.0d0 / ((c + ri) ** 2 + re)
      g = g * a1
      h(k + 1) = 0.0d0
      h(k + 2) = -1.0d0
      h(k + 3) = 2.0d0 * ((ri - c) * (ri + c) + re) * a1
      h(k + 4) = ((ri - c) ** 2 + re) * a1
      k = k + 4
    enddo
  enddo
  gn = g
  if(n .ne. 2 * int(n / 2)) then
    wpc = c ** 2 * cos(sigma_fh - sigma_fl) * clh
    wmc = - c ** 2 * cos(sigma_fh + sigma_fl) * clh
    a1 = 1.0d0 / (wpc + c)
    gn = g * c * a1
    h(k + 1) = 0.0d0
    h(k + 2) = -1.0d0
    h(k + 3) = 2.0d0 * wmc * a1
    h(k + 4) = (wpc - c) * a1
  endif

  return
end subroutine calc_bpf_coef

