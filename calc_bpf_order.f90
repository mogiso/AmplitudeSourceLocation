subroutine calc_bpf_order(fl, fh, fs, ap, as, sample, m, n, c)
  use constants, only : pi
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Butterworth band-pass filter by Saito(1978)                
  !!! fl: low frequency cut-off (Hz) not normalized              
  !!! fh: high frequency cut-off (Hz) not normalized             
  !!! fs: stopped frequency (Hz) not normalized                  
  !!! Ap: maximum attenuation in fp (usually 0.5)               
  !!! As: minimum attenuation in fs (usually 5.0)               
  !!! sample: sampling interval (sec)                            
  !!! m: order of filter                                         
  !!! n: order of butterworth function                           
  !!! c: filter coefficient                                      
  !!!                                                            
  !!! Usage:                                                     
  !!! 1. call calc_bpf_order(fl, fh, fs, ap, as, sample, m, n, c)
  !!! 2. allocate filter coefficient allocate(h(4 * m))          
  !!! 3. call calc_bpf_coef(fl, fh, sample, m, n, h, c, gn)      
  !!! 4. call tandem(data, data, ndata, h, m, nml)               
  !!!
  !!! Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!! Reference: Saito, M. (1978) An automatic design algorithm for band selective recursive
  !!!                             digital filters, Geophysical Exploration, 31(4), 240-263 (In Japanese)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8, intent(IN) :: fl,fh, fs, ap, as, sample
  integer, intent(OUT) :: m, n
  real*8, intent(OUT) :: c

  real*8 :: sigma_fl, sigma_fh, sigma_fs, op, os

  sigma_fl = pi * fl * sample
  sigma_fh = pi * fh * sample
  sigma_fs = abs(fs * sample) * pi

  op = sin(sigma_fh - sigma_fl) / (cos(sigma_fl) * cos(sigma_fh))
  os = abs(tan(sigma_fs) - tan(sigma_fh) * tan(sigma_fl) / tan(sigma_fs))

  n = max(2, ifix(abs(real(log(as / ap) / log(op / os))) + 0.5))
  m = n
  c = sqrt(exp(log(ap * as) / dble(n)) / (op * os))

  return
end subroutine calc_bpf_order

