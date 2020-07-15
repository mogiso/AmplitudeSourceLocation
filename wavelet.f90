subroutine kupper_wavelet(ndata, sample, t0, m0, tdiff, wavelet)
  use nrtype,    only : dp, fp
  use constants, only : pi
  implicit none
  !!Kupper wavelet wavelet(t) = 9.0 / 16.0 * pi * m0 / t0 (sin(pi / t0 * t) - 1.0 / 3.0 * sin(3.0 * pi / t0 * t))

  integer, intent(in)         :: ndata            !!number of data for output wavelet
  real(kind = dp), intent(in) :: sample           !!sampling period (s)
  real(kind = dp), intent(in) :: t0, m0           !!t0: characteristic source duration time (s), m0: seismic moment (amplitude)
  real(kind = dp), intent(in) :: tdiff            !!time for wavelet delay (s)
  real(kind = dp), intent(out) :: wavelet(ndata)

  integer :: i, tdiff_int
  real(kind = dp) :: time

  tdiff_int = int(tdiff / sample)
  
  wavelet(1 : ndata) = 0.0_dp
  do i = 1, ndata
    time = sample * dble(i - 1)
    if(time .ge. 0.0_dp .and. time .le. t0) then
      wavelet(i) = 9.0_dp / 16.0_dp * pi * m0 / t0 * &
      &  (sin(pi / t0 * time) - 1.0_dp / 3.0_dp * sin(3.0_dp * pi / t0 * time)) !!From Takemura et al
      !wavelet(i) = m0 * 3.0_dp / 4.0_dp * pi / t0 * sin(pi * time / t0) ** 3  !!From OpenSWPC document
    endif
  enddo

  !!delay the initial time
  do i = ndata, 1, -1
    if(i - tdiff_int .ge. 1) then
      wavelet(i) = wavelet(i - tdiff_int)
    else
      wavelet(i) = 0.0_dp
    endif
  enddo

  return
end subroutine kupper_wavelet

subroutine boxcar_wavelet(ndata, sample, t0, amp, tdiff, wavelet)
  use nrtype,    only : dp, fp
  implicit none
  !!Boxcar function wavelet(t) = 1 / t0 (0<=t<=t0)

  integer, intent(in) :: ndata                     !!number of data for output wavelet
  real(kind = dp), intent(in) :: sample            !!sampling period (s)
  real(kind = dp), intent(in) :: t0, amp           !!t0: duration (s), amp: (amplitude)
  real(kind = dp), intent(in) :: tdiff             !!time for wavelet delay (s)
  real(kind = dp), intent(out) :: wavelet(ndata)

  integer :: i, tdiff_int
  real(kind = dp) :: time

  tdiff_int = int(tdiff / sample)
  
  wavelet(1 : ndata) = 0.0_dp
  do i = 1, ndata
    time = sample * dble(i - 1)
    if(time .ge. 0.0_dp .and. time .le. t0) then
      wavelet(i) = amp / t0
    endif
  enddo

  !!delay the initial time
  do i = ndata, 1, -1
    if(i - tdiff_int .ge. 1) then
      wavelet(i) = wavelet(i - tdiff_int)
    else
      wavelet(i) = 0.0_dp
    endif
  enddo

  return
end subroutine boxcar_wavelet

subroutine triangle_wavelet(ndata, sample, t0, amp, tdiff, wavelet)
  use nrtype,    only : dp, fp
  implicit none
  !!triangle function wavelet(t) = 4 * t / t0**2 (0<=t<=t0/2), -4(t - t0) / t0**2

  integer, intent(in) :: ndata                   !!number of data for output wavelet
  real(kind = dp), intent(in) :: sample            !!sampling period (s)
  real(kind = dp), intent(in) :: t0, amp            !!t0: dulation (s), amp: (amplitude)
  real(kind = dp), intent(in) :: tdiff             !!time for wavelet delay (s)
  real(kind = dp), intent(out) :: wavelet(ndata)

  integer :: i, tdiff_int
  real(kind = dp) :: time

  tdiff_int = int(tdiff / sample)
  
  wavelet(1 : ndata) = 0.0_dp
  do i = 1, ndata
    time = sample * dble(i - 1)
    if(time .ge. 0.0_dp .and. time .le. t0 / 2.0_dp) then
      wavelet(i) = amp * 4.0_dp * time / (t0 ** 2)
    elseif(time .gt. t0 / 2.0_dp .and. time .le. t0) then
      wavelet(i) = amp * (-4.0_dp * (time - t0) / (t0 ** 2))
    endif
  enddo

  !!delay the wavelet
  do i = ndata, 1, -1
    if(i - tdiff_int .ge. 1) then
      wavelet(i) = wavelet(i - tdiff_int)
    else
      wavelet(i) = 0.0_dp
    endif
  enddo

  return
end subroutine triangle_wavelet

subroutine cosine_wavelet(ndata, sample, t0, amp, tdiff, wavelet)
  use nrtype,    only : dp, fp
  use constants, only : pi
  implicit none
  !!cosine function wavelet(t) = 1 / t0 * (1 - cos(2*pi*t/t0))

  integer, intent(in) :: ndata            !!number of data for output wavelet
  real(kind = dp), intent(in) :: sample            !!sampling period (s)
  real(kind = dp), intent(in) :: t0, amp           !!t0: duration (s), amp: (amplitude)
  real(kind = dp), intent(in) :: tdiff             !!time for wavelet delay (s)
  real(kind = dp), intent(out) :: wavelet(ndata)

  integer :: i, tdiff_int
  real(kind = dp) :: time

  tdiff_int = int(tdiff / sample)
  
  wavelet(1 : ndata) = 0.0_dp
  do i = 1, ndata
    time = sample * dble(i - 1)
    if(time .ge. 0.0_dp .and. time .le. t0) then
      wavelet(i) = amp / t0 * (1.0_dp - cos(2.0_dp * pi * time / t0))
    endif
  enddo

  !!delay the wavelet
  do i = ndata, 1, -1
    if(i - tdiff_int .ge. 1) then
      wavelet(i) = wavelet(i - tdiff_int)
    else
      wavelet(i) = 0.0_dp
    endif
  enddo

  return
end subroutine cosine_wavelet


subroutine ricker_wavelet(ndata, sample, fc, amp, tdiff, wavelet)
  use nrtype,    only : dp, fp
  use constants, only : pi
  implicit none

  !!Ricker wavelet (shifted)
  !!wavelet(t) = -(2 * pi ** 2 * fc ** 2 * (t - 1.25 / fc) ** 2 - 1) * exp(-pi ** 2 * fc ** 2 * (t - 1.25 / fc) ** 2)
  !!Meaning of shifted, refer to the GMS wavepage

  integer, intent(in) :: ndata          !!number of data for output wavelet
  real(kind = dp), intent(in) :: sample          !!sampling period (s)
  real(kind = dp), intent(in) :: fc, amp         !!fc: characteristic frequency (inverse of pulse width)
  real(kind = dp), intent(in) :: tdiff             !!time for wavelet delay
  real(kind = dp), intent(out) :: wavelet(ndata)

  integer :: i, tdiff_int
  real(kind = dp) :: time
  wavelet(1 : ndata) = 0.0_dp
  do i = 1, ndata
    time = sample * dble(i - 1)
    if(pi ** 2 * fc ** 2 * (time - 1.25_dp / fc) ** 2 .le. 30.0_dp) then
      wavelet(i) = (2.0_dp * pi ** 2 * fc ** 2 * (time - 1.25_dp / fc) ** 2 - 1.0_dp) &
      &  * exp(-pi ** 2 * fc ** 2 * (time - 1.25_dp / fc) ** 2)
    endif
    wavelet(i) = wavelet(i) * (-amp)
  enddo
  tdiff_int = int(tdiff / sample)

  !!delay the wavelet
  do i = ndata, 1, -1
    if(i - tdiff_int .ge. 1) then
      wavelet(i) = wavelet(i - tdiff_int)
    else
      wavelet(i) = 0.0_dp
    endif
  enddo

  return
end subroutine ricker_wavelet

subroutine smoothed_ramp_function(ndata, sample, fc, amp, tdiff, wavelet)
  use nrtype,    only : dp
  implicit none

  !!Ramp function (shifted) Refer to the GMS wavepage for meaning of shifted
  !!f(t) = 2fc(1 - (tanh(4fc(t-1/fc))) ** 2)

  integer, intent(in) :: ndata          !!number of data for output wavelet
  real(kind = dp), intent(in) :: sample          !!sampling period (s)
  real(kind = dp), intent(in) :: fc, amp         !!fc: characteristic frequency (inverse of pulse width)
  real(kind = dp), intent(in) :: tdiff             !!time for wavelet delay
  real(kind = dp), intent(out) :: wavelet(ndata)

  integer :: i, tdiff_int
  real(kind = dp) :: time

  wavelet(1 : ndata) = 0.0_dp
  do i = 1, ndata
    time = sample * dble(i - 1)
    wavelet(i) = 2.0_dp * fc &
    &  * (1.0_dp &
    &  - ((exp(4.0_dp * fc * (time - 1.0_dp / fc)) - (exp(-4.0_dp * fc * (time - 1.0_dp / fc)))) &
    &  / (exp(4.0_dp * fc * (time - 1.0_dp / fc)) + (exp(-4.0_dp * fc * (time - 1.0_dp / fc))))) ** 2)
    wavelet(i) = wavelet(i) * amp
  enddo

  !!delay the wavelet
  tdiff_int = int(tdiff / sample)
  do i = ndata, 1, -1
    if(i - tdiff_int .ge. 1) then
      wavelet(i) = wavelet(i - tdiff_int)
    else
      wavelet(i) = 0.0_dp
    endif
  enddo
  
  return
end subroutine smoothed_ramp_function

