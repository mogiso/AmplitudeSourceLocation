subroutine tandem1(x, y, n, h, m, nml)
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! nml => 0: Normal filtering
  !!! nml < 0: Reverse filtering
  !!! x(n) : input data
  !!! y(n) : output(filtered) data
  !!!        work well if same name with input data
  !!! h(4 * m) : filter coefficient calculated by calc_[bhl]pf_coef.f90
  !!! m : filter order calclated by calc_[bhl]_pf_order.f90
  !!!
  !!! Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!! Copyright: (c) Masashi Ogiso 2020
  !!! License  : MIT License (https://opensource.org/licenses/MIT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, intent(IN) :: n, m, nml
  real*8, intent(IN) :: x(n), h(4 * m)
  real*8, intent(OUT) :: y(n)
  
  integer :: i, j
  real*8 :: tmp_h(4)


  do j = 1, 4
    tmp_h(j) = h(j)
  enddo

  call recfil1(x, y, n, tmp_h, nml)

  if(m .gt. 1) then
    do i = 2, m
      do j = 1, 4
        tmp_h(j) = h(4 * (i - 1) + j)
      enddo
      call recfil1(y, y, n, tmp_h, nml)
    enddo
  endif

  return
end subroutine tandem1

subroutine recfil1(x, y, n, h, nml)
  implicit none
  integer, intent(IN) :: n, nml
  real*8, intent(IN) :: x(n), h(4)
  real*8, intent(OUT) :: y(n)

  integer :: j, jd, i
  real*8 :: a, aa, b, bb, u1, u2, u3, v1, v2, v3

  if(nml .ge. 0) then
    j = 1
    jd = 1
  else
    j = n
    jd = -1
  endif

  a = h(1)
  aa = h(2)
  b = h(3)
  bb = h(4)
  u1 = 0.0d0
  u2 = 0.0d0
  v1 = 0.0d0
  v2 = 0.0d0

  do i = 1, n
    u3 = u2
    u2 = u1
    u1 = x(j)
    v3 = v2
    v2 = v1
    v1 = u1 + a * u2 + aa * u3 - b * v2 - bb * v3
    y(j) = v1
    j = j + jd
  enddo

  return
end subroutine recfil1

