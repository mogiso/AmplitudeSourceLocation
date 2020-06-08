!!nrtype.F90
!!modules that give precision of floating points
!!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
!!Copyright: (c) Masashi Ogiso 2020
!!License  : MIT License https://opensource.org/licenses/MIT

module nrtype
  implicit none
  private

  integer, public, parameter :: sp = selected_real_kind(6)
  integer, public, parameter :: dp = selected_real_kind(15)

#ifdef DOUBLE
  integer, public, parameter :: fp = dp
#else
  integer, public, parameter :: fp = sp
#endif

end module nrtype

