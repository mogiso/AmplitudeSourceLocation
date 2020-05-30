!!nrtype.F90
!!modules that give precision of floating points
!!Last update: 2019-11-19 11:07:50

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

