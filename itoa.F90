!!itoa.F90
!!convert integer to character
!!max_len should be smaller than 9
!!Last update: 2019-11-21 08:52:01 Add explanations

module itoa

  implicit none
  private

  public :: int_to_char

  contains

  subroutine int_to_char(index, max_len, index_char)
    integer, intent(in) :: index, max_len
    character(len = *), intent(out) :: index_char
    integer :: i

    character(len = 5) :: cfmt

    cfmt = '(ix)'
    write(cfmt(3 : 3), '(i1)') max_len
    write(index_char, cfmt) index
    do i = 1, max_len
      if(index_char(i : i) .eq. " ") index_char(i : i) = "0"
    enddo

    return
  end subroutine int_to_char
end module itoa
    
