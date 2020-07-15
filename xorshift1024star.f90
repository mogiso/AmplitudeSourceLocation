! Written in 2015-2016 by Shun Sakuraba
! 
! To the extent possible under law, the author has dedicated all copyright
! and related and neighboring rights to this software to the public domain
! worldwide. This software is distributed without any warranty.
! See <http://creativecommons.org/publicdomain/zero/1.0/>.
! 
! Vigna's xorshift1024* pseudorandom generator.
! (Sebastiano Vigna. An experimental exploration of Marsaglia's xorshift generators, scrambled. CoRR, abs/1402.6246, 2014.)
! xorshift1024* is a pseudorandom generator with a reasonable speed and a good state space size. This is a standard choice of the generator.
!
! This file is mostly based on F95 standard despite extentioned f90.
! This is because .f95 is unrecognized on ifort 16.
!
! XXX: Because Fortran standard does not specify integer model,
! int64 overflow is undefined behaviour.
! Thus, there might be warnings issued while compiling this program. (Or worse, it may not work)
! A workaround is using int(int(a, kind=16) "op" b, kind=8) trick, which works well with gfortran,
! but commercial compilers do not accept 128-bit arithmetics (incl. Intel, Cray, and Fujitsu).

module xorshift1024star
  implicit none
  
  ! random number state
  public :: xorshift1024star_state
  ! state initialization functions
  public :: state_init_full, state_init, state_init_self
  ! global state initialization functions
  public :: rand_init, rand_init_self

  ! jump state by 2^512 (overall state space is 2^1024 - 1)
  public :: state_jump
  ! jump global state
  public :: rand_jump

  ! Draw integer from range, or uniform random number
  public :: draw_integer, draw_uniform
  ! Draw integer from range, or uniform random number, from global rng state
  public :: rand_integer, rand_uniform

  type xorshift1024star_state
     integer(8) :: s(0:15)
     integer :: ptr
  end type xorshift1024star_state

  type(xorshift1024star_state) :: global_state
  ! Generated from splitmix64(0)
  data global_state%s &
       /-2152535657050944081_8, &! 0xE220A8397B1DCDAF
        -2204311364474961569_8, &! 0xE168B67E321D015F
         -564803942837962399_8, &! 0xF82969CE73673961
        -6927212692750830023_8, &! 0x9FDD98C205027639
         7587201670112224958_8, &! 0x694B27B8624A26BE
        -7737210624469382739_8, &! 0x949FE7EE381FC9AD
         8132190365596779910_8, &! 0x70DB580D23400586
         5952206671755968406_8, &! 0x529A7C6E8BB97F96
         -243991858305861371_8, &! 0xFC9D2AB295977905
         6084237800721128618_8, &! 0x546F8E0B48FC9CAA
        -9131961608637716695_8, &! 0x8144C14CD6A4E729
        -4970413959199997670_8, &! 0xBB058AD6A73DC11A
        -7880794127290781191_8, &! 0x92A1CB7ED6D931F9
         9150292177230689029_8, &! 0x7EFC5E413EDEF305
        -5669289616600066314_8, &! 0xB152A32181BE16F6
          214460927578729249_8/  ! 0x02F9EB13CE87EB21
  data global_state%ptr /0/
contains

  subroutine state_init_full(state, s)
    implicit none
    type(xorshift1024star_state), intent(out) :: state
    integer(8), intent(in) :: s(16)
    integer :: i
    integer(8) :: seed

    if(all(s(:) == 0_8)) then
       ! initialize with some non-zero values
       seed = 0
       do i = 0, 15
          state%s(i) = splitmix64(seed)
       end do
    else
       state%s(0:15) = s(1:16)
    endif
    state%ptr = 0
  end subroutine state_init_full

  subroutine state_jump(state)
    implicit none
    type(xorshift1024star_state), intent(out) :: state
    integer(8) :: jump_table(16)
    integer(8) :: t(0:15)
    integer(8) :: tmp
    integer :: i, j, b
    data jump_table / -8924956236279331811_8,& !0x84242F96ECA9C41D
                      -6645523562763818923_8,& !0xA3C65B8776F96855
                       6572057659653707831_8,& !0x5B34A39F070B5837
                       4938671967096216094_8,& !0x4489AFFCE4F31A1E
                       3458459993260519232_8,& !0x2FFEEB0A48316F40
                      -2581239258607468510_8,& !0xDC2D9891FE68C022
                       3916182429352585840_8,& !0x3659132BB12FEA70
                      -6142490363719071048_8,& !0xAAC17D8EFA43CAB8
                      -4266174017505289453_8,& !0xC4CB815590989B13
                       6839126324828817723_8,& !0x5EE975283D71C93B
                       7572038374137779520_8,& !0x691548C86C1BD540
                       8723688107328792229_8,& !0x7910C41D10A1E6A5
                        819591658532496040_8,& !0x0B5FC64563B3E2A8
                        324108011427370141_8,& !0x047F7684E9FC949D
                      -5075132425047734838_8,& !0xB99181F2D8F685CA
                       2902007988922235075_8/  !0x284600E3F30E38C3
    
    t(:) = 0

    do i = 1, 16
       do b = 0, 63
          if(btest(jump_table(i), b)) then
             do j = 0, 15
                t(j) = ieor(t(j), state%s(iand(j + state%ptr, 15)))
             end do
             tmp = draw_integer8(state)
          end if
       end do
    end do
    state%s = t
  end subroutine state_jump

  subroutine rand_jump()
    implicit none
    call state_jump(global_state)
  end subroutine rand_jump

  integer(8) function draw_integer8(state)
    implicit none
    type(xorshift1024star_state), intent(inout) :: state
    
    integer(8) :: s0, s1, r
    integer(8), parameter :: spreader = 1181783497276652981_8
    
    ! get state variables
    s0 = state%s(state%ptr)
    state%ptr = iand(state%ptr + 1, 15)
    s1 = state%s(state%ptr)

    ! xorshift
    s1 = ieor(s1, ishft(s1, 31))
    s1 = ieor(s1, ishft(s1, -11))
    s0 = ieor(s0, ishft(s1, -30))

    r = ieor(s0, s1)
    state%s(state%ptr) = r

    draw_integer8 = r * spreader
  end function draw_integer8
  
  ! Following are utility functions
  
  subroutine state_init(state, seed)
    implicit none
    type(xorshift1024star_state), intent(out) :: state
    integer, intent(in) :: seed
    integer(8) :: seed64
    integer(8) :: s(16)
    integer :: i

    seed64 = seed
    ! use splimix64 for initialization
    do i = 1, 16
       s(i) = splitmix64(seed64)
    end do

    call state_init_full(state, s)
  end subroutine state_init

  subroutine rand_init(seed)
    implicit none
    integer, intent(in) :: seed
    call state_init(global_state, seed)
  end subroutine rand_init

  subroutine state_init_self(state)
    implicit none
    type(xorshift1024star_state), intent(out) :: state
    integer :: current_time(8)
    !call date_and_time(VALUES=current_time)
    call system_clock(current_time(4))
    if (current_time(4) == -huge(0)) then
       stop "xorshift1024star:state_init_self: too exotic system! (UNIX time not available)"
    endif
    call state_init(state, current_time(4))
  end subroutine state_init_self

  subroutine rand_init_self()
    call state_init_self(global_state)
  end subroutine rand_init_self
  
  ! rmax must be positive
  integer function draw_integer(state, rmax)
    implicit none
    type(xorshift1024star_state), intent(inout) :: state
    integer, intent(in) :: rmax
    integer(8) :: rmax8
    integer(8), parameter :: rmask = 9223372036854775807_8 ! (1 << 63) - 1
    integer(8) :: rnd, qmax, q

    rmax8 = int(rmax, 8)
    ! real maximum is (rmask + 1) / rmax8, but it is impossible to calculate it in 64-bit arithmetic.
    ! Instead we evaluate conservatively. Even if the compiler uses integer(8) as a default integer, 
    ! rmax is maxed at ((1 << 63) - 1), thus qmax >= 1 is ensured.
    qmax = rmask / rmax8
    do
       ! mask to convert to positive integer
       rnd = iand(draw_integer8(state), rmask)
       ! Now both are positive, divide to check whether it is ok to go
       q = rnd / rmax8
       if(q < qmax) then
          draw_integer = int(mod(rnd, rmax8), kind=kind(rmax))
          exit
       endif
       ! otherwise repeat the last step to ensure equidistribution
    end do
  end function draw_integer

  integer(8) function rand_integer(rmax)
    implicit none
    integer, intent(in) :: rmax
    rand_integer = draw_integer(global_state, rmax)
  end function rand_integer

  real(8) function draw_uniform(state)
    implicit none
    type(xorshift1024star_state), intent(inout) :: state
    integer(8) :: rnd 

    ! 1.0 / (1 << 53)
    real(8), parameter :: multiplier = 1.0d0 / 9007199254740992d0

    rnd = draw_integer8(state)
    
    ! 53-bit, divided by 2^53
    draw_uniform = real(ishft(rnd, -11), kind=8) * multiplier
  end function draw_uniform

  real(8) function rand_uniform()
    implicit none
    rand_uniform = draw_uniform(global_state)
  end function rand_uniform

  ! SplitMix64 implementation for random number seed initialization
  integer(8) function splitmix64(seed)
    implicit none
    integer(8), intent(inout) :: seed
    ! As of Fortran95 BOZ literals are available only in DATA stmts.
    ! Also, gfortran does not allow integer(8) greater than 2^63...
    integer(8), parameter :: step = -7046029254386353131_8 ! 0x9E3779B97F4A7C15
    integer(8), parameter :: mix1 = -4658895280553007687_8 ! 0xBF58476D1CE4E5B9
    integer(8), parameter :: mix2 = -7723592293110705685_8 ! 0x94D049BB133111EB

    ! gfortran may issue warning at the following lines.
    ! This is because splitmix64 assumes uint64 wrap-around, which is undefined in F90/95 std.
    ! Even though there are warnings, AFAIK, generated assembler codes are ones as expected.
    seed = seed + step
    splitmix64 = seed
    splitmix64 = ieor(splitmix64, ishft(splitmix64, -30)) * mix1
    splitmix64 = ieor(splitmix64, ishft(splitmix64, -27)) * mix2
    splitmix64 = ieor(splitmix64, ishft(splitmix64, -31))
  end function splitmix64

end module xorshift1024star



