!***************************************************************************************************
!>
!  Module for LSQR kinds and parameters
!
!### History
!  * Jacob Williams : 8 Nov 2019 : created module

    module lsqr_kinds

    use nrtype, only : fp
    !use iso_fortran_env, only: wp => real64  ! double precision kinds

    implicit none

    public

    ! parameters:
    real(fp),parameter :: zero = 0.0_fp
    real(fp),parameter :: one  = 1.0_fp

!***************************************************************************************************
    end module lsqr_kinds
!***************************************************************************************************
