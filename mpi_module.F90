module mpi_module
  !! module contains common variables and subroutines related to openmpi
  !!Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
  !!Copyright: (c) Masashi Ogiso 2022
  !!License  : MIT License (https://opensource.org/licenses/MIT)

  use nrtype, only : fp, sp, dp
  use mpi

  implicit none
  private
  integer :: mpi_ierr, mpi_errcode, dum_array_size, mpi_ilen

  integer, public :: mpi_rank, mpi_size
  character(len = 16), public :: mpi_hostname

  public :: mpi_ini_rank_size, mpifinalize, mpiabort, &
  &         mpi_bcast_int, mpi_bcast_fp, mpi_bcast_dp, &
  &         mpi_bcast_int_1d, mpi_bcast_int_2d, &
  &         mpi_bcast_fp_1d, mpi_bcast_fp_2d, &
  &         mpi_bcast_fp_3d, mpi_bcast_fp_4d, &
  &         mpi_bcast_dp_1d, mpi_bcast_dp_2d, &
  &         mpi_bcast_dp_3d, mpi_bcast_logical_1d, &
  &         mpi_bcast_logical_2d, &
  &         mpi_get_hostname, mpi_barrier_sub, &
  &         mpi_allreduce_fp_1d, mpi_allreduce_dp_1d, &
  &         mpi_allreduce_fp_2d, mpi_allreduce_dp_2d, &
  &         mpi_allreduce_fp_3d, mpi_allreduce_dp_3d, &
  &         mpi_allreduce_sp_3d, &
  &         mpi_allreduce_fp_4d, mpi_allreduce_dp_4d, &
  &         mpi_allreduce_dp_5d, mpi_allreduce_dp_6d, &
  &         mpi_allreduce_int_1d, mpi_allreduce_int_3d, &
  &         mpi_allreduce_int_4d, mpi_allreduce_fp_min, &
  &         mpi_allreduce_lor_1d

contains

#ifdef DOUBLE
#define MPI_REAL_KIND MPI_REAL8
#else
#define MPI_REAL_KIND MPI_REAL
#endif

!!MPI_INIT, MPI_COMM_SIZE, MPI_COMM_RANK
  subroutine mpi_ini_rank_size

    call mpi_init(mpi_ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_ierr)
    call mpi_comm_size(MPI_COMM_WORLD, mpi_size, mpi_ierr)

    return
  end subroutine mpi_ini_rank_size

!!MPI_FINALIZE
  subroutine mpifinalize
    call mpi_finalize(mpi_ierr)
    return
  end subroutine mpifinalize

!!MPI_Abort
  subroutine mpiabort
    call mpi_abort(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
    return
  end subroutine mpiabort

!!MPI_BCAST for integer
!!single variable
  subroutine mpi_bcast_int(dum_variable, sendrank)
    integer, intent(inout) :: dum_variable
    integer, intent(in) :: sendrank
    
    call mpi_bcast(dum_variable, 1, MPI_INTEGER, sendrank, MPI_COMM_WORLD, mpi_ierr)
  
    return
  end subroutine mpi_bcast_int

  subroutine mpi_bcast_fp(dum_variable, sendrank)
    real(kind = fp), intent(inout) :: dum_variable
    integer, intent(in) :: sendrank
    
    call mpi_bcast(dum_variable, 1, MPI_REAL_KIND, sendrank, MPI_COMM_WORLD, mpi_ierr)
  
    return
  end subroutine mpi_bcast_fp

  subroutine mpi_bcast_dp(dum_variable, sendrank)
    real(kind = dp), intent(inout) :: dum_variable
    integer, intent(in) :: sendrank
    
    call mpi_bcast(dum_variable, 1, MPI_REAL8, sendrank, MPI_COMM_WORLD, mpi_ierr)
  
    return
  end subroutine mpi_bcast_dp

!!input should be array 1d, 2d
  subroutine mpi_bcast_int_1d(dum_array_int, sendrank)
    integer, intent(inout) :: dum_array_int(:)
    integer, intent(in) :: sendrank
    
    dum_array_size = size(dum_array_int)
    call mpi_bcast(dum_array_int, dum_array_size, MPI_INTEGER, sendrank, MPI_COMM_WORLD, mpi_ierr)
  
    return
  end subroutine mpi_bcast_int_1d

  subroutine mpi_bcast_int_2d(dum_array_int, sendrank)
    integer, intent(inout) :: dum_array_int(:, :)
    integer, intent(in) :: sendrank
    
    dum_array_size = size(dum_array_int)
    call mpi_bcast(dum_array_int, dum_array_size, MPI_INTEGER, sendrank, MPI_COMM_WORLD, mpi_ierr)
  
    return
  end subroutine mpi_bcast_int_2d

!!MPI_BCAST for real
!!input should be array 1d, 2d
  subroutine mpi_bcast_fp_1d(dum_array_real, sendrank)
    real(fp), intent(inout) :: dum_array_real(:)
    integer, intent(in) :: sendrank

    dum_array_size = size(dum_array_real)
    call mpi_bcast(dum_array_real, dum_array_size, MPI_REAL_KIND, sendrank, MPI_COMM_WORLD, mpi_ierr)
  
    return
  end subroutine mpi_bcast_fp_1d

  subroutine mpi_bcast_fp_2d(dum_array_real, sendrank)
    real(fp), intent(inout) :: dum_array_real(:, :)
    integer, intent(in) :: sendrank

    dum_array_size = size(dum_array_real)
    call mpi_bcast(dum_array_real, dum_array_size, MPI_REAL_KIND, sendrank, MPI_COMM_WORLD, mpi_ierr)
  
    return
  end subroutine mpi_bcast_fp_2d

  subroutine mpi_bcast_fp_3d(dum_array_real, sendrank)
    real(fp), intent(inout) :: dum_array_real(:, :, :)
    integer, intent(in) :: sendrank

    dum_array_size = size(dum_array_real)
    call mpi_bcast(dum_array_real, dum_array_size, MPI_REAL_KIND, sendrank, MPI_COMM_WORLD, mpi_ierr)
  
    return
  end subroutine mpi_bcast_fp_3d

  subroutine mpi_bcast_fp_4d(dum_array_real, sendrank)
    real(fp), intent(inout) :: dum_array_real(:, :, :, :)
    integer, intent(in) :: sendrank

    dum_array_size = size(dum_array_real)
    call mpi_bcast(dum_array_real, dum_array_size, MPI_REAL_KIND, sendrank, MPI_COMM_WORLD, mpi_ierr)
  
    return
  end subroutine mpi_bcast_fp_4d

!!MPI_BCAST for double
!!input should be array 1d, 2d
  subroutine mpi_bcast_dp_1d(dum_array_real, sendrank)
    real(dp), intent(inout) :: dum_array_real(:)
    integer, intent(in) :: sendrank

    dum_array_size = size(dum_array_real)
    call mpi_bcast(dum_array_real, dum_array_size, MPI_REAL8, sendrank, MPI_COMM_WORLD, mpi_ierr)
  
    return
  end subroutine mpi_bcast_dp_1d

  subroutine mpi_bcast_dp_2d(dum_array_real, sendrank)
    real(dp), intent(inout) :: dum_array_real(:, :)
    integer, intent(in) :: sendrank

    dum_array_size = size(dum_array_real)
    call mpi_bcast(dum_array_real, dum_array_size, MPI_REAL8, sendrank, MPI_COMM_WORLD, mpi_ierr)
  
    return
  end subroutine mpi_bcast_dp_2d

  subroutine mpi_bcast_dp_3d(dum_array_real, sendrank)
    real(dp), intent(inout) :: dum_array_real(:, :, :)
    integer, intent(in) :: sendrank

    dum_array_size = size(dum_array_real)
    call mpi_bcast(dum_array_real, dum_array_size, MPI_REAL8, sendrank, MPI_COMM_WORLD, mpi_ierr)
  
    return
  end subroutine mpi_bcast_dp_3d

  subroutine mpi_bcast_logical_1d(dum_array_real, sendrank)
    logical, intent(inout) :: dum_array_real(:)
    integer, intent(in) :: sendrank

    dum_array_size = size(dum_array_real)
    call mpi_bcast(dum_array_real, dum_array_size, MPI_LOGICAL, sendrank, MPI_COMM_WORLD, mpi_ierr)
  
    return
  end subroutine mpi_bcast_logical_1d

  subroutine mpi_bcast_logical_2d(dum_array_real, sendrank)
    logical, intent(inout) :: dum_array_real(:, :)
    integer, intent(in) :: sendrank

    dum_array_size = size(dum_array_real)
    call mpi_bcast(dum_array_real, dum_array_size, MPI_LOGICAL, sendrank, MPI_COMM_WORLD, mpi_ierr)
  
    return
  end subroutine mpi_bcast_logical_2d

!!MPI_GET_PROCESSOR_NAME
  subroutine mpi_get_hostname
    call mpi_get_processor_name(mpi_hostname, mpi_ilen, mpi_ierr)
    return
  end subroutine mpi_get_hostname

!!MPI_BARRIER
  subroutine mpi_barrier_sub
    call mpi_barrier(MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_barrier_sub

!!MPI_ALLREDUCE for same send/recv array with MPI_SUM
  subroutine mpi_allreduce_fp_1d(dum_array_real)
    real(fp), intent(inout) :: dum_array_real(:)
    
    dum_array_size = size(dum_array_real)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_real, dum_array_size, MPI_REAL_KIND, MPI_SUM, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_fp_1d

  subroutine mpi_allreduce_dp_1d(dum_array_real)
    real(dp), intent(inout) :: dum_array_real(:)
    
    dum_array_size = size(dum_array_real)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_real, dum_array_size, MPI_REAL8, MPI_SUM, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_dp_1d

  subroutine mpi_allreduce_fp_2d(dum_array_real)
    real(fp), intent(inout) :: dum_array_real(:, :)
    
    dum_array_size = size(dum_array_real)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_real, dum_array_size, MPI_REAL_KIND, MPI_SUM, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_fp_2d

  subroutine mpi_allreduce_dp_2d(dum_array_real)
    real(dp), intent(inout) :: dum_array_real(:, :)
    
    dum_array_size = size(dum_array_real)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_real, dum_array_size, MPI_REAL8, MPI_SUM, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_dp_2d

  subroutine mpi_allreduce_fp_3d(dum_array_real)
    real(fp), intent(inout) :: dum_array_real(:, :, :)
    
    dum_array_size = size(dum_array_real)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_real, dum_array_size, MPI_REAL_KIND, MPI_SUM, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_fp_3d

  subroutine mpi_allreduce_dp_3d(dum_array_real)
    real(dp), intent(inout) :: dum_array_real(:, :, :)
    
    dum_array_size = size(dum_array_real)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_real, dum_array_size, MPI_REAL8, MPI_SUM, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_dp_3d

  subroutine mpi_allreduce_sp_3d(dum_array_real)
    real(sp), intent(inout) :: dum_array_real(:, :, :)
    
    dum_array_size = size(dum_array_real)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_real, dum_array_size, MPI_REAL, MPI_SUM, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_sp_3d

  subroutine mpi_allreduce_fp_4d(dum_array_real)
    real(fp), intent(inout) :: dum_array_real(:, :, :, :)
    
    dum_array_size = size(dum_array_real)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_real, dum_array_size, MPI_REAL_KIND, MPI_SUM, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_fp_4d

  subroutine mpi_allreduce_dp_4d(dum_array_real)
    real(dp), intent(inout) :: dum_array_real(:, :, :, :)
    
    dum_array_size = size(dum_array_real)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_real, dum_array_size, MPI_REAL8, MPI_SUM, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_dp_4d

  subroutine mpi_allreduce_dp_5d(dum_array_real)
    real(dp), intent(inout) :: dum_array_real(:, :, :, :, :)
    
    dum_array_size = size(dum_array_real)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_real, dum_array_size, MPI_REAL8, MPI_SUM, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_dp_5d

  subroutine mpi_allreduce_dp_6d(dum_array_real)
    real(dp), intent(inout) :: dum_array_real(:, :, :, :, :, :)
    
    dum_array_size = size(dum_array_real)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_real, dum_array_size, MPI_REAL8, MPI_SUM, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_dp_6d

  subroutine mpi_allreduce_int_1d(dum_array_int)
    integer, intent(inout) :: dum_array_int(:)
    
    dum_array_size = size(dum_array_int)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_int, dum_array_size, MPI_INTEGER, MPI_SUM, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_int_1d

  subroutine mpi_allreduce_int_3d(dum_array_int)
    integer, intent(inout) :: dum_array_int(:, :, :)
    
    dum_array_size = size(dum_array_int)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_int, dum_array_size, MPI_INTEGER, MPI_SUM, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_int_3d

  subroutine mpi_allreduce_int_4d(dum_array_int)
    integer, intent(inout) :: dum_array_int(:, :, :, :)
    
    dum_array_size = size(dum_array_int)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_int, dum_array_size, MPI_INTEGER, MPI_SUM, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_int_4d

  subroutine mpi_allreduce_lor_1d(dum_array_logical)
    logical, intent(inout) :: dum_array_logical(:)

    dum_array_size = size(dum_array_logical)
    call mpi_allreduce(MPI_IN_PLACE, dum_array_logical, dum_array_size, MPI_LOGICAL, MPI_LOR, &
    &                  MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_lor_1d

  subroutine mpi_allreduce_fp_min(dum_val)
    real(fp), intent(inout) :: dum_val

    call mpi_allreduce(MPI_IN_PLACE, dum_val, 1, MPI_REAL_KIND, MPI_MIN, MPI_COMM_WORLD, mpi_ierr)
    return
  end subroutine mpi_allreduce_fp_min



end module mpi_module

