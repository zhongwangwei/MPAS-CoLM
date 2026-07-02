#include <define.h>

MODULE MOD_SPMD_Task

!-----------------------------------------------------------------------------------------
! !DESCRIPTION:
!
!    SPMD refers to "Single PROGRAM/Multiple Data" parallelization.
!
!    MPAS owns MPI initialization and decomposition. Every MPAS rank stays active
!    as a CoLM compute rank, with rank 0 used only for logging/error aggregation.
!
!    CoLM element ownership is supplied by MPAS cell ownership. Patch/PFT state
!    remains internal to CoLM and is mapped back to the owning element/cell by CoLM.
!
!  Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE

   include 'mpif.h'

   integer, parameter :: p_root = 0

   logical :: p_is_root = .false.
   logical :: p_is_active = .false.
   logical :: p_is_compute = .false.
   logical :: p_is_history_task = .false.

   integer :: p_comm_glb_plus = MPI_COMM_NULL
   integer :: p_iam_glb_plus = -1

   ! MPAS-owned communicator aliases used by legacy CoLM helper APIs.
   integer :: p_comm_glb = MPI_COMM_NULL
   integer :: p_iam_glb = -1
   integer :: p_np_glb = 0

   integer :: p_comm_group = MPI_COMM_NULL
   integer :: p_iam_group = -1
   integer :: p_np_group = 0

   integer :: p_my_group = 0
   integer :: p_address_root = p_root

   integer :: p_comm_active = MPI_COMM_NULL
   integer :: p_iam_active = -1
   integer :: p_np_active = 0

   integer, allocatable :: p_itis_active (:)
   integer, allocatable :: p_address_active (:)

   integer :: p_comm_compute = MPI_COMM_NULL
   integer :: p_iam_compute = -1
   integer :: p_np_compute = 0

   integer, allocatable :: p_itis_compute (:)
   integer, allocatable :: p_address_compute (:)

   integer :: p_address_history_task = -1

   integer :: p_stat (MPI_STATUS_SIZE)
   integer :: p_err = 0

   integer, PUBLIC, parameter :: mpi_tag_size = 1
   integer, PUBLIC, parameter :: mpi_tag_mesg = 2
   integer, PUBLIC, parameter :: mpi_tag_data = 3

   integer  :: MPI_INULL_P(1)
   logical  :: MPI_LNULL_P(1)
   real(r8) :: MPI_RNULL_P(1)

   integer, parameter :: MesgMaxSize = 4194304 ! 4MB

   PUBLIC :: spmd_init
   PUBLIC :: spmd_exit
   PUBLIC :: divide_processes_into_groups
   PUBLIC :: spmd_assign_history_task

CONTAINS

   !-----------------------------------------
   SUBROUTINE spmd_init (MyComm_r)

   IMPLICIT NONE
   integer, intent(in), optional :: MyComm_r
   logical mpi_inited
   integer :: iproc

      CALL MPI_INITIALIZED (mpi_inited, p_err)

      IF ( .not. mpi_inited ) THEN
         CALL mpi_init (p_err)
      ENDIF

      IF (present(MyComm_r)) THEN
         CALL MPI_Comm_dup (MyComm_r, p_comm_glb, p_err)
      ELSE
         CALL MPI_Comm_dup (MPI_COMM_WORLD, p_comm_glb, p_err)
      ENDIF

      CALL mpi_comm_rank (p_comm_glb, p_iam_glb, p_err)
      CALL mpi_comm_size (p_comm_glb, p_np_glb,  p_err)

      p_address_root = p_root
      p_is_root = (p_iam_glb == p_address_root)
      p_is_active = .true.
      p_is_compute = .true.
      p_is_history_task = .false.

      p_comm_group = p_comm_glb
      p_iam_group = p_iam_glb
      p_np_group = p_np_glb
      p_my_group = 0

      p_comm_active = p_comm_glb
      p_iam_active = p_iam_glb
      p_np_active = p_np_glb

      p_comm_compute = p_comm_glb
      p_iam_compute = p_iam_glb
      p_np_compute = p_np_glb

      p_comm_glb_plus = MPI_COMM_NULL
      p_iam_glb_plus = -1
      p_address_history_task = -1

      IF (allocated(p_itis_active       )) deallocate (p_itis_active       )
      IF (allocated(p_address_active    )) deallocate (p_address_active    )
      IF (allocated(p_itis_compute   )) deallocate (p_itis_compute   )
      IF (allocated(p_address_compute)) deallocate (p_address_compute)

      allocate (p_itis_active        (0:p_np_glb-1))
      allocate (p_address_active     (0:p_np_glb-1))
      allocate (p_itis_compute    (0:p_np_glb-1))
      allocate (p_address_compute (0:p_np_glb-1))

      DO iproc = 0, p_np_glb-1
         p_itis_active(iproc) = iproc
         p_address_active(iproc) = iproc
         p_itis_compute(iproc) = iproc
         p_address_compute(iproc) = iproc
      ENDDO

   END SUBROUTINE spmd_init

   !-----------------------------------------
   SUBROUTINE divide_processes_into_groups (ngrp)

   IMPLICIT NONE
   integer, intent(in) :: ngrp

      IF (ngrp <= 0) THEN
         CALL mpi_abort (p_comm_glb, 1, p_err)
      ENDIF

      ! MPAS owns the process decomposition. Keep every MPAS rank active.
      p_comm_group = p_comm_glb
      p_iam_group = p_iam_glb
      p_np_group = p_np_glb
      p_my_group = 0

   END SUBROUTINE divide_processes_into_groups

   !-----------------------------------------
   SUBROUTINE spmd_exit

      IF (allocated(p_itis_active       )) deallocate (p_itis_active       )
      IF (allocated(p_address_active    )) deallocate (p_address_active    )
      IF (allocated(p_itis_compute   )) deallocate (p_itis_compute   )
      IF (allocated(p_address_compute)) deallocate (p_address_compute)

      IF (p_comm_glb /= MPI_COMM_NULL) THEN
         CALL mpi_barrier (p_comm_glb, p_err)
         CALL MPI_Comm_free (p_comm_glb, p_err)
      ENDIF

      p_comm_glb = MPI_COMM_NULL
      p_comm_group = MPI_COMM_NULL
      p_comm_active = MPI_COMM_NULL
      p_comm_compute = MPI_COMM_NULL

   END SUBROUTINE spmd_exit

   ! ----- -----
   SUBROUTINE spmd_assign_history_task ()

   IMPLICIT NONE

      CALL CoLM_stop('MPAS embedded CoLM does not support a standalone history MPI task.')

   END SUBROUTINE spmd_assign_history_task

   ! -- STOP all processes --
   SUBROUTINE CoLM_stop (mesg)

   IMPLICIT NONE
   character(len=*), optional :: mesg
   logical :: mpi_inited

      IF (present(mesg)) write(*,*) trim(mesg)

      CALL MPI_INITIALIZED (mpi_inited, p_err)
      IF (mpi_inited .and. p_comm_glb /= MPI_COMM_NULL) THEN
         CALL mpi_abort (p_comm_glb, 1, p_err)
      ENDIF

      STOP

   END SUBROUTINE CoLM_stop

END MODULE MOD_SPMD_Task
