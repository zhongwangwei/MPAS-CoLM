#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeTimeVars
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Time Variables in gridded hydrological processes.
!
! Created by Shupeng Zhang, Oct 2025
!-------------------------------------------------------------------------------------

   USE MOD_Precision
#ifdef GridRiverLakeSediment
   USE MOD_Grid_RiverLakeSediment, only: write_sediment_restart
#endif
   IMPLICIT NONE

   ! -- state variables --
   real(r8), allocatable :: wdsrf_ucat (:) ! river or lake water depth [m]
   real(r8), allocatable :: veloc_riv  (:) ! river velocity            [m/s]
   real(r8), allocatable :: momen_riv  (:) ! unit river momentum       [m^2/s]
   real(r8), allocatable :: volresv    (:) ! reservoir water volume    [m^3]

   ! -- restart file path (saved for deferred sediment restart read) --
   character(len=512) :: gridriver_restart_file = ''

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_GridRiverLakeTimeVars
   PUBLIC :: deallocate_GridRiverLakeTimeVars

   PUBLIC :: read_GridRiverLakeTimeVars
   PUBLIC :: write_GridRiverLakeTimeVars

CONTAINS

   SUBROUTINE allocate_GridRiverLakeTimeVars

   USE MOD_SPMD_Task
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   USE MOD_Grid_Reservoir,        only: numresv
   IMPLICIT NONE

	      IF (p_is_compute) THEN

	         allocate (wdsrf_ucat (numucat))
	         allocate (veloc_riv  (numucat))
	         allocate (momen_riv  (numucat))
	         allocate (volresv    (numresv))

	      ENDIF

   END SUBROUTINE allocate_GridRiverLakeTimeVars


   SUBROUTINE READ_GridRiverLakeTimeVars (file_restart)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Vector_ReadWrite
#ifdef MPAS_EMBEDDED_COLM
   USE MOD_NetCDFSerial, only: ncio_read_indexed_serial
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_ucid
   USE MOD_Grid_Reservoir,        only: numresv, totalnumresv, resv_global_index
#else
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_data_address
   USE MOD_Grid_Reservoir,        only: numresv, resv_data_address, totalnumresv
#endif
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

      gridriver_restart_file = trim(file_restart)

#ifdef MPAS_EMBEDDED_COLM
      IF (p_is_compute .and. numucat > 0) THEN
         CALL ncio_read_indexed_serial (file_restart, 'wdsrf_ucat', ucat_ucid, wdsrf_ucat)
         CALL ncio_read_indexed_serial (file_restart, 'veloc_riv',  ucat_ucid, veloc_riv )
      ENDIF
#else
      CALL vector_read_and_scatter (file_restart, wdsrf_ucat, numucat, 'wdsrf_ucat', ucat_data_address)
      CALL vector_read_and_scatter (file_restart, veloc_riv,  numucat, 'veloc_riv',  ucat_data_address)
#endif

      IF (DEF_Reservoir_Method > 0) THEN
         IF (totalnumresv > 0) THEN
#ifdef MPAS_EMBEDDED_COLM
            IF (p_is_compute .and. numresv > 0) THEN
               CALL ncio_read_indexed_serial (file_restart, 'volresv', resv_global_index, volresv)
            ENDIF
#else
            CALL vector_read_and_scatter (file_restart, volresv, numresv, 'volresv', resv_data_address)
#endif
         ENDIF
      ENDIF

      ! Note: sediment restart is read separately in grid_sediment_read_restart,
      ! called from grid_riverlake_flow_init after sediment module is initialized.

   END SUBROUTINE READ_GridRiverLakeTimeVars


   SUBROUTINE WRITE_GridRiverLakeTimeVars (file_restart)

   USE MOD_SPMD_Task
   USE MOD_Namelist
#ifndef MPAS_EMBEDDED_COLM
   USE MOD_NetCDFSerial
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address
   USE MOD_Grid_Reservoir,        only: numresv, totalnumresv, resv_data_address
#endif
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

#ifdef MPAS_EMBEDDED_COLM
      CALL write_gridriver_restart_mpas_embedded (file_restart)
#else
      IF (p_is_root) THEN
         CALL ncio_create_file (trim(file_restart))
         CALL ncio_define_dimension(file_restart, 'ucatch', totalnumucat)
      ENDIF

      CALL vector_gather_and_write (&
         wdsrf_ucat, numucat, totalnumucat, ucat_data_address, file_restart, 'wdsrf_ucat', 'ucatch')

      CALL vector_gather_and_write (&
         veloc_riv, numucat, totalnumucat, ucat_data_address, file_restart, 'veloc_riv', 'ucatch')

      IF (DEF_Reservoir_Method > 0) THEN
         IF (totalnumresv > 0) THEN

            IF (p_is_root) CALL ncio_define_dimension(file_restart, 'reservoir', totalnumresv)

            CALL vector_gather_and_write (&
               volresv, numresv, totalnumresv, resv_data_address, file_restart, 'volresv', 'reservoir')
         ENDIF
      ENDIF
#endif

#ifdef GridRiverLakeSediment
      IF (DEF_USE_SEDIMENT) THEN
         CALL write_sediment_restart(file_restart)
      ENDIF
#endif

   END SUBROUTINE WRITE_GridRiverLakeTimeVars

#ifdef MPAS_EMBEDDED_COLM
   SUBROUTINE write_gridriver_restart_mpas_embedded (file_restart)

   USE mpi, only: MPI_INFO_NULL, MPI_OFFSET_KIND
   USE pnetcdf
   USE MOD_SPMD_Task, only: p_comm_compute
   USE MOD_Namelist, only: DEF_Reservoir_Method
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_ucid
   USE MOD_Grid_Reservoir,        only: numresv, totalnumresv, resv_global_index
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

   integer :: ierr
   integer :: ncid
   integer :: dim_ucatch
   integer :: dim_reservoir
   integer :: var_wdsrf
   integer :: var_veloc
   integer :: var_volresv
   logical :: write_reservoir

      write_reservoir = DEF_Reservoir_Method > 0 .and. totalnumresv > 0

      ierr = nf90mpi_create(p_comm_compute, trim(file_restart), &
         IOR(NF90_CLOBBER, NF90_64BIT_OFFSET), MPI_INFO_NULL, ncid)
      CALL pnetcdf_check(ierr, 'create', file_restart)

      ierr = nf90mpi_def_dim(ncid, 'ucatch', int(totalnumucat, MPI_OFFSET_KIND), dim_ucatch)
      CALL pnetcdf_check(ierr, 'define ucatch dimension', file_restart)

      ierr = nf90mpi_def_var(ncid, 'wdsrf_ucat', NF90_DOUBLE, (/dim_ucatch/), var_wdsrf)
      CALL pnetcdf_check(ierr, 'define wdsrf_ucat', file_restart)

      ierr = nf90mpi_def_var(ncid, 'veloc_riv', NF90_DOUBLE, (/dim_ucatch/), var_veloc)
      CALL pnetcdf_check(ierr, 'define veloc_riv', file_restart)

      IF (write_reservoir) THEN
         ierr = nf90mpi_def_dim(ncid, 'reservoir', int(totalnumresv, MPI_OFFSET_KIND), dim_reservoir)
         CALL pnetcdf_check(ierr, 'define reservoir dimension', file_restart)

         ierr = nf90mpi_def_var(ncid, 'volresv', NF90_DOUBLE, (/dim_reservoir/), var_volresv)
         CALL pnetcdf_check(ierr, 'define volresv', file_restart)
      ENDIF

      ierr = nf90mpi_enddef(ncid)
      CALL pnetcdf_check(ierr, 'end define mode', file_restart)

      ierr = nf90mpi_begin_indep_data(ncid)
      CALL pnetcdf_check(ierr, 'begin independent data mode', file_restart)

      CALL pnetcdf_write_real8_points(ncid, var_wdsrf, ucat_ucid, wdsrf_ucat, numucat, &
         'wdsrf_ucat', file_restart)
      CALL pnetcdf_write_real8_points(ncid, var_veloc, ucat_ucid, veloc_riv, numucat, &
         'veloc_riv', file_restart)

      IF (write_reservoir) THEN
         CALL pnetcdf_write_real8_points(ncid, var_volresv, resv_global_index, volresv, numresv, &
            'volresv', file_restart)
      ENDIF

      ierr = nf90mpi_end_indep_data(ncid)
      CALL pnetcdf_check(ierr, 'end independent data mode', file_restart)

      ierr = nf90mpi_close(ncid)
      CALL pnetcdf_check(ierr, 'close', file_restart)

   END SUBROUTINE write_gridriver_restart_mpas_embedded

   SUBROUTINE pnetcdf_write_real8_points(ncid, varid, index, data, ndata, varname, filename)

   USE mpi, only: MPI_OFFSET_KIND
   USE pnetcdf, only: nf90mpi_put_var
   USE MOD_SPMD_Task, only: CoLM_stop
   IMPLICIT NONE

   integer, intent(in) :: ncid
   integer, intent(in) :: varid
   integer, intent(in) :: index(:)
   real(r8), intent(in) :: data(:)
   integer, intent(in) :: ndata
   character(len=*), intent(in) :: varname
   character(len=*), intent(in) :: filename

   integer :: ierr
   integer :: i
   integer(kind=MPI_OFFSET_KIND) :: start(1)
   real(r8) :: value

      IF (ndata > size(index) .or. ndata > size(data)) THEN
         CALL CoLM_stop('PnetCDF indexed write size mismatch for '//trim(varname))
      ENDIF

      DO i = 1, ndata
         IF (index(i) < 1) THEN
            CALL CoLM_stop('PnetCDF indexed write invalid index for '//trim(varname))
         ENDIF

         start(1) = int(index(i), MPI_OFFSET_KIND)
         value = data(i)
         ierr = nf90mpi_put_var(ncid, varid, value, start=start)
         CALL pnetcdf_check(ierr, 'write '//trim(varname), filename)
      ENDDO

   END SUBROUTINE pnetcdf_write_real8_points

   SUBROUTINE pnetcdf_check(status, action, filename)

   USE pnetcdf, only: NF90_NOERR, nf90mpi_strerror
   USE MOD_SPMD_Task, only: CoLM_stop
   IMPLICIT NONE

   integer, intent(in) :: status
   character(len=*), intent(in) :: action
   character(len=*), intent(in) :: filename

      IF (status /= NF90_NOERR) THEN
         write(*,'(A)') 'PnetCDF error during '//trim(action)//' for '//trim(filename)//': ' &
            //trim(nf90mpi_strerror(status))
         CALL CoLM_stop()
      ENDIF

   END SUBROUTINE pnetcdf_check
#endif

   SUBROUTINE deallocate_GridRiverLakeTimeVars

   IMPLICIT NONE

      IF (allocated (wdsrf_ucat)) deallocate (wdsrf_ucat)
      IF (allocated (veloc_riv )) deallocate (veloc_riv )
      IF (allocated (momen_riv )) deallocate (momen_riv )
      IF (allocated (volresv   )) deallocate (volresv   )

   END SUBROUTINE deallocate_GridRiverLakeTimeVars

END MODULE MOD_Grid_RiverLakeTimeVars
#endif
