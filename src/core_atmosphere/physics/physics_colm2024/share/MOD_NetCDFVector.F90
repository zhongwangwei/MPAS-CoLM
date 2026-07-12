#include <define.h>

!----------------------------------------------------------------------------------
! !DESCRIPTION:
!
!    High-level Subroutines to read and write variables in files with netCDF format.
!
!    CoLM read and write netCDF files mainly in three ways:
!    1. Serial: read and write data by a single process;
!    2. Vector: read/write data associated with CoLM pixelsets
!    3. Block : read blocked data by IO
!               Notice: input file is a single file.
!
!    This MODULE CONTAINS subroutines of "2. Vector".
!
!    Two implementations can be used,
!    1) "MOD_NetCDFVectorBlk.F90":
!       A vector is saved in separated files, each associated with a block.
!       READ/WRITE are fast in this way and compression can be used.
!       However, there may be too many files, especially when blocks are small.
!       CHOOSE this implementation by "#undef VectorInOneFile" in include/define.h
!    2) "MOD_NetCDFVectorOne.F90":
!       A vector is saved in one file.
!       READ/WRITE may be slow in this way.
!       CHOOSE this implementation by "#define VectorInOneFile" in include/define.h
!
!  Created by Shupeng Zhang, May 2023
!----------------------------------------------------------------------------------

MODULE MOD_NetCDFVector


   USE MOD_DataType
#ifdef MPAS_EMBEDDED_COLM
   USE mpi, only: MPI_Barrier
   USE MOD_MPAS_MPI, only: CoLM_stop
#endif
   IMPLICIT NONE

#ifdef MPAS_EMBEDDED_COLM
   type :: distributed_vector_map_type
      character(len=512) :: filename = ''
      integer :: iblk = 0
      integer :: jblk = 0
      integer :: istt = 0
      integer :: iend = -1
      integer*8, allocatable :: eindex(:)
      integer, allocatable :: ipxstt(:)
      integer, allocatable :: ipxend(:)
      integer, allocatable :: settyp(:)
      integer, allocatable :: source_rank(:)
      integer, allocatable :: source_pos(:)
      integer, allocatable :: source_ranks(:)
      integer, allocatable :: source_lengths(:)
   END type distributed_vector_map_type

   type(distributed_vector_map_type), allocatable, save :: distributed_vector_maps(:)
   character(len=512), allocatable, save :: validated_distributed_markers(:)
   integer, allocatable, save :: distributed_marker_ranks(:)
#endif

   ! PUBLIC subroutines

   INTERFACE ncio_read_vector
      MODULE procedure ncio_read_vector_logical_1d
      MODULE procedure ncio_read_vector_int32_1d
      MODULE procedure ncio_read_vector_int64_1d
      MODULE procedure ncio_read_vector_real8_1d
      MODULE procedure ncio_read_vector_real8_2d
      MODULE procedure ncio_read_vector_real8_3d
      MODULE procedure ncio_read_vector_real8_4d
      MODULE procedure ncio_read_vector_real8_5d
   END INTERFACE ncio_read_vector

   PUBLIC :: ncio_create_file_vector
   PUBLIC :: ncio_define_dimension_vector
   PUBLIC :: ncio_begin_distributed_write
   PUBLIC :: ncio_complete_distributed_write
   PUBLIC :: ncio_reset_distributed_validation

   INTERFACE ncio_write_vector
      MODULE procedure ncio_write_vector_logical_1d
      MODULE procedure ncio_write_vector_int32_1d
      MODULE procedure ncio_write_vector_int32_3d
      MODULE procedure ncio_write_vector_int64_1d
      MODULE procedure ncio_write_vector_real8_1d
      MODULE procedure ncio_write_vector_real8_2d
      MODULE procedure ncio_write_vector_real8_3d
      MODULE procedure ncio_write_vector_real8_4d
      MODULE procedure ncio_write_vector_real8_5d
   END INTERFACE ncio_write_vector

#ifdef MPAS_EMBEDDED_COLM
   INTERFACE read_distributed_vector_block
      MODULE procedure read_distributed_int32_1d
      MODULE procedure read_distributed_int64_1d
      MODULE procedure read_distributed_int8_1d
      MODULE procedure read_distributed_real8_1d
      MODULE procedure read_distributed_real8_2d
      MODULE procedure read_distributed_real8_3d
      MODULE procedure read_distributed_real8_4d
      MODULE procedure read_distributed_real8_5d
   END INTERFACE read_distributed_vector_block
#endif

CONTAINS

   !---------------------------------------------------------
   SUBROUTINE get_filename_vector_block (filename, iblk, jblk, fileblock, for_write, use_srcpos, distributed_restart)

   USE MOD_Block, only: get_filename_block
#ifdef MPAS_EMBEDDED_COLM
   USE MOD_MPAS_MPI, only: mpas_rank, CoLM_stop
#endif
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   integer, intent(in) :: iblk, jblk
   character(len=*), intent(out) :: fileblock
   logical, intent(in), optional :: for_write
   logical, intent(inout), optional :: use_srcpos
   logical, intent(out), optional :: distributed_restart

#ifdef MPAS_EMBEDDED_COLM
   character(len=512) :: complete_marker
   character(len=256) :: fileblock_rank
   logical :: marker_exists, rank_file_exists, distributed_file
#endif

      CALL get_filename_block (filename, iblk, jblk, fileblock)
      IF (present(distributed_restart)) distributed_restart = .false.

#ifdef MPAS_EMBEDDED_COLM
      CALL add_rank_suffix(fileblock, mpas_rank, fileblock_rank)

      distributed_file = .false.
      IF (present(for_write)) distributed_file = for_write
      IF ((.not. distributed_file) .and. present(use_srcpos)) THEN
         complete_marker = trim(filename)//'.mpas_complete'
         inquire(file=trim(complete_marker), exist=marker_exists)
         IF (marker_exists) THEN
            distributed_file = .true.
         ELSE
            inquire(file=trim(fileblock_rank), exist=rank_file_exists)
            IF (rank_file_exists) THEN
               CALL CoLM_stop('MPAS embedded CoLM found an incomplete distributed restart without its completion marker: '// &
                              trim(fileblock_rank))
            ENDIF
         ENDIF
      ENDIF
      IF (present(for_write)) THEN
         IF (distributed_file) fileblock = fileblock_rank
      ELSEIF (distributed_file) THEN
         ! A distributed restart can be repartitioned.  The caller selects the
         ! old rank shards through the stable pixelset identity map below.
         fileblock = fileblock_rank
         IF (present(use_srcpos)) use_srcpos = .false.
      ENDIF
      IF (present(distributed_restart)) distributed_restart = distributed_file
#endif

	   END SUBROUTINE get_filename_vector_block

	   !---------------------------------------------------------
	   SUBROUTINE add_rank_suffix(fileblock, rank, fileblock_rank)

	   IMPLICIT NONE

	   character(len=*), intent(in) :: fileblock
	   integer, intent(in) :: rank
	   character(len=*), intent(out) :: fileblock_rank

	   character(len=32) :: rank_suffix
	   integer :: idot

	      write(rank_suffix,'("_mpasr",I8.8)') rank
	      idot = len_trim(fileblock)
	      DO WHILE (idot > 0)
	         IF (fileblock(idot:idot) == '.') EXIT
	         idot = idot - 1
	      ENDDO

	      IF (idot > 0) THEN
	         fileblock_rank = fileblock(1:idot-1) // trim(rank_suffix) // fileblock(idot:len_trim(fileblock))
	      ELSE
	         fileblock_rank = trim(fileblock) // trim(rank_suffix)
	      ENDIF

	   END SUBROUTINE add_rank_suffix

	   !---------------------------------------------------------
	   SUBROUTINE read_distributed_marker(complete_marker, marker_ranks)

#ifdef MPAS_EMBEDDED_COLM
	   USE MOD_MPAS_MPI, only: CoLM_stop
#endif
	   IMPLICIT NONE

	   character(len=*), intent(in) :: complete_marker
	   integer, intent(out) :: marker_ranks

#ifdef MPAS_EMBEDDED_COLM
	   character(len=512) :: line
	   character(len=512), allocatable :: markers_new(:)
	   integer, allocatable :: marker_ranks_new(:)
	   integer :: i, iunit, ios, nvalidated, separator

	      IF (allocated(validated_distributed_markers)) THEN
	         DO i = 1, size(validated_distributed_markers)
	            IF (trim(validated_distributed_markers(i)) == trim(complete_marker)) THEN
	               marker_ranks = distributed_marker_ranks(i)
	               RETURN
	            ENDIF
	         ENDDO
	      ENDIF

	      open(newunit=iunit, file=trim(complete_marker), status='old', action='read', iostat=ios)
	      IF (ios /= 0) CALL CoLM_stop('Cannot read MPAS embedded CoLM restart marker: '//trim(complete_marker))
	      read(iunit,'(A)',iostat=ios) line
	      close(iunit)
	      separator = index(line, ':')
	      IF (ios /= 0 .or. separator < 1 .or. separator >= len_trim(line)) THEN
	         CALL CoLM_stop('Invalid MPAS embedded CoLM restart marker: '//trim(complete_marker))
	      ENDIF
	      read(line(separator+1:),*,iostat=ios) marker_ranks
	      IF (ios /= 0 .or. marker_ranks < 1) THEN
	         CALL CoLM_stop('Invalid MPI rank count in MPAS embedded CoLM restart marker: '//trim(complete_marker))
	      ENDIF
	      IF (allocated(validated_distributed_markers)) THEN
	         nvalidated = size(validated_distributed_markers)
	         allocate(markers_new(nvalidated + 1))
	         allocate(marker_ranks_new(nvalidated + 1))
	         markers_new(1:nvalidated) = validated_distributed_markers
	         marker_ranks_new(1:nvalidated) = distributed_marker_ranks
	         markers_new(nvalidated + 1) = complete_marker
	         marker_ranks_new(nvalidated + 1) = marker_ranks
	         CALL move_alloc(markers_new, validated_distributed_markers)
	         CALL move_alloc(marker_ranks_new, distributed_marker_ranks)
	      ELSE
	         allocate(validated_distributed_markers(1))
	         allocate(distributed_marker_ranks(1))
	         validated_distributed_markers(1) = complete_marker
	         distributed_marker_ranks(1) = marker_ranks
	      ENDIF
#endif

	   END SUBROUTINE read_distributed_marker

	   !---------------------------------------------------------
	   SUBROUTINE remove_rank_distributed_shards(filename)

#ifdef MPAS_EMBEDDED_COLM
	   USE MOD_Block, only: gblock, get_filename_block
	   USE MOD_MPAS_MPI, only: mpas_rank, CoLM_stop
#endif
	   IMPLICIT NONE

	   character(len=*), intent(in) :: filename

#ifdef MPAS_EMBEDDED_COLM
	   character(len=512) :: fileblock, fileblock_rank
	   integer :: iblk, jblk, iunit, iostat_open, iostat_close
	   logical :: shard_exists

	      IF (gblock%nxblk < 1 .or. gblock%nyblk < 1) THEN
	         CALL CoLM_stop('Cannot clean distributed CoLM restart shards before block metadata is initialized.')
	      ENDIF

	      DO jblk = 1, gblock%nyblk
	         DO iblk = 1, gblock%nxblk
	            CALL get_filename_block(filename, iblk, jblk, fileblock)
	            CALL add_rank_suffix(fileblock, mpas_rank, fileblock_rank)
	            inquire(file=trim(fileblock_rank), exist=shard_exists)
	            IF (.not. shard_exists) CYCLE

	            open(newunit=iunit, file=trim(fileblock_rank), status='old', action='readwrite', iostat=iostat_open)
	            IF (iostat_open /= 0) THEN
	               CALL CoLM_stop('Cannot open stale MPAS embedded CoLM restart shard for removal: '// &
	                              trim(fileblock_rank))
	            ENDIF
	            close(iunit, status='delete', iostat=iostat_close)
	            IF (iostat_close /= 0) THEN
	               CALL CoLM_stop('Cannot remove stale MPAS embedded CoLM restart shard: '//trim(fileblock_rank))
	            ENDIF
	         ENDDO
	      ENDDO
#endif

	   END SUBROUTINE remove_rank_distributed_shards

	   !---------------------------------------------------------
	   SUBROUTINE ncio_begin_distributed_write(filename)

#ifdef MPAS_EMBEDDED_COLM
	   USE MOD_MPAS_MPI, only: mpas_is_root, mpas_comm, mpas_mpi_ierr, CoLM_stop, mpas_mpi_check
#endif
	   IMPLICIT NONE

	   character(len=*), intent(in) :: filename

#ifdef MPAS_EMBEDDED_COLM
	   character(len=512) :: complete_marker
	   integer :: iunit, iostat_open, iostat_close
	   logical :: marker_exists

	      complete_marker = trim(filename)//'.mpas_complete'
	      IF (mpas_is_root) THEN
	         inquire(file=trim(complete_marker), exist=marker_exists)
	         IF (marker_exists) THEN
	            open(newunit=iunit, file=trim(complete_marker), status='old', action='readwrite', iostat=iostat_open)
	            IF (iostat_open /= 0) CALL CoLM_stop('Cannot replace MPAS embedded CoLM restart marker: '// &
	                                                trim(complete_marker))
	            close(iunit, status='delete', iostat=iostat_close)
	            IF (iostat_close /= 0) CALL CoLM_stop('Cannot remove stale MPAS embedded CoLM restart marker: '// &
	                                                  trim(complete_marker))
	         ENDIF
	      ENDIF
	      CALL mpi_barrier(mpas_comm, mpas_mpi_ierr)
	      CALL mpas_mpi_check('distributed restart marker removal')
	      CALL remove_rank_distributed_shards(filename)
	      CALL mpi_barrier(mpas_comm, mpas_mpi_ierr)
	      CALL mpas_mpi_check('distributed restart shard cleanup')
#endif

	   END SUBROUTINE ncio_begin_distributed_write

	   !---------------------------------------------------------
	   SUBROUTINE ncio_complete_distributed_write(filename)

#ifdef MPAS_EMBEDDED_COLM
	   USE MOD_MPAS_MPI, only: mpas_is_root, mpas_comm, mpas_mpi_ierr, mpas_size, CoLM_stop, mpas_mpi_check
#endif
	   IMPLICIT NONE

	   character(len=*), intent(in) :: filename

#ifdef MPAS_EMBEDDED_COLM
	   character(len=512) :: complete_marker
	   integer :: iunit, iostat_open, iostat_close

	      CALL mpi_barrier(mpas_comm, mpas_mpi_ierr)
	      CALL mpas_mpi_check('distributed restart shard completion')
	      complete_marker = trim(filename)//'.mpas_complete'
	      IF (mpas_is_root) THEN
	         open(newunit=iunit, file=trim(complete_marker), status='replace', action='write', iostat=iostat_open)
	         IF (iostat_open /= 0) CALL CoLM_stop('Cannot create MPAS embedded CoLM restart marker: '// &
	                                             trim(complete_marker))
	         write(iunit,'(A,I0)') 'MPAS ranks: ', mpas_size
	         write(iunit,'(A)') 'Format: distributed-pixelset-v1'
	         close(iunit, iostat=iostat_close)
	         IF (iostat_close /= 0) CALL CoLM_stop('Cannot close MPAS embedded CoLM restart marker: '// &
	                                               trim(complete_marker))
	      ENDIF
	      CALL mpi_barrier(mpas_comm, mpas_mpi_ierr)
	      CALL mpas_mpi_check('distributed restart marker publication')
#endif

	   END SUBROUTINE ncio_complete_distributed_write

	   !---------------------------------------------------------
	   SUBROUTINE ncio_reset_distributed_validation()

	      IF (allocated(distributed_vector_maps)) deallocate(distributed_vector_maps)
	      IF (allocated(validated_distributed_markers)) deallocate(validated_distributed_markers)
	      IF (allocated(distributed_marker_ranks)) deallocate(distributed_marker_ranks)

	   END SUBROUTINE ncio_reset_distributed_validation

	   !---------------------------------------------------------
	   SUBROUTINE prepare_distributed_vector_map(filename, pixelset, iblk, jblk, map_index)

	   USE MOD_NetCDFSerial
	   USE MOD_Block, only: get_filename_block
	   USE MOD_Pixelset
	   USE MOD_MPAS_MPI, only: CoLM_stop
	   IMPLICIT NONE

	   character(len=*), intent(in) :: filename
	   type(pixelset_type), intent(in) :: pixelset
	   integer, intent(in) :: iblk, jblk
	   integer, intent(out) :: map_index

#ifdef MPAS_EMBEDDED_COLM
	   type(distributed_vector_map_type), allocatable :: maps_new(:)
	   character(len=256) :: fileblock, source_file
	   character(len=512) :: complete_marker
	   integer :: istt, iend, nlocal, marker_ranks, irank, isource, ilocal
	   integer :: imap, nmaps, nused, iused
	   integer*8, allocatable :: file_eindex(:)
	   integer, allocatable :: file_ipxstt(:), file_ipxend(:), file_settyp(:)
	   integer, allocatable :: rank_lengths(:)
	   logical, allocatable :: used_rank(:)
	   logical :: source_exists

	      istt = pixelset%vecgs%vstt(iblk,jblk)
	      iend = pixelset%vecgs%vend(iblk,jblk)
	      nlocal = iend - istt + 1
	      IF (nlocal < 1 .or. istt < 1 .or. iend > pixelset%nset) THEN
	         CALL CoLM_stop('Invalid local pixelset block while preparing distributed restart mapping.')
	      ENDIF
	      IF (.not. allocated(pixelset%eindex) .or. .not. allocated(pixelset%ipxstt) .or. &
	          .not. allocated(pixelset%ipxend) .or. .not. allocated(pixelset%settyp)) THEN
	         CALL CoLM_stop('Current MPAS embedded CoLM pixelset lacks restart identity metadata.')
	      ENDIF
	      IF (size(pixelset%eindex) /= pixelset%nset .or. size(pixelset%ipxstt) /= pixelset%nset .or. &
	          size(pixelset%ipxend) /= pixelset%nset .or. size(pixelset%settyp) /= pixelset%nset) THEN
	         CALL CoLM_stop('Current MPAS embedded CoLM pixelset restart identity lengths are inconsistent.')
	      ENDIF
	      CALL validate_identity_vector(pixelset%eindex(istt:iend), pixelset%ipxstt(istt:iend), &
	                                    pixelset%ipxend(istt:iend), pixelset%settyp(istt:iend), &
	                                    'current MPAS decomposition')

	      IF (allocated(distributed_vector_maps)) THEN
	         DO imap = 1, size(distributed_vector_maps)
	            IF (trim(distributed_vector_maps(imap)%filename) /= trim(filename) .or. &
	                distributed_vector_maps(imap)%iblk /= iblk .or. &
	                distributed_vector_maps(imap)%jblk /= jblk .or. &
	                distributed_vector_maps(imap)%istt /= istt .or. &
	                distributed_vector_maps(imap)%iend /= iend) CYCLE
	            IF (all(distributed_vector_maps(imap)%eindex == pixelset%eindex(istt:iend)) .and. &
	                all(distributed_vector_maps(imap)%ipxstt == pixelset%ipxstt(istt:iend)) .and. &
	                all(distributed_vector_maps(imap)%ipxend == pixelset%ipxend(istt:iend)) .and. &
	                all(distributed_vector_maps(imap)%settyp == pixelset%settyp(istt:iend))) THEN
	               map_index = imap
	               RETURN
	            ENDIF
	         ENDDO
	      ENDIF

	      complete_marker = trim(filename)//'.mpas_complete'
	      CALL read_distributed_marker(complete_marker, marker_ranks)
	      CALL get_filename_block(filename, iblk, jblk, fileblock)
	      allocate(rank_lengths(marker_ranks), used_rank(marker_ranks))
	      rank_lengths = 0
	      used_rank = .false.

	      IF (allocated(distributed_vector_maps)) THEN
	         nmaps = size(distributed_vector_maps)
	      ELSE
	         nmaps = 0
	      ENDIF
	      allocate(maps_new(nmaps + 1))
	      IF (nmaps > 0) maps_new(1:nmaps) = distributed_vector_maps
	      map_index = nmaps + 1
	      maps_new(map_index)%filename = filename
	      maps_new(map_index)%iblk = iblk
	      maps_new(map_index)%jblk = jblk
	      maps_new(map_index)%istt = istt
	      maps_new(map_index)%iend = iend
	      maps_new(map_index)%eindex = pixelset%eindex(istt:iend)
	      maps_new(map_index)%ipxstt = pixelset%ipxstt(istt:iend)
	      maps_new(map_index)%ipxend = pixelset%ipxend(istt:iend)
	      maps_new(map_index)%settyp = pixelset%settyp(istt:iend)
	      allocate(maps_new(map_index)%source_rank(nlocal))
	      allocate(maps_new(map_index)%source_pos(nlocal))
	      maps_new(map_index)%source_rank = -1
	      maps_new(map_index)%source_pos = -1

	      DO irank = 0, marker_ranks - 1
	         CALL add_rank_suffix(fileblock, irank, source_file)
	         inquire(file=trim(source_file), exist=source_exists)
	         IF (.not. source_exists) CYCLE
	         IF (.not. ncio_var_exist(source_file, 'mpas_eindex', readflag=.false.) .or. &
	             .not. ncio_var_exist(source_file, 'mpas_ipxstt', readflag=.false.) .or. &
	             .not. ncio_var_exist(source_file, 'mpas_ipxend', readflag=.false.) .or. &
	             .not. ncio_var_exist(source_file, 'mpas_settyp', readflag=.false.)) THEN
	            CALL CoLM_stop('Distributed MPAS-CoLM restart shard lacks pixelset identity metadata: '// &
	                           trim(source_file))
	         ENDIF
	         CALL ncio_read_serial(source_file, 'mpas_eindex', file_eindex)
	         CALL ncio_read_serial(source_file, 'mpas_ipxstt', file_ipxstt)
	         CALL ncio_read_serial(source_file, 'mpas_ipxend', file_ipxend)
	         CALL ncio_read_serial(source_file, 'mpas_settyp', file_settyp)
	         IF (size(file_eindex) < 1 .or. size(file_ipxstt) /= size(file_eindex) .or. &
	             size(file_ipxend) /= size(file_eindex) .or. size(file_settyp) /= size(file_eindex)) THEN
	            CALL CoLM_stop('Distributed MPAS-CoLM restart shard has inconsistent identity metadata: '// &
	                           trim(source_file))
	         ENDIF
	         CALL validate_identity_vector(file_eindex, file_ipxstt, file_ipxend, file_settyp, source_file)
	         rank_lengths(irank + 1) = size(file_eindex)

	         DO isource = 1, size(file_eindex)
	            ilocal = find_identity(pixelset%eindex(istt:iend), pixelset%ipxstt(istt:iend), &
	                                   pixelset%ipxend(istt:iend), pixelset%settyp(istt:iend), &
	                                   file_eindex(isource), file_ipxstt(isource), &
	                                   file_ipxend(isource), file_settyp(isource))
	            IF (ilocal > 0) THEN
	               IF (maps_new(map_index)%source_rank(ilocal) >= 0) THEN
	                  CALL CoLM_stop('Duplicate pixelset identity across distributed MPAS-CoLM restart shards.')
	               ENDIF
	               maps_new(map_index)%source_rank(ilocal) = irank
	               maps_new(map_index)%source_pos(ilocal) = isource
	               used_rank(irank + 1) = .true.
	            ENDIF
	         ENDDO

	         deallocate(file_eindex, file_ipxstt, file_ipxend, file_settyp)
	      ENDDO

	      IF (any(maps_new(map_index)%source_rank < 0)) THEN
	         ilocal = minloc(maps_new(map_index)%source_rank, dim=1)
	         write(*,'(A,I0,A,I0,A,I0,A,I0)') 'Missing distributed restart identity: element=', &
	            maps_new(map_index)%eindex(ilocal), ', ipxstt=', maps_new(map_index)%ipxstt(ilocal), &
	            ', ipxend=', maps_new(map_index)%ipxend(ilocal), ', settyp=', maps_new(map_index)%settyp(ilocal)
	         CALL CoLM_stop('Distributed MPAS-CoLM restart does not cover the current pixelset decomposition.')
	      ENDIF

	      nused = count(used_rank)
	      allocate(maps_new(map_index)%source_ranks(nused))
	      allocate(maps_new(map_index)%source_lengths(nused))
	      iused = 0
	      DO irank = 0, marker_ranks - 1
	         IF (.not. used_rank(irank + 1)) CYCLE
	         iused = iused + 1
	         maps_new(map_index)%source_ranks(iused) = irank
	         maps_new(map_index)%source_lengths(iused) = rank_lengths(irank + 1)
	      ENDDO

	      deallocate(rank_lengths, used_rank)
	      CALL move_alloc(maps_new, distributed_vector_maps)
#else
	      map_index = 0
#endif

	   END SUBROUTINE prepare_distributed_vector_map

	   !---------------------------------------------------------
	   SUBROUTINE validate_identity_vector(eindex, ipxstt, ipxend, settyp, source_name)

	   USE MOD_MPAS_MPI, only: CoLM_stop
	   IMPLICIT NONE

	   integer*8, intent(in) :: eindex(:)
	   integer, intent(in) :: ipxstt(:), ipxend(:), settyp(:)
	   character(len=*), intent(in) :: source_name
	   integer :: istart, iend, i, j

	      IF (size(ipxstt) /= size(eindex) .or. size(ipxend) /= size(eindex) .or. &
	          size(settyp) /= size(eindex)) THEN
	         CALL CoLM_stop('Inconsistent pixelset identity vector lengths in '//trim(source_name)//'.')
	      ENDIF
	      IF (size(eindex) > 1) THEN
	         IF (any(eindex(2:) < eindex(:size(eindex)-1))) THEN
	            CALL CoLM_stop('Pixelset identities are not ordered by element in '//trim(source_name)//'.')
	         ENDIF
	      ENDIF

	      istart = 1
	      DO WHILE (istart <= size(eindex))
	         iend = istart
	         DO WHILE (iend < size(eindex))
	            IF (eindex(iend + 1) /= eindex(istart)) EXIT
	            iend = iend + 1
	         ENDDO
	         DO i = istart, iend - 1
	            DO j = i + 1, iend
	               IF (ipxstt(i) == ipxstt(j) .and. ipxend(i) == ipxend(j) .and. &
	                   settyp(i) == settyp(j)) THEN
	                  CALL CoLM_stop('Duplicate pixelset identity in '//trim(source_name)//'.')
	               ENDIF
	            ENDDO
	         ENDDO
	         istart = iend + 1
	      ENDDO

	   END SUBROUTINE validate_identity_vector

	   !---------------------------------------------------------
	   INTEGER FUNCTION find_identity(eindex, ipxstt, ipxend, settyp, target_eindex, &
	                                  target_ipxstt, target_ipxend, target_settyp)

	   IMPLICIT NONE

	   integer*8, intent(in) :: eindex(:), target_eindex
	   integer, intent(in) :: ipxstt(:), ipxend(:), settyp(:)
	   integer, intent(in) :: target_ipxstt, target_ipxend, target_settyp
	   integer :: lo, hi, mid, ifirst, ilast, i

	      find_identity = 0
	      lo = 1
	      hi = size(eindex)
	      DO WHILE (lo <= hi)
	         mid = lo + (hi - lo) / 2
	         IF (eindex(mid) < target_eindex) THEN
	            lo = mid + 1
	         ELSEIF (eindex(mid) > target_eindex) THEN
	            hi = mid - 1
	         ELSE
	            ifirst = mid
	            ilast = mid
	            DO WHILE (ifirst > 1)
	               IF (eindex(ifirst - 1) /= target_eindex) EXIT
	               ifirst = ifirst - 1
	            ENDDO
	            DO WHILE (ilast < size(eindex))
	               IF (eindex(ilast + 1) /= target_eindex) EXIT
	               ilast = ilast + 1
	            ENDDO
	            DO i = ifirst, ilast
	               IF (ipxstt(i) == target_ipxstt .and. ipxend(i) == target_ipxend .and. &
	                   settyp(i) == target_settyp) THEN
	                  find_identity = i
	                  RETURN
	               ENDIF
	            ENDDO
	            RETURN
	         ENDIF
	      ENDDO

	   END FUNCTION find_identity

	   !---------------------------------------------------------
	   SUBROUTINE get_distributed_source(map_index, source_index, source_file, source_rank, source_length)

	   USE MOD_Block, only: get_filename_block
	   USE MOD_MPAS_MPI, only: CoLM_stop
	   IMPLICIT NONE

	   integer, intent(in) :: map_index, source_index
	   character(len=*), intent(out) :: source_file
	   integer, intent(out) :: source_rank, source_length

#ifdef MPAS_EMBEDDED_COLM
	   character(len=256) :: fileblock

	      IF (.not. allocated(distributed_vector_maps) .or. map_index < 1 .or. &
	          map_index > size(distributed_vector_maps)) THEN
	         CALL CoLM_stop('Invalid distributed MPAS-CoLM restart map index.')
	      ENDIF
	      IF (source_index < 1 .or. source_index > size(distributed_vector_maps(map_index)%source_ranks)) THEN
	         CALL CoLM_stop('Invalid distributed MPAS-CoLM restart source index.')
	      ENDIF
	      source_rank = distributed_vector_maps(map_index)%source_ranks(source_index)
	      source_length = distributed_vector_maps(map_index)%source_lengths(source_index)
	      CALL get_filename_block(distributed_vector_maps(map_index)%filename, &
	                              distributed_vector_maps(map_index)%iblk, &
	                              distributed_vector_maps(map_index)%jblk, fileblock)
	      CALL add_rank_suffix(fileblock, source_rank, source_file)
#else
	      source_file = ''
	      source_rank = -1
	      source_length = 0
#endif

	   END SUBROUTINE get_distributed_source

#ifdef MPAS_EMBEDDED_COLM
	   !---------------------------------------------------------
	   SUBROUTINE read_distributed_int32_1d(filename, dataname, pixelset, iblk, jblk, rdata, found, defval)

	   USE MOD_NetCDFSerial
	   USE MOD_Pixelset, only: pixelset_type
	   IMPLICIT NONE

	   character(len=*), intent(in) :: filename, dataname
	   type(pixelset_type), intent(in) :: pixelset
	   integer, intent(in) :: iblk, jblk
	   integer, allocatable, intent(out) :: rdata(:)
	   logical, intent(out) :: found
	   integer, intent(in), optional :: defval

	   integer :: map_index, isource, source_rank, source_length, ilocal
	   character(len=256) :: source_file
	   integer, allocatable :: source_data(:)

	      CALL prepare_distributed_vector_map(filename, pixelset, iblk, jblk, map_index)
	      allocate(rdata(size(distributed_vector_maps(map_index)%source_rank)))
	      IF (present(defval)) rdata = defval
	      found = .false.
	      DO isource = 1, size(distributed_vector_maps(map_index)%source_ranks)
	         CALL get_distributed_source(map_index, isource, source_file, source_rank, source_length)
	         IF (ncio_var_exist(source_file, dataname, readflag=.false.)) THEN
	            CALL ncio_read_serial(source_file, dataname, source_data)
	            IF (size(source_data) /= source_length) THEN
	               CALL CoLM_stop('Distributed vector length mismatch for '//trim(dataname)//' in '//trim(source_file)//'.')
	            ENDIF
	            found = .true.
	            DO ilocal = 1, size(rdata)
	               IF (distributed_vector_maps(map_index)%source_rank(ilocal) == source_rank) THEN
	                  rdata(ilocal) = source_data(distributed_vector_maps(map_index)%source_pos(ilocal))
	               ENDIF
	            ENDDO
	            deallocate(source_data)
	         ELSEIF (.not. present(defval)) THEN
	            CALL ncio_vector_stop_missing_block(filename, dataname, source_file)
	         ENDIF
	      ENDDO

	   END SUBROUTINE read_distributed_int32_1d

	   !---------------------------------------------------------
	   SUBROUTINE read_distributed_int64_1d(filename, dataname, pixelset, iblk, jblk, rdata, found, defval)

	   USE MOD_NetCDFSerial
	   USE MOD_Pixelset, only: pixelset_type
	   IMPLICIT NONE

	   character(len=*), intent(in) :: filename, dataname
	   type(pixelset_type), intent(in) :: pixelset
	   integer, intent(in) :: iblk, jblk
	   integer*8, allocatable, intent(out) :: rdata(:)
	   logical, intent(out) :: found
	   integer*8, intent(in), optional :: defval

	   integer :: map_index, isource, source_rank, source_length, ilocal
	   character(len=256) :: source_file
	   integer*8, allocatable :: source_data(:)

	      CALL prepare_distributed_vector_map(filename, pixelset, iblk, jblk, map_index)
	      allocate(rdata(size(distributed_vector_maps(map_index)%source_rank)))
	      IF (present(defval)) rdata = defval
	      found = .false.
	      DO isource = 1, size(distributed_vector_maps(map_index)%source_ranks)
	         CALL get_distributed_source(map_index, isource, source_file, source_rank, source_length)
	         IF (ncio_var_exist(source_file, dataname, readflag=.false.)) THEN
	            CALL ncio_read_serial(source_file, dataname, source_data)
	            IF (size(source_data) /= source_length) THEN
	               CALL CoLM_stop('Distributed vector length mismatch for '//trim(dataname)//' in '//trim(source_file)//'.')
	            ENDIF
	            found = .true.
	            DO ilocal = 1, size(rdata)
	               IF (distributed_vector_maps(map_index)%source_rank(ilocal) == source_rank) THEN
	                  rdata(ilocal) = source_data(distributed_vector_maps(map_index)%source_pos(ilocal))
	               ENDIF
	            ENDDO
	            deallocate(source_data)
	         ELSEIF (.not. present(defval)) THEN
	            CALL ncio_vector_stop_missing_block(filename, dataname, source_file)
	         ENDIF
	      ENDDO

	   END SUBROUTINE read_distributed_int64_1d

	   !---------------------------------------------------------
	   SUBROUTINE read_distributed_int8_1d(filename, dataname, pixelset, iblk, jblk, rdata, found, defval)

	   USE MOD_NetCDFSerial
	   USE MOD_Pixelset, only: pixelset_type
	   IMPLICIT NONE

	   character(len=*), intent(in) :: filename, dataname
	   type(pixelset_type), intent(in) :: pixelset
	   integer, intent(in) :: iblk, jblk
	   integer(1), allocatable, intent(out) :: rdata(:)
	   logical, intent(out) :: found
	   logical, intent(in), optional :: defval

	   integer :: map_index, isource, source_rank, source_length, ilocal
	   character(len=256) :: source_file
	   integer(1), allocatable :: source_data(:)

	      CALL prepare_distributed_vector_map(filename, pixelset, iblk, jblk, map_index)
	      allocate(rdata(size(distributed_vector_maps(map_index)%source_rank)))
	      IF (present(defval)) THEN
	         IF (defval) THEN
	            rdata = 1
	         ELSE
	            rdata = 0
	         ENDIF
	      ENDIF
	      found = .false.
	      DO isource = 1, size(distributed_vector_maps(map_index)%source_ranks)
	         CALL get_distributed_source(map_index, isource, source_file, source_rank, source_length)
	         IF (ncio_var_exist(source_file, dataname, readflag=.false.)) THEN
	            CALL ncio_read_serial(source_file, dataname, source_data)
	            IF (size(source_data) /= source_length) THEN
	               CALL CoLM_stop('Distributed vector length mismatch for '//trim(dataname)//' in '//trim(source_file)//'.')
	            ENDIF
	            IF (any((source_data /= int(0,kind=kind(source_data))) .and. &
	                    (source_data /= int(1,kind=kind(source_data))))) THEN
	               CALL CoLM_stop('Invalid logical vector value in '//trim(source_file)//'.')
	            ENDIF
	            found = .true.
	            DO ilocal = 1, size(rdata)
	               IF (distributed_vector_maps(map_index)%source_rank(ilocal) == source_rank) THEN
	                  rdata(ilocal) = source_data(distributed_vector_maps(map_index)%source_pos(ilocal))
	               ENDIF
	            ENDDO
	            deallocate(source_data)
	         ELSEIF (.not. present(defval)) THEN
	            CALL ncio_vector_stop_missing_block(filename, dataname, source_file)
	         ENDIF
	      ENDDO

	   END SUBROUTINE read_distributed_int8_1d

	   !---------------------------------------------------------
	   SUBROUTINE read_distributed_real8_1d(filename, dataname, pixelset, iblk, jblk, rdata, found, defval)

	   USE MOD_Precision
	   USE MOD_NetCDFSerial
	   USE MOD_Pixelset, only: pixelset_type
	   IMPLICIT NONE

	   character(len=*), intent(in) :: filename, dataname
	   type(pixelset_type), intent(in) :: pixelset
	   integer, intent(in) :: iblk, jblk
	   real(r8), allocatable, intent(out) :: rdata(:)
	   logical, intent(out) :: found
	   real(r8), intent(in), optional :: defval

	   integer :: map_index, isource, source_rank, source_length, ilocal
	   character(len=256) :: source_file
	   real(r8), allocatable :: source_data(:)

	      CALL prepare_distributed_vector_map(filename, pixelset, iblk, jblk, map_index)
	      allocate(rdata(size(distributed_vector_maps(map_index)%source_rank)))
	      IF (present(defval)) rdata = defval
	      found = .false.
	      DO isource = 1, size(distributed_vector_maps(map_index)%source_ranks)
	         CALL get_distributed_source(map_index, isource, source_file, source_rank, source_length)
	         IF (ncio_var_exist(source_file, dataname, readflag=.false.)) THEN
	            CALL ncio_read_serial(source_file, dataname, source_data)
	            IF (size(source_data) /= source_length) THEN
	               CALL CoLM_stop('Distributed vector length mismatch for '//trim(dataname)//' in '//trim(source_file)//'.')
	            ENDIF
	            found = .true.
	            DO ilocal = 1, size(rdata)
	               IF (distributed_vector_maps(map_index)%source_rank(ilocal) == source_rank) THEN
	                  rdata(ilocal) = source_data(distributed_vector_maps(map_index)%source_pos(ilocal))
	               ENDIF
	            ENDDO
	            deallocate(source_data)
	         ELSEIF (.not. present(defval)) THEN
	            CALL ncio_vector_stop_missing_block(filename, dataname, source_file)
	         ENDIF
	      ENDDO

	   END SUBROUTINE read_distributed_real8_1d

	   !---------------------------------------------------------
	   SUBROUTINE read_distributed_real8_2d(filename, dataname, ndim1, pixelset, iblk, jblk, rdata, found, defval)

	   USE MOD_Precision
	   USE MOD_NetCDFSerial
	   USE MOD_Pixelset, only: pixelset_type
	   IMPLICIT NONE

	   character(len=*), intent(in) :: filename, dataname
	   integer, intent(in) :: ndim1, iblk, jblk
	   type(pixelset_type), intent(in) :: pixelset
	   real(r8), allocatable, intent(out) :: rdata(:,:)
	   logical, intent(out) :: found
	   real(r8), intent(in), optional :: defval

	   integer :: map_index, isource, source_rank, source_length, ilocal
	   character(len=256) :: source_file
	   real(r8), allocatable :: source_data(:,:)

	      CALL prepare_distributed_vector_map(filename, pixelset, iblk, jblk, map_index)
	      allocate(rdata(ndim1,size(distributed_vector_maps(map_index)%source_rank)))
	      IF (present(defval)) rdata = defval
	      found = .false.
	      DO isource = 1, size(distributed_vector_maps(map_index)%source_ranks)
	         CALL get_distributed_source(map_index, isource, source_file, source_rank, source_length)
	         IF (ncio_var_exist(source_file, dataname, readflag=.false.)) THEN
	            CALL ncio_read_serial(source_file, dataname, source_data)
	            CALL validate_vector_fixed_dimensions((/size(source_data,1)/), (/ndim1/), source_file, dataname)
	            IF (size(source_data,2) /= source_length) THEN
	               CALL CoLM_stop('Distributed vector length mismatch for '//trim(dataname)//' in '//trim(source_file)//'.')
	            ENDIF
	            found = .true.
	            DO ilocal = 1, size(rdata,2)
	               IF (distributed_vector_maps(map_index)%source_rank(ilocal) == source_rank) THEN
	                  rdata(:,ilocal) = source_data(:,distributed_vector_maps(map_index)%source_pos(ilocal))
	               ENDIF
	            ENDDO
	            deallocate(source_data)
	         ELSEIF (.not. present(defval)) THEN
	            CALL ncio_vector_stop_missing_block(filename, dataname, source_file)
	         ENDIF
	      ENDDO

	   END SUBROUTINE read_distributed_real8_2d

	   !---------------------------------------------------------
	   SUBROUTINE read_distributed_real8_3d(filename, dataname, ndim1, ndim2, pixelset, iblk, jblk, &
	                                       rdata, found, defval)

	   USE MOD_Precision
	   USE MOD_NetCDFSerial
	   USE MOD_Pixelset, only: pixelset_type
	   IMPLICIT NONE

	   character(len=*), intent(in) :: filename, dataname
	   integer, intent(in) :: ndim1, ndim2, iblk, jblk
	   type(pixelset_type), intent(in) :: pixelset
	   real(r8), allocatable, intent(out) :: rdata(:,:,:)
	   logical, intent(out) :: found
	   real(r8), intent(in), optional :: defval

	   integer :: map_index, isource, source_rank, source_length, ilocal
	   character(len=256) :: source_file
	   real(r8), allocatable :: source_data(:,:,:)

	      CALL prepare_distributed_vector_map(filename, pixelset, iblk, jblk, map_index)
	      allocate(rdata(ndim1,ndim2,size(distributed_vector_maps(map_index)%source_rank)))
	      IF (present(defval)) rdata = defval
	      found = .false.
	      DO isource = 1, size(distributed_vector_maps(map_index)%source_ranks)
	         CALL get_distributed_source(map_index, isource, source_file, source_rank, source_length)
	         IF (ncio_var_exist(source_file, dataname, readflag=.false.)) THEN
	            CALL ncio_read_serial(source_file, dataname, source_data)
	            CALL validate_vector_fixed_dimensions((/size(source_data,1),size(source_data,2)/), &
	                                                   (/ndim1,ndim2/), source_file, dataname)
	            IF (size(source_data,3) /= source_length) THEN
	               CALL CoLM_stop('Distributed vector length mismatch for '//trim(dataname)//' in '//trim(source_file)//'.')
	            ENDIF
	            found = .true.
	            DO ilocal = 1, size(rdata,3)
	               IF (distributed_vector_maps(map_index)%source_rank(ilocal) == source_rank) THEN
	                  rdata(:,:,ilocal) = source_data(:,:,distributed_vector_maps(map_index)%source_pos(ilocal))
	               ENDIF
	            ENDDO
	            deallocate(source_data)
	         ELSEIF (.not. present(defval)) THEN
	            CALL ncio_vector_stop_missing_block(filename, dataname, source_file)
	         ENDIF
	      ENDDO

	   END SUBROUTINE read_distributed_real8_3d

	   !---------------------------------------------------------
	   SUBROUTINE read_distributed_real8_4d(filename, dataname, ndim1, ndim2, ndim3, pixelset, iblk, jblk, &
	                                       rdata, found, defval)

	   USE MOD_Precision
	   USE MOD_NetCDFSerial
	   USE MOD_Pixelset, only: pixelset_type
	   IMPLICIT NONE

	   character(len=*), intent(in) :: filename, dataname
	   integer, intent(in) :: ndim1, ndim2, ndim3, iblk, jblk
	   type(pixelset_type), intent(in) :: pixelset
	   real(r8), allocatable, intent(out) :: rdata(:,:,:,:)
	   logical, intent(out) :: found
	   real(r8), intent(in), optional :: defval

	   integer :: map_index, isource, source_rank, source_length, ilocal
	   character(len=256) :: source_file
	   real(r8), allocatable :: source_data(:,:,:,:)

	      CALL prepare_distributed_vector_map(filename, pixelset, iblk, jblk, map_index)
	      allocate(rdata(ndim1,ndim2,ndim3,size(distributed_vector_maps(map_index)%source_rank)))
	      IF (present(defval)) rdata = defval
	      found = .false.
	      DO isource = 1, size(distributed_vector_maps(map_index)%source_ranks)
	         CALL get_distributed_source(map_index, isource, source_file, source_rank, source_length)
	         IF (ncio_var_exist(source_file, dataname, readflag=.false.)) THEN
	            CALL ncio_read_serial(source_file, dataname, source_data)
	            CALL validate_vector_fixed_dimensions((/size(source_data,1),size(source_data,2), &
	                                                     size(source_data,3)/), &
	                                                   (/ndim1,ndim2,ndim3/), source_file, dataname)
	            IF (size(source_data,4) /= source_length) THEN
	               CALL CoLM_stop('Distributed vector length mismatch for '//trim(dataname)//' in '//trim(source_file)//'.')
	            ENDIF
	            found = .true.
	            DO ilocal = 1, size(rdata,4)
	               IF (distributed_vector_maps(map_index)%source_rank(ilocal) == source_rank) THEN
	                  rdata(:,:,:,ilocal) = source_data(:,:,:,distributed_vector_maps(map_index)%source_pos(ilocal))
	               ENDIF
	            ENDDO
	            deallocate(source_data)
	         ELSEIF (.not. present(defval)) THEN
	            CALL ncio_vector_stop_missing_block(filename, dataname, source_file)
	         ENDIF
	      ENDDO

	   END SUBROUTINE read_distributed_real8_4d

	   !---------------------------------------------------------
	   SUBROUTINE read_distributed_real8_5d(filename, dataname, ndim1, ndim2, ndim3, ndim4, pixelset, iblk, jblk, &
	                                       rdata, found, defval)

	   USE MOD_Precision
	   USE MOD_NetCDFSerial
	   USE MOD_Pixelset, only: pixelset_type
	   IMPLICIT NONE

	   character(len=*), intent(in) :: filename, dataname
	   integer, intent(in) :: ndim1, ndim2, ndim3, ndim4, iblk, jblk
	   type(pixelset_type), intent(in) :: pixelset
	   real(r8), allocatable, intent(out) :: rdata(:,:,:,:,:)
	   logical, intent(out) :: found
	   real(r8), intent(in), optional :: defval

	   integer :: map_index, isource, source_rank, source_length, ilocal
	   character(len=256) :: source_file
	   real(r8), allocatable :: source_data(:,:,:,:,:)

	      CALL prepare_distributed_vector_map(filename, pixelset, iblk, jblk, map_index)
	      allocate(rdata(ndim1,ndim2,ndim3,ndim4,size(distributed_vector_maps(map_index)%source_rank)))
	      IF (present(defval)) rdata = defval
	      found = .false.
	      DO isource = 1, size(distributed_vector_maps(map_index)%source_ranks)
	         CALL get_distributed_source(map_index, isource, source_file, source_rank, source_length)
	         IF (ncio_var_exist(source_file, dataname, readflag=.false.)) THEN
	            CALL ncio_read_serial(source_file, dataname, source_data)
	            CALL validate_vector_fixed_dimensions((/size(source_data,1),size(source_data,2), &
	                                                     size(source_data,3),size(source_data,4)/), &
	                                                   (/ndim1,ndim2,ndim3,ndim4/), source_file, dataname)
	            IF (size(source_data,5) /= source_length) THEN
	               CALL CoLM_stop('Distributed vector length mismatch for '//trim(dataname)//' in '//trim(source_file)//'.')
	            ENDIF
	            found = .true.
	            DO ilocal = 1, size(rdata,5)
	               IF (distributed_vector_maps(map_index)%source_rank(ilocal) == source_rank) THEN
	                  rdata(:,:,:,:,ilocal) = source_data(:,:,:,:,distributed_vector_maps(map_index)%source_pos(ilocal))
	               ENDIF
	            ENDDO
	            deallocate(source_data)
	         ELSEIF (.not. present(defval)) THEN
	            CALL ncio_vector_stop_missing_block(filename, dataname, source_file)
	         ENDIF
	      ENDDO

	   END SUBROUTINE read_distributed_real8_5d
#endif

	   !---------------------------------------------------------
	   SUBROUTINE ncio_vector_stop_missing_block (filename, dataname, fileblock)

	   USE MOD_MPAS_MPI, only: CoLM_stop
	   IMPLICIT NONE

	   character(len=*), intent(in) :: filename
	   character(len=*), intent(in) :: dataname
	   character(len=*), intent(in) :: fileblock

	      write(*,*) 'Error : required vector data '//trim(dataname) &
	         //' in '//trim(filename)//' is missing from block file '//trim(fileblock)//'.'
	      CALL CoLM_stop ()

	   END SUBROUTINE ncio_vector_stop_missing_block

	   !---------------------------------------------------------
	   SUBROUTINE validate_vector_block_data(pixelset, istt, iend, source_length, &
	         use_srcpos, fileblock, dataname)

	   USE MOD_Pixelset, only: pixelset_type
	   USE MOD_MPAS_MPI, only: CoLM_stop
	   IMPLICIT NONE

	   type(pixelset_type), intent(in) :: pixelset
	   integer, intent(in) :: istt, iend, source_length
	   logical, intent(in) :: use_srcpos
	   character(len=*), intent(in) :: fileblock, dataname

	   integer :: expected_length

	      IF (istt < 1 .or. iend < istt .or. iend > pixelset%nset) THEN
	         CALL CoLM_stop('Invalid local vector bounds while reading '//trim(dataname)// &
	                        ' from '//trim(fileblock)//'.')
	      ENDIF

	      expected_length = iend - istt + 1
	      IF (use_srcpos) THEN
	         IF (.not. allocated(pixelset%srcpos)) THEN
	            CALL CoLM_stop('Missing source-position map while reading '//trim(dataname)// &
	                           ' from '//trim(fileblock)//'.')
	         ENDIF
	         IF (size(pixelset%srcpos) /= pixelset%nset) THEN
	            CALL CoLM_stop('Invalid source-position map length while reading '//trim(dataname)// &
	                           ' from '//trim(fileblock)//'.')
	         ENDIF
	         IF (source_length < 1 .or. minval(pixelset%srcpos(istt:iend)) < 1 .or. &
	             maxval(pixelset%srcpos(istt:iend)) > source_length) THEN
	            CALL CoLM_stop('Source-position map exceeds vector data while reading '//trim(dataname)// &
	                           ' from '//trim(fileblock)//'.')
	         ENDIF
	      ELSEIF (source_length /= expected_length) THEN
	         CALL CoLM_stop('Rank-local vector length mismatch while reading '//trim(dataname)// &
	                        ' from '//trim(fileblock)//'.')
	      ENDIF

	   END SUBROUTINE validate_vector_block_data

	   !---------------------------------------------------------
	   SUBROUTINE validate_vector_fixed_dimensions(actual, expected, fileblock, dataname)

	   USE MOD_MPAS_MPI, only: CoLM_stop
	   IMPLICIT NONE

	   integer, intent(in) :: actual(:), expected(:)
	   character(len=*), intent(in) :: fileblock, dataname

	      IF (size(actual) /= size(expected)) THEN
	         CALL CoLM_stop('Non-vector rank mismatch while reading '//trim(dataname)// &
	                        ' from '//trim(fileblock)//'.')
	      ENDIF
	      IF (any(actual /= expected)) THEN
	         CALL CoLM_stop('Non-vector dimension mismatch while reading '//trim(dataname)// &
	                        ' from '//trim(fileblock)//'.')
	      ENDIF

	   END SUBROUTINE validate_vector_fixed_dimensions

	   !---------------------------------------------------------
	   SUBROUTINE ncio_read_vector_int32_1d ( &
         filename, dataname, pixelset, rdata, defval)

   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   type(pixelset_type), intent(in) :: pixelset

   integer, allocatable, intent(inout) :: rdata (:)
   integer, intent(in), optional :: defval

   ! Local variables
	   integer :: iblkgrp, iblk, jblk, istt, iend, iset
	   character(len=256) :: fileblock
	   integer, allocatable :: sbuff(:)
			   logical :: any_data_exists, block_has_data, use_srcpos, distributed_restart

      IF (.true.) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (pixelset%nset))
         ENDIF
      ENDIF

      any_data_exists = .false.

      IF (.true.) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            use_srcpos = allocated(pixelset%srcpos)
            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, use_srcpos = use_srcpos, &
                                            distributed_restart = distributed_restart)

	            IF (distributed_restart) THEN
	               CALL read_distributed_vector_block(filename, dataname, pixelset, iblk, jblk, &
	                                                  sbuff, block_has_data, defval)
	            ELSE
	               allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))
	               block_has_data = .false.
			               IF (ncio_var_exist(fileblock,dataname,readflag=.false.)) THEN
			                  CALL ncio_read_serial (fileblock, dataname, sbuff)
			                  block_has_data = .true.
			               ELSEIF (present(defval)) THEN
			                  sbuff(:) = defval
			               ELSEIF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
			                  CALL ncio_vector_stop_missing_block (filename, dataname, fileblock)
			               ENDIF
	            ENDIF
	            any_data_exists = any_data_exists .or. block_has_data

	            istt = pixelset%vecgs%vstt(iblk,jblk)
	            iend = pixelset%vecgs%vend(iblk,jblk)
	            IF (block_has_data .and. .not. distributed_restart) &
	               CALL validate_vector_block_data(pixelset, istt, iend, size(sbuff), &
	                                                               use_srcpos, fileblock, dataname)
			            IF (use_srcpos .and. block_has_data .and. .not. distributed_restart) THEN
	               DO iset = istt, iend
	                  rdata(iset) = sbuff(pixelset%srcpos(iset))
	               ENDDO
	            ELSE
	               rdata(istt:iend) = sbuff
	            ENDIF

            deallocate (sbuff)

         ENDDO

         IF (pixelset%nset > 0 .and. .not. any_data_exists) THEN
            IF (ncio_vector_report_missing(.not. present(defval))) THEN
               IF (.not. present(defval)) THEN
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found.'
                  CALL CoLM_stop ()
               ELSE
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found, default value is used.'
               ENDIF
            ENDIF
         ENDIF
      ENDIF


   END SUBROUTINE ncio_read_vector_int32_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_int64_1d ( &
         filename, dataname, pixelset, rdata, defval)

   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   type(pixelset_type), intent(in) :: pixelset

   integer*8, allocatable, intent(inout) :: rdata (:)
   integer*8, intent(in), optional :: defval

   ! Local variables
	   integer :: iblkgrp, iblk, jblk, istt, iend, iset
	   character(len=256) :: fileblock
	   integer*8, allocatable :: sbuff(:)
			   logical :: any_data_exists, block_has_data, use_srcpos, distributed_restart

      IF (.true.) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (pixelset%nset))
         ENDIF
      ENDIF

      any_data_exists = .false.

      IF (.true.) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            use_srcpos = allocated(pixelset%srcpos)
            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, use_srcpos = use_srcpos, &
                                            distributed_restart = distributed_restart)

	            IF (distributed_restart) THEN
	               CALL read_distributed_vector_block(filename, dataname, pixelset, iblk, jblk, &
	                                                  sbuff, block_has_data, defval)
	            ELSE
	               allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))
	               block_has_data = .false.
			               IF (ncio_var_exist(fileblock,dataname,readflag=.false.)) THEN
			                  CALL ncio_read_serial (fileblock, dataname, sbuff)
			                  block_has_data = .true.
			               ELSEIF (present(defval)) THEN
			                  sbuff(:) = defval
			               ELSEIF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
			                  CALL ncio_vector_stop_missing_block (filename, dataname, fileblock)
			               ENDIF
	            ENDIF
	            any_data_exists = any_data_exists .or. block_has_data

	            istt = pixelset%vecgs%vstt(iblk,jblk)
	            iend = pixelset%vecgs%vend(iblk,jblk)
	            IF (block_has_data .and. .not. distributed_restart) &
	               CALL validate_vector_block_data(pixelset, istt, iend, size(sbuff), &
	                                                               use_srcpos, fileblock, dataname)
			            IF (use_srcpos .and. block_has_data .and. .not. distributed_restart) THEN
	               DO iset = istt, iend
	                  rdata(iset) = sbuff(pixelset%srcpos(iset))
	               ENDDO
	            ELSE
	               rdata(istt:iend) = sbuff
	            ENDIF

            deallocate (sbuff)

         ENDDO

         IF (pixelset%nset > 0 .and. .not. any_data_exists) THEN
            IF (ncio_vector_report_missing(.not. present(defval))) THEN
               IF (.not. present(defval)) THEN
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found.'
                  CALL CoLM_stop ()
               ELSE
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found, default value is used.'
               ENDIF
            ENDIF
         ENDIF
      ENDIF


   END SUBROUTINE ncio_read_vector_int64_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_logical_1d (filename, dataname, pixelset, rdata, &
         defval)

   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   type(pixelset_type), intent(in) :: pixelset

   logical, allocatable, intent(inout) :: rdata (:)
   logical, intent(in), optional :: defval

   ! Local variables
	   integer :: iblkgrp, iblk, jblk, istt, iend, iset
	   character(len=256) :: fileblock
	   integer(1), allocatable :: sbuff(:)
			   logical :: any_data_exists, block_has_data, use_srcpos, distributed_restart

      IF (.true.) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (pixelset%nset))
         ENDIF
      ENDIF

      any_data_exists = .false.

      IF (.true.) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            use_srcpos = allocated(pixelset%srcpos)
            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, use_srcpos = use_srcpos, &
                                            distributed_restart = distributed_restart)

	            IF (distributed_restart) THEN
	               CALL read_distributed_vector_block(filename, dataname, pixelset, iblk, jblk, &
	                                                  sbuff, block_has_data, defval)
	            ELSE
	               allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))
	               block_has_data = .false.
		               IF (ncio_var_exist(fileblock,dataname,readflag=.false.)) THEN
		                  CALL ncio_read_serial (fileblock, dataname, sbuff)
		                  block_has_data = .true.
			               ELSEIF (present(defval)) THEN
		                  IF (defval) THEN
		                     sbuff(:) = 1
		                  ELSE
		                     sbuff(:) = 0
		                  ENDIF
		               ELSEIF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
		                  CALL ncio_vector_stop_missing_block (filename, dataname, fileblock)
		               ENDIF
	            ENDIF
	            any_data_exists = any_data_exists .or. block_has_data

	            istt = pixelset%vecgs%vstt(iblk,jblk)
	            iend = pixelset%vecgs%vend(iblk,jblk)
	            IF (block_has_data .and. .not. distributed_restart) &
	               CALL validate_vector_block_data(pixelset, istt, iend, size(sbuff), &
	                                                               use_srcpos, fileblock, dataname)
			            IF (use_srcpos .and. block_has_data .and. .not. distributed_restart) THEN
	               DO iset = istt, iend
	                  rdata(iset) = (sbuff(pixelset%srcpos(iset)) == int(1, kind=kind(sbuff)))
	               ENDDO
	            ELSE
	               rdata(istt:iend) = (sbuff == int(1, kind=kind(sbuff)))
	            ENDIF

            deallocate (sbuff)

         ENDDO

         IF (pixelset%nset > 0 .and. .not. any_data_exists) THEN
            IF (ncio_vector_report_missing(.not. present(defval))) THEN
               IF (.not. present(defval)) THEN
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found.'
                  CALL CoLM_stop ()
               ELSE
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found, default value is used.'
               ENDIF
            ENDIF
         ENDIF
      ENDIF


   END SUBROUTINE ncio_read_vector_logical_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_real8_1d (filename, dataname, pixelset, rdata, &
         defval)

   USE MOD_Precision
   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   type(pixelset_type), intent(in) :: pixelset

   real(r8), allocatable, intent(inout) :: rdata (:)
   real(r8), intent(in), optional :: defval

   ! Local variables
	   integer :: iblkgrp, iblk, jblk, istt, iend, iset
	   character(len=256) :: fileblock
	   real(r8), allocatable :: sbuff(:)
			   logical :: any_data_exists, block_has_data, use_srcpos, distributed_restart

      IF (.true.) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (pixelset%nset))
         ENDIF
      ENDIF

      any_data_exists = .false.

      IF (.true.) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            use_srcpos = allocated(pixelset%srcpos)
            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, use_srcpos = use_srcpos, &
                                            distributed_restart = distributed_restart)

	            IF (distributed_restart) THEN
	               CALL read_distributed_vector_block(filename, dataname, pixelset, iblk, jblk, &
	                                                  sbuff, block_has_data, defval)
	            ELSE
	               allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))
	               block_has_data = .false.
			               IF (ncio_var_exist(fileblock,dataname,readflag=.false.)) THEN
			                  CALL ncio_read_serial (fileblock, dataname, sbuff)
			                  block_has_data = .true.
			               ELSEIF (present(defval)) THEN
			                  sbuff(:) = defval
			               ELSEIF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
			                  CALL ncio_vector_stop_missing_block (filename, dataname, fileblock)
			               ENDIF
	            ENDIF
	            any_data_exists = any_data_exists .or. block_has_data

	            istt = pixelset%vecgs%vstt(iblk,jblk)
	            iend = pixelset%vecgs%vend(iblk,jblk)
	            IF (block_has_data .and. .not. distributed_restart) &
	               CALL validate_vector_block_data(pixelset, istt, iend, size(sbuff), &
	                                                               use_srcpos, fileblock, dataname)
			            IF (use_srcpos .and. block_has_data .and. .not. distributed_restart) THEN
	               DO iset = istt, iend
	                  rdata(iset) = sbuff(pixelset%srcpos(iset))
	               ENDDO
	            ELSE
	               rdata(istt:iend) = sbuff
	            ENDIF

            deallocate (sbuff)

         ENDDO

         IF (pixelset%nset > 0 .and. .not. any_data_exists) THEN
            IF (ncio_vector_report_missing(.not. present(defval))) THEN
               IF (.not. present(defval)) THEN
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found.'
                  CALL CoLM_stop ()
               ELSE
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found, default value is used.'
               ENDIF
            ENDIF
         ENDIF
      ENDIF


   END SUBROUTINE ncio_read_vector_real8_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_real8_2d ( &
         filename, dataname, ndim1, pixelset, rdata, defval)

   USE MOD_Precision
   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   integer, intent(in) :: ndim1
   type(pixelset_type), intent(in) :: pixelset

   real(r8), allocatable, intent(inout) :: rdata (:,:)
   real(r8), intent(in), optional :: defval

   ! Local variables
	   integer :: iblkgrp, iblk, jblk, istt, iend, iset
	   character(len=256) :: fileblock
	   real(r8), allocatable :: sbuff(:,:)
		   logical :: any_data_exists, block_has_data, use_srcpos, distributed_restart

      IF (.true.) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (ndim1, pixelset%nset))
         ENDIF
      ENDIF

      any_data_exists = .false.

      IF (.true.) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            use_srcpos = allocated(pixelset%srcpos)
            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, use_srcpos = use_srcpos, &
                                            distributed_restart = distributed_restart)

	            IF (distributed_restart) THEN
	               CALL read_distributed_vector_block(filename, dataname, ndim1, pixelset, iblk, jblk, &
	                                                  sbuff, block_has_data, defval)
	            ELSE
	               allocate (sbuff (ndim1, pixelset%vecgs%vlen(iblk,jblk)))
	               block_has_data = .false.
		               IF (ncio_var_exist(fileblock,dataname,readflag=.false.)) THEN
		                  CALL ncio_read_serial (fileblock, dataname, sbuff)
		                  block_has_data = .true.
		               ELSEIF (present(defval)) THEN
		                  sbuff(:,:) = defval
		               ELSEIF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
		                  CALL ncio_vector_stop_missing_block (filename, dataname, fileblock)
		               ENDIF
	            ENDIF
	            any_data_exists = any_data_exists .or. block_has_data

	            istt = pixelset%vecgs%vstt(iblk,jblk)
	            iend = pixelset%vecgs%vend(iblk,jblk)
	            IF (block_has_data .and. .not. distributed_restart) THEN
	               CALL validate_vector_fixed_dimensions((/size(sbuff,1)/), (/ndim1/), fileblock, dataname)
	               CALL validate_vector_block_data(pixelset, istt, iend, size(sbuff,2), &
	                                               use_srcpos, fileblock, dataname)
	            ENDIF
		            IF (use_srcpos .and. block_has_data .and. .not. distributed_restart) THEN
	               DO iset = istt, iend
	                  rdata(:,iset) = sbuff(:,pixelset%srcpos(iset))
	               ENDDO
	            ELSE
	               rdata(:,istt:iend) = sbuff
	            ENDIF

            deallocate (sbuff)

         ENDDO

         IF (pixelset%nset > 0 .and. .not. any_data_exists) THEN
            IF (ncio_vector_report_missing(.not. present(defval))) THEN
               IF (.not. present(defval)) THEN
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found.'
                  CALL CoLM_stop ()
               ELSE
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found, default value is used.'
               ENDIF
            ENDIF
         ENDIF
      ENDIF


   END SUBROUTINE ncio_read_vector_real8_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_real8_3d ( &
         filename, dataname, ndim1, ndim2, pixelset, rdata, defval)

   USE MOD_Precision
   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   integer, intent(in) :: ndim1, ndim2
   type(pixelset_type), intent(in) :: pixelset

   real(r8), allocatable, intent(inout) :: rdata (:,:,:)
   real(r8), intent(in), optional :: defval

   ! Local variables
	   integer :: iblkgrp, iblk, jblk, istt, iend, iset
	   character(len=256) :: fileblock
	   real(r8), allocatable :: sbuff(:,:,:)
		   logical :: any_data_exists, block_has_data, use_srcpos, distributed_restart

      IF (.true.) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (ndim1,ndim2, pixelset%nset))
         ENDIF
      ENDIF

      any_data_exists = .false.

      IF (.true.) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            use_srcpos = allocated(pixelset%srcpos)
            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, use_srcpos = use_srcpos, &
                                            distributed_restart = distributed_restart)

	            IF (distributed_restart) THEN
	               CALL read_distributed_vector_block(filename, dataname, ndim1, ndim2, pixelset, iblk, jblk, &
	                                                  sbuff, block_has_data, defval)
	            ELSE
	               allocate (sbuff (ndim1,ndim2, pixelset%vecgs%vlen(iblk,jblk)))
	               block_has_data = .false.
		               IF (ncio_var_exist(fileblock,dataname,readflag=.false.)) THEN
		                  CALL ncio_read_serial (fileblock, dataname, sbuff)
		                  block_has_data = .true.
		               ELSEIF (present(defval)) THEN
		                  sbuff(:,:,:) = defval
		               ELSEIF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
		                  CALL ncio_vector_stop_missing_block (filename, dataname, fileblock)
		               ENDIF
	            ENDIF
	            any_data_exists = any_data_exists .or. block_has_data

	            istt = pixelset%vecgs%vstt(iblk,jblk)
	            iend = pixelset%vecgs%vend(iblk,jblk)
	            IF (block_has_data .and. .not. distributed_restart) THEN
	               CALL validate_vector_fixed_dimensions((/size(sbuff,1),size(sbuff,2)/), &
	                                                     (/ndim1,ndim2/), fileblock, dataname)
	               CALL validate_vector_block_data(pixelset, istt, iend, size(sbuff,3), &
	                                               use_srcpos, fileblock, dataname)
	            ENDIF
		            IF (use_srcpos .and. block_has_data .and. .not. distributed_restart) THEN
	               DO iset = istt, iend
	                  rdata(:,:,iset) = sbuff(:,:,pixelset%srcpos(iset))
	               ENDDO
	            ELSE
	               rdata(:,:,istt:iend) = sbuff
	            ENDIF

            deallocate (sbuff)

         ENDDO

         IF (pixelset%nset > 0 .and. .not. any_data_exists) THEN
            IF (ncio_vector_report_missing(.not. present(defval))) THEN
               IF (.not. present(defval)) THEN
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found.'
                  CALL CoLM_stop ()
               ELSE
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found, default value is used.'
               ENDIF
            ENDIF
         ENDIF
      ENDIF


   END SUBROUTINE ncio_read_vector_real8_3d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_real8_4d ( &
         filename, dataname, ndim1, ndim2, ndim3, pixelset, rdata, defval)

   USE MOD_Precision
   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   integer, intent(in) :: ndim1, ndim2, ndim3
   type(pixelset_type), intent(in) :: pixelset

   real(r8), allocatable, intent(inout) :: rdata (:,:,:,:)
   real(r8), intent(in), optional :: defval

   ! Local variables
	   integer :: iblkgrp, iblk, jblk, istt, iend, iset
	   character(len=256) :: fileblock
	   real(r8), allocatable :: sbuff(:,:,:,:)
		   logical :: any_data_exists, block_has_data, use_srcpos, distributed_restart

      IF (.true.) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (ndim1,ndim2,ndim3, pixelset%nset))
         ENDIF
      ENDIF

      any_data_exists = .false.

      IF (.true.) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            use_srcpos = allocated(pixelset%srcpos)
            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, use_srcpos = use_srcpos, &
                                            distributed_restart = distributed_restart)

	            IF (distributed_restart) THEN
	               CALL read_distributed_vector_block(filename, dataname, ndim1, ndim2, ndim3, &
	                                                  pixelset, iblk, jblk, sbuff, block_has_data, defval)
	            ELSE
	               allocate (sbuff (ndim1,ndim2,ndim3, pixelset%vecgs%vlen(iblk,jblk)))
	               block_has_data = .false.
		               IF (ncio_var_exist(fileblock,dataname,readflag=.false.)) THEN
		                  CALL ncio_read_serial (fileblock, dataname, sbuff)
		                  block_has_data = .true.
		               ELSEIF (present(defval)) THEN
		                  sbuff(:,:,:,:) = defval
		               ELSEIF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
		                  CALL ncio_vector_stop_missing_block (filename, dataname, fileblock)
		               ENDIF
	            ENDIF
	            any_data_exists = any_data_exists .or. block_has_data

	            istt = pixelset%vecgs%vstt(iblk,jblk)
	            iend = pixelset%vecgs%vend(iblk,jblk)
	            IF (block_has_data .and. .not. distributed_restart) THEN
	               CALL validate_vector_fixed_dimensions((/size(sbuff,1),size(sbuff,2),size(sbuff,3)/), &
	                                                     (/ndim1,ndim2,ndim3/), fileblock, dataname)
	               CALL validate_vector_block_data(pixelset, istt, iend, size(sbuff,4), &
	                                               use_srcpos, fileblock, dataname)
	            ENDIF
		            IF (use_srcpos .and. block_has_data .and. .not. distributed_restart) THEN
	               DO iset = istt, iend
	                  rdata(:,:,:,iset) = sbuff(:,:,:,pixelset%srcpos(iset))
	               ENDDO
	            ELSE
	               rdata(:,:,:,istt:iend) = sbuff
	            ENDIF

            deallocate (sbuff)

         ENDDO

         IF (pixelset%nset > 0 .and. .not. any_data_exists) THEN
            IF (ncio_vector_report_missing(.not. present(defval))) THEN
               IF (.not. present(defval)) THEN
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found.'
                  CALL CoLM_stop ()
               ELSE
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found, default value is used.'
               ENDIF
            ENDIF
         ENDIF
      ENDIF


   END SUBROUTINE ncio_read_vector_real8_4d


   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_real8_5d ( &
         filename, dataname, ndim1, ndim2, ndim3, ndim4, pixelset, rdata, defval)

   USE MOD_Precision
   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   integer, intent(in) :: ndim1, ndim2, ndim3, ndim4
   type(pixelset_type), intent(in) :: pixelset

   real(r8), allocatable, intent(inout) :: rdata (:,:,:,:,:)
   real(r8), intent(in), optional :: defval

   ! Local variables
	   integer :: iblkgrp, iblk, jblk, istt, iend, iset
	   character(len=256) :: fileblock
	   real(r8), allocatable :: sbuff(:,:,:,:,:)
		   logical :: any_data_exists, block_has_data, use_srcpos, distributed_restart

      IF (.true.) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (ndim1,ndim2,ndim3,ndim4, pixelset%nset))
         ENDIF
      ENDIF

      any_data_exists = .false.

      IF (.true.) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            use_srcpos = allocated(pixelset%srcpos)
            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, use_srcpos = use_srcpos, &
                                            distributed_restart = distributed_restart)

	            IF (distributed_restart) THEN
	               CALL read_distributed_vector_block(filename, dataname, ndim1, ndim2, ndim3, ndim4, &
	                                                  pixelset, iblk, jblk, sbuff, block_has_data, defval)
	            ELSE
	               allocate (sbuff (ndim1,ndim2,ndim3,ndim4, pixelset%vecgs%vlen(iblk,jblk)))
	               block_has_data = .false.
		               IF (ncio_var_exist(fileblock,dataname,readflag=.false.)) THEN
		                  CALL ncio_read_serial (fileblock, dataname, sbuff)
		                  block_has_data = .true.
		               ELSEIF (present(defval)) THEN
		                  sbuff(:,:,:,:,:) = defval
		               ELSEIF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
		                  CALL ncio_vector_stop_missing_block (filename, dataname, fileblock)
		               ENDIF
	            ENDIF
	            any_data_exists = any_data_exists .or. block_has_data

	            istt = pixelset%vecgs%vstt(iblk,jblk)
	            iend = pixelset%vecgs%vend(iblk,jblk)
	            IF (block_has_data .and. .not. distributed_restart) THEN
	               CALL validate_vector_fixed_dimensions((/size(sbuff,1),size(sbuff,2),size(sbuff,3), &
	                                                       size(sbuff,4)/), &
	                                                     (/ndim1,ndim2,ndim3,ndim4/), fileblock, dataname)
	               CALL validate_vector_block_data(pixelset, istt, iend, size(sbuff,5), &
	                                               use_srcpos, fileblock, dataname)
	            ENDIF
		            IF (use_srcpos .and. block_has_data .and. .not. distributed_restart) THEN
	               DO iset = istt, iend
	                  rdata(:,:,:,:,iset) = sbuff(:,:,:,:,pixelset%srcpos(iset))
	               ENDDO
	            ELSE
	               rdata(:,:,:,:,istt:iend) = sbuff
	            ENDIF

            deallocate (sbuff)

         ENDDO

         IF (pixelset%nset > 0 .and. .not. any_data_exists) THEN
            IF (ncio_vector_report_missing(.not. present(defval))) THEN
               IF (.not. present(defval)) THEN
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found.'
                  CALL CoLM_stop ()
               ELSE
                  write(*,*) 'Warning : restart data '//trim(dataname) &
                     //' in '//trim(filename)//' not found, default value is used.'
               ENDIF
            ENDIF
         ENDIF
      ENDIF


   END SUBROUTINE ncio_read_vector_real8_5d


   !---------------------------------------------------------
   SUBROUTINE ncio_create_file_vector (filename, pixelset)

   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   type(pixelset_type), intent(in) :: pixelset

   ! Local variables
   integer :: iblkgrp, iblk, jblk
#ifdef MPAS_EMBEDDED_COLM
   integer :: istt, iend
#endif
   character(len=256) :: fileblock

      IF (.true.) THEN
         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, for_write = .true.)
            CALL ncio_create_file (fileblock)
#ifdef MPAS_EMBEDDED_COLM
	            istt = pixelset%vecgs%vstt(iblk,jblk)
	            iend = pixelset%vecgs%vend(iblk,jblk)
	            IF (iend >= istt) THEN
               IF (.not. allocated(pixelset%eindex) .or. .not. allocated(pixelset%ipxstt) .or. &
                   .not. allocated(pixelset%ipxend) .or. .not. allocated(pixelset%settyp)) THEN
                  write(*,'(A)') 'Error: MPAS embedded CoLM cannot write vector state without pixelset identity metadata.'
	                  CALL CoLM_stop()
	               ENDIF
	               CALL ncio_define_dimension(fileblock, 'mpas_local_set', iend - istt + 1)
	               CALL ncio_write_serial(fileblock, 'mpas_eindex', pixelset%eindex(istt:iend), 'mpas_local_set')
               CALL ncio_write_serial(fileblock, 'mpas_ipxstt', pixelset%ipxstt(istt:iend), 'mpas_local_set')
               CALL ncio_write_serial(fileblock, 'mpas_ipxend', pixelset%ipxend(istt:iend), 'mpas_local_set')
               CALL ncio_write_serial(fileblock, 'mpas_settyp', pixelset%settyp(istt:iend), 'mpas_local_set')
            ENDIF
#endif

         ENDDO
      ENDIF

   END SUBROUTINE ncio_create_file_vector

   !---------------------------------------------------------
   SUBROUTINE ncio_define_dimension_vector (filename, pixelset, dimname, dimlen)

   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*),    intent(in) :: filename
   type(pixelset_type), intent(in) :: pixelset
   character(len=*), intent(in)  :: dimname
   integer, intent(in), optional :: dimlen

   ! Local variables
   integer :: iblkgrp, iblk, jblk
   character(len=256) :: fileblock
   logical :: fexists

      IF (.true.) THEN
         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, for_write = .true.)
            inquire (file=trim(fileblock), exist=fexists)
            IF (.not. fexists) THEN
               CALL ncio_create_file (fileblock)
            ENDIF

            IF (present(dimlen)) THEN
               CALL ncio_define_dimension (fileblock, trim(dimname), dimlen)
            ELSE
               CALL ncio_define_dimension (fileblock, trim(dimname), &
                  pixelset%vecgs%vlen(iblk,jblk))
            ENDIF

         ENDDO
      ENDIF

   END SUBROUTINE ncio_define_dimension_vector

	   !---------------------------------------------------------
	   SUBROUTINE validate_vector_write_layout(pixelset, vector_size, dataname, actual_dims, expected_dims)

	   USE MOD_Pixelset, only: pixelset_type
	   USE MOD_MPAS_MPI, only: CoLM_stop
	   IMPLICIT NONE

	   type(pixelset_type), intent(in) :: pixelset
	   integer, intent(in) :: vector_size
	   character(len=*), intent(in) :: dataname
	   integer, intent(in), optional :: actual_dims(:)
	   integer, intent(in), optional :: expected_dims(:)

	   integer :: iblkgrp
	   integer :: iblk
	   integer :: jblk
	   integer :: istt
	   integer :: iend
	   integer :: vlen
	   logical, allocatable :: covered(:)

	      IF (pixelset%nset < 0 .or. vector_size /= pixelset%nset) THEN
	         CALL CoLM_stop('Invalid vector length while writing '//trim(dataname)//'.')
	      ENDIF
	      IF (present(actual_dims) .neqv. present(expected_dims)) THEN
	         CALL CoLM_stop('Incomplete fixed-dimension metadata while writing '//trim(dataname)//'.')
	      ENDIF
	      IF (present(actual_dims) .and. present(expected_dims)) THEN
	         IF (size(actual_dims) /= size(expected_dims) .or. any(actual_dims /= expected_dims) .or. &
	             any(actual_dims <= 0)) THEN
	            CALL CoLM_stop('Fixed dimensions do not match the data while writing '//trim(dataname)//'.')
	         ENDIF
	      ENDIF

	      IF (pixelset%nblkgrp < 0) THEN
	         CALL CoLM_stop('Invalid block count while writing '//trim(dataname)//'.')
	      ENDIF
	      IF (pixelset%nset == 0) THEN
	         IF (pixelset%nblkgrp /= 0) CALL CoLM_stop('Empty vector has nonempty block layout while writing '//trim(dataname)//'.')
	         RETURN
	      ENDIF
	      IF (pixelset%nblkgrp < 1 .or. .not. allocated(pixelset%xblkgrp) .or. &
	          .not. allocated(pixelset%yblkgrp) .or. .not. allocated(pixelset%vecgs%vlen) .or. &
	          .not. allocated(pixelset%vecgs%vstt) .or. .not. allocated(pixelset%vecgs%vend)) THEN
	         CALL CoLM_stop('Missing vector block layout while writing '//trim(dataname)//'.')
	      ENDIF
	      IF (size(pixelset%xblkgrp) /= pixelset%nblkgrp .or. &
	          size(pixelset%yblkgrp) /= pixelset%nblkgrp) THEN
	         CALL CoLM_stop('Inconsistent vector block list while writing '//trim(dataname)//'.')
	      ENDIF
	      IF (any(shape(pixelset%vecgs%vlen) /= shape(pixelset%vecgs%vstt)) .or. &
	          any(shape(pixelset%vecgs%vlen) /= shape(pixelset%vecgs%vend))) THEN
	         CALL CoLM_stop('Inconsistent vector block arrays while writing '//trim(dataname)//'.')
	      ENDIF

	      allocate(covered(pixelset%nset))
	      covered = .false.
	      DO iblkgrp = 1, pixelset%nblkgrp
	         iblk = pixelset%xblkgrp(iblkgrp)
	         jblk = pixelset%yblkgrp(iblkgrp)
	         IF (iblk < lbound(pixelset%vecgs%vlen,1) .or. iblk > ubound(pixelset%vecgs%vlen,1) .or. &
	             jblk < lbound(pixelset%vecgs%vlen,2) .or. jblk > ubound(pixelset%vecgs%vlen,2)) THEN
	            CALL CoLM_stop('Vector block index is out of range while writing '//trim(dataname)//'.')
	         ENDIF
	         vlen = pixelset%vecgs%vlen(iblk,jblk)
	         istt = pixelset%vecgs%vstt(iblk,jblk)
	         iend = pixelset%vecgs%vend(iblk,jblk)
	         IF (vlen < 1 .or. istt < 1 .or. iend < istt .or. iend > pixelset%nset .or. &
	             iend - istt + 1 /= vlen) THEN
	            CALL CoLM_stop('Invalid vector block bounds while writing '//trim(dataname)//'.')
	         ENDIF
	         IF (any(covered(istt:iend))) THEN
	            CALL CoLM_stop('Overlapping vector blocks while writing '//trim(dataname)//'.')
	         ENDIF
	         covered(istt:iend) = .true.
	      ENDDO
	      IF (.not. all(covered)) CALL CoLM_stop('Vector block layout has gaps while writing '//trim(dataname)//'.')
	      deallocate(covered)

	   END SUBROUTINE validate_vector_write_layout

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_int32_1d ( &
         filename, dataname, dimname, pixelset, wdata, compress_level)

   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: dimname
   type(pixelset_type), intent(in) :: pixelset
   integer, intent(in) :: wdata (:)

   integer, intent(in), optional :: compress_level

   ! Local variables
   integer :: iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: fileblock
   integer, allocatable :: rbuff(:)

      IF (.true.) THEN
	         CALL validate_vector_write_layout(pixelset, size(wdata), dataname)

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff (pixelset%vecgs%vlen(iblk,jblk)))
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rbuff = wdata(istt:iend)

            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, for_write = .true.)

            IF (present(compress_level)) THEN
               CALL ncio_write_serial (fileblock, dataname, rbuff, dimname, &
                  compress = compress_level)
            ELSE
               CALL ncio_write_serial (fileblock, dataname, rbuff, dimname)
            ENDIF

            deallocate (rbuff)

         ENDDO

      ENDIF


   END SUBROUTINE ncio_write_vector_int32_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_logical_1d ( &
         filename, dataname, dimname, pixelset, wdata, compress_level)

   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: dimname
   type(pixelset_type), intent(in) :: pixelset
   logical, intent(in) :: wdata (:)

   integer, intent(in), optional :: compress_level

   ! Local variables
   integer :: iblkgrp, iblk, jblk, istt, iend, i
   character(len=256) :: fileblock
   integer(1), allocatable :: rbuff(:)

      IF (.true.) THEN
	         CALL validate_vector_write_layout(pixelset, size(wdata), dataname)

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff (pixelset%vecgs%vlen(iblk,jblk)))
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            DO i = istt, iend
               IF(wdata(i))THEN
                  rbuff(i-istt+1) = 1
               ELSE
                  rbuff(i-istt+1) = 0
               ENDIF
            ENDDO

            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, for_write = .true.)

            IF (present(compress_level)) THEN
               CALL ncio_write_serial (fileblock, dataname, rbuff, dimname, &
                  compress = compress_level)
            ELSE
               CALL ncio_write_serial (fileblock, dataname, rbuff, dimname)
            ENDIF

            deallocate (rbuff)

         ENDDO

      ENDIF


   END SUBROUTINE ncio_write_vector_logical_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_int32_3d ( &
         filename, dataname, dim1name, ndim1, dim2name, ndim2, &
         dim3name, pixelset, wdata, compress_level)

   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: dim1name, dim2name, dim3name
   type(pixelset_type), intent(in) :: pixelset
   integer, intent(in) :: ndim1, ndim2
   integer, intent(in) :: wdata (:,:,:)

   integer, intent(in), optional :: compress_level

   ! Local variables
   integer :: iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: fileblock
   integer, allocatable :: rbuff(:,:,:)

      IF (.true.) THEN
	         CALL validate_vector_write_layout(pixelset, size(wdata,3), dataname, &
	                                           (/size(wdata,1), size(wdata,2)/), (/ndim1, ndim2/))

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff (ndim1,ndim2,pixelset%vecgs%vlen(iblk,jblk)))
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rbuff = wdata(:,:,istt:iend)

            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, for_write = .true.)

            IF (present(compress_level)) THEN
               CALL ncio_write_serial (fileblock, dataname, rbuff, &
                  dim1name, dim2name, dim3name, compress = compress_level)
            ELSE
               CALL ncio_write_serial (fileblock, dataname, rbuff, &
                  dim1name, dim2name, dim3name)
            ENDIF

            deallocate (rbuff)

         ENDDO

      ENDIF


   END SUBROUTINE ncio_write_vector_int32_3d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_int64_1d ( &
         filename, dataname, dimname, pixelset, wdata, compress_level)

   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: dimname
   type(pixelset_type), intent(in) :: pixelset
   integer*8, intent(in) :: wdata (:)

   integer, intent(in), optional :: compress_level

   ! Local variables
   integer :: iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: fileblock
   integer*8, allocatable :: rbuff(:)

      IF (.true.) THEN
	         CALL validate_vector_write_layout(pixelset, size(wdata), dataname)

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff (pixelset%vecgs%vlen(iblk,jblk)))
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rbuff = wdata(istt:iend)

            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, for_write = .true.)

            IF (present(compress_level)) THEN
               CALL ncio_write_serial (fileblock, dataname, rbuff, dimname, &
                  compress = compress_level)
            ELSE
               CALL ncio_write_serial (fileblock, dataname, rbuff, dimname)
            ENDIF

            deallocate (rbuff)

         ENDDO

      ENDIF


   END SUBROUTINE ncio_write_vector_int64_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_real8_1d ( &
         filename, dataname, dimname, pixelset, wdata, compress_level)

   USE MOD_Precision
   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: dimname
   type(pixelset_type), intent(in) :: pixelset
   real(r8), intent(in) :: wdata (:)

   integer, intent(in), optional :: compress_level

   ! Local variables
   integer :: iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: fileblock
   real(r8), allocatable :: rbuff(:)

      IF (.true.) THEN
	         CALL validate_vector_write_layout(pixelset, size(wdata), dataname)

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff (pixelset%vecgs%vlen(iblk,jblk)))
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rbuff = wdata(istt:iend)

            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, for_write = .true.)
            IF (present(compress_level)) THEN
               CALL ncio_write_serial (fileblock, dataname, rbuff, &
                  dimname, compress = compress_level)
            ELSE
               CALL ncio_write_serial (fileblock, dataname, rbuff, dimname)
            ENDIF

            deallocate (rbuff)

         ENDDO

      ENDIF


   END SUBROUTINE ncio_write_vector_real8_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_real8_2d ( &
         filename, dataname, dim1name, ndim1, &
         dim2name, pixelset, wdata, compress_level)

   USE MOD_Precision
   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: dim1name, dim2name
   integer,  intent(in) :: ndim1
   type(pixelset_type), intent(in) :: pixelset
   real(r8), intent(in) :: wdata (:,:)

   integer,  intent(in), optional :: compress_level

   ! Local variables
   integer :: iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: fileblock
   real(r8), allocatable :: rbuff(:,:)

      IF (.true.) THEN
	         CALL validate_vector_write_layout(pixelset, size(wdata,2), dataname, &
	                                           (/size(wdata,1)/), (/ndim1/))

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff (ndim1, pixelset%vecgs%vlen(iblk,jblk)))
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rbuff = wdata(:,istt:iend)

            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, for_write = .true.)

            IF (present(compress_level)) THEN
               CALL ncio_write_serial (fileblock, dataname, rbuff, &
                  dim1name, dim2name, compress = compress_level)
            ELSE
               CALL ncio_write_serial (fileblock, dataname, rbuff, &
                  dim1name, dim2name)
            ENDIF

            deallocate (rbuff)

         ENDDO

      ENDIF


   END SUBROUTINE ncio_write_vector_real8_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_real8_3d ( &
         filename, dataname, dim1name, ndim1, dim2name, ndim2, &
         dim3name, pixelset, wdata, compress_level)

   USE MOD_Precision
   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: dim1name, dim2name, dim3name
   type(pixelset_type), intent(in) :: pixelset
   integer,  intent(in) :: ndim1, ndim2
   real(r8), intent(in) :: wdata (:,:,:)

   integer,  intent(in), optional :: compress_level

   ! Local variables
   integer :: iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: fileblock
   real(r8), allocatable :: rbuff(:,:,:)

      IF (.true.) THEN
	         CALL validate_vector_write_layout(pixelset, size(wdata,3), dataname, &
	                                           (/size(wdata,1), size(wdata,2)/), (/ndim1, ndim2/))

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff (ndim1, ndim2, pixelset%vecgs%vlen(iblk,jblk)))
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rbuff = wdata(:,:,istt:iend)

            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, for_write = .true.)
            IF (present(compress_level)) THEN
               CALL ncio_write_serial (fileblock, dataname, rbuff, &
                  dim1name, dim2name, dim3name, compress = compress_level)
            ELSE
               CALL ncio_write_serial (fileblock, dataname, rbuff, &
                  dim1name, dim2name, dim3name)
            ENDIF

            deallocate (rbuff)

         ENDDO

      ENDIF


   END SUBROUTINE ncio_write_vector_real8_3d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_real8_4d ( &
         filename, dataname, dim1name, ndim1, dim2name, ndim2, dim3name, ndim3, &
         dim4name, pixelset, wdata, compress_level)

   USE MOD_Precision
   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: dim1name, dim2name, dim3name, dim4name
   integer,  intent(in) :: ndim1, ndim2, ndim3
   type(pixelset_type), intent(in) :: pixelset
   real(r8), intent(in) :: wdata (:,:,:,:)

   integer,  intent(in), optional :: compress_level

   ! Local variables
   integer :: iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: fileblock
   real(r8), allocatable :: rbuff(:,:,:,:)

      IF (.true.) THEN
	         CALL validate_vector_write_layout(pixelset, size(wdata,4), dataname, &
	                                           (/size(wdata,1), size(wdata,2), size(wdata,3)/), &
	                                           (/ndim1, ndim2, ndim3/))

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff (ndim1, ndim2, ndim3, pixelset%vecgs%vlen(iblk,jblk)))
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rbuff = wdata(:,:,:,istt:iend)

            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, for_write = .true.)
            IF (present(compress_level)) THEN
               CALL ncio_write_serial (fileblock, dataname, rbuff, &
                  dim1name, dim2name, dim3name, dim4name, compress = compress_level)
            ELSE
               CALL ncio_write_serial (fileblock, dataname, rbuff, &
                  dim1name, dim2name, dim3name, dim4name)
            ENDIF

            deallocate (rbuff)

         ENDDO

      ENDIF


   END SUBROUTINE ncio_write_vector_real8_4d


   !------------------------------------------------
   SUBROUTINE ncio_write_vector_real8_5d ( &
         filename, dataname, dim1name, ndim1, dim2name, ndim2, &
         dim3name, ndim3, dim4name, ndim4, dim5name, pixelset, wdata, compress_level)

   USE MOD_Precision
   USE MOD_NetCDFSerial
   USE MOD_MPAS_MPI
   USE MOD_Block
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: dim1name, dim2name, dim3name, dim4name, dim5name
   type(pixelset_type), intent(in) :: pixelset
   integer,  intent(in) :: ndim1, ndim2, ndim3, ndim4
   real(r8), intent(in) :: wdata (:,:,:,:,:)

   integer,  intent(in), optional :: compress_level

   ! Local variables
   integer :: iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: fileblock
   real(r8), allocatable :: rbuff(:,:,:,:,:)

      IF (.true.) THEN
	         CALL validate_vector_write_layout(pixelset, size(wdata,5), dataname, &
	                                           (/size(wdata,1), size(wdata,2), size(wdata,3), size(wdata,4)/), &
	                                           (/ndim1, ndim2, ndim3, ndim4/))

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff (ndim1, ndim2, ndim3, ndim4, pixelset%vecgs%vlen(iblk,jblk)))
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rbuff = wdata(:,:,:,:,istt:iend)

            CALL get_filename_vector_block (filename, iblk, jblk, fileblock, for_write = .true.)
            IF (present(compress_level)) THEN
               CALL ncio_write_serial (fileblock, dataname, rbuff, &
                  dim1name, dim2name, dim3name, dim4name, dim5name, compress = compress_level)
            ELSE
               CALL ncio_write_serial (fileblock, dataname, rbuff, &
                  dim1name, dim2name, dim3name, dim4name, dim5name)
            ENDIF

            deallocate (rbuff)

         ENDDO

      ENDIF


   END SUBROUTINE ncio_write_vector_real8_5d
   !------------------------------------------------

   LOGICAL FUNCTION ncio_vector_report_missing(mandatory)

   USE MOD_MPAS_MPI, only: mpas_rank, mpas_root
   IMPLICIT NONE

   logical, intent(in) :: mandatory

#ifdef MPAS_EMBEDDED_COLM
      ncio_vector_report_missing = mandatory .or. (mpas_rank == mpas_root)
#else
      ncio_vector_report_missing = (mpas_rank == mpas_root)
#endif

   END FUNCTION ncio_vector_report_missing


END MODULE MOD_NetCDFVector
