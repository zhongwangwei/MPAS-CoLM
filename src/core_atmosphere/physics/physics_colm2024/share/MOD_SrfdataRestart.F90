#include <define.h>

MODULE MOD_SrfdataRestart
!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    This module includes subroutines to read/write data of mesh and pixelsets.
!
!  Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   IMPLICIT NONE

   ! ----- subroutines -----
   PUBLIC :: mesh_save_to_file
   PUBLIC :: mesh_load_from_file

   PUBLIC :: pixelset_save_to_file
   PUBLIC :: pixelset_load_from_file

CONTAINS

   ! -----------------------
   SUBROUTINE mesh_save_to_file (dir_landdata, lc_year)

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_Mesh
   USE MOD_Block
   USE MOD_Utils
   IMPLICIT NONE

   character(len=*), intent(in) :: dir_landdata
   integer         , intent(in) :: lc_year

   ! Local variables
   character(len=256) :: filename, fileblock, cyear
   integer :: ie, je, nelm, totlen, tothis, iblk, jblk, owner_rank, i
   integer,   allocatable :: nelm_rank(:), ndsp_rank(:)
   integer*8, allocatable :: elmindx(:)
   integer,   allocatable :: npxlall(:)
   integer,   allocatable :: elmpixels(:,:)
   real(r8),  allocatable :: lon(:), lat(:)

   integer :: nsend, nrecv, ndone, ndsp

#ifdef MPAS_EMBEDDED_COLM
      CALL CoLM_stop('MPAS embedded CoLM does not support standalone mesh_save_to_file.')
#endif

      ! add parameter input for time year
      write(cyear,'(i4.4)') lc_year
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_root) THEN
         write(*,*) 'Saving land elements ...'
         CALL system('mkdir -p ' // trim(dir_landdata) // '/mesh/' // trim(cyear))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      filename = trim(dir_landdata) // '/mesh/' //trim(cyear) // '/mesh.nc'

      DO jblk = 1, gblock%nyblk
         DO iblk = 1, gblock%nxblk

#ifdef USEMPI
            IF (p_is_compute) THEN
               IF (gblock%pio(iblk,jblk) == p_address_active(p_my_group)) THEN
#endif
                  nelm = 0
                  totlen = 0
                  DO ie = 1, numelm
                     IF ((mesh(ie)%xblk == iblk) .and. (mesh(ie)%yblk == jblk)) THEN
                        nelm = nelm + 1
                        totlen = totlen + mesh(ie)%npxl
                     ENDIF
                  ENDDO

                  IF (nelm > 0) THEN

                     allocate (elmindx (nelm))
                     allocate (npxlall (nelm))
                     allocate (elmpixels (2,totlen))

                     je = 0
                     ndsp = 0
                     DO ie = 1, numelm
                        IF ((mesh(ie)%xblk == iblk) .and. (mesh(ie)%yblk == jblk)) THEN
                           je = je + 1
                           elmindx(je) = mesh(ie)%indx
                           npxlall(je) = mesh(ie)%npxl

                           elmpixels(1,ndsp+1:ndsp+npxlall(je)) = mesh(ie)%ilon
                           elmpixels(2,ndsp+1:ndsp+npxlall(je)) = mesh(ie)%ilat

                           ndsp = ndsp + npxlall(je)
                        ENDIF
                     ENDDO
                  ENDIF

#ifdef USEMPI
                  CALL mpi_gather (nelm, 1, MPI_INTEGER, &
                     MPI_INULL_P, 1, MPI_INTEGER, p_root, p_comm_group, p_err)

                  CALL mpi_gatherv (elmindx, nelm, MPI_INTEGER8, &
                     MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER8, & ! unused on non-root ranks
                     p_root, p_comm_group, p_err)

                  CALL mpi_gatherv (npxlall, nelm, MPI_INTEGER, &
                     MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER, & ! unused on non-root ranks
                     p_root, p_comm_group, p_err)

                  ndone = 0
                  DO WHILE (ndone < totlen)
                     nsend = max(min(totlen-ndone, MesgMaxSize/8), 1)
                     CALL mpi_send (nsend, 1, &
                        MPI_INTEGER, p_root, mpi_tag_size, p_comm_group, p_err)
                     CALL mpi_send (elmpixels(:,ndone+1:ndone+nsend), 2*nsend, &
                        MPI_INTEGER, p_root, mpi_tag_data, p_comm_group, p_err)
                     ndone = ndone + nsend
                  ENDDO
               ENDIF
            ENDIF
#endif

#ifdef USEMPI
            IF (p_is_active) THEN
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

                  allocate (nelm_rank (0:p_np_group-1))
                  nelm_rank(0) = 0
                  CALL mpi_gather (MPI_IN_PLACE, 0, MPI_INTEGER, &
                     nelm_rank, 1, MPI_INTEGER, p_root, p_comm_group, p_err)

                  nelm = sum(nelm_rank)

                  allocate (ndsp_rank(0:p_np_group-1))
                  ndsp_rank(0) = 0
                  DO owner_rank = 1, p_np_group-1
                     ndsp_rank(owner_rank) = ndsp_rank(owner_rank-1) + nelm_rank(owner_rank-1)
                  ENDDO

                  allocate (elmindx (nelm))
                  CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER8, &
                     elmindx, nelm_rank(0:), ndsp_rank(0:), MPI_INTEGER8, &
                     p_root, p_comm_group, p_err)

                  allocate (npxlall (nelm))
                  CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
                     npxlall, nelm_rank(0:), ndsp_rank(0:), MPI_INTEGER, &
                     p_root, p_comm_group, p_err)

                  totlen = sum(npxlall)
                  allocate (elmpixels (2, totlen))

                  ndone = 0
                  DO owner_rank = 1, p_np_group-1

                     ndsp   = ndsp_rank(owner_rank)
                     tothis = ndone + sum(npxlall(ndsp+1:ndsp+nelm_rank(owner_rank)))

                     DO WHILE (ndone < tothis)

                        CALL mpi_recv (nrecv, 1, &
                           MPI_INTEGER, owner_rank, mpi_tag_size, p_comm_group, p_stat, p_err)
                        CALL mpi_recv (elmpixels(:,ndone+1:ndone+nrecv), 2*nrecv, &
                           MPI_INTEGER, owner_rank, mpi_tag_data, p_comm_group, p_stat, p_err)

                        ndone = ndone + nrecv
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
#endif

            IF (p_is_active) THEN
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (nelm > 0) THEN
                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     CALL ncio_create_file (fileblock)

                     CALL ncio_define_dimension (fileblock, 'element',nelm)
                     CALL ncio_define_dimension (fileblock, 'ncoor',  2   )
                     CALL ncio_define_dimension (fileblock, 'pixel',  totlen)

                     CALL ncio_write_serial (fileblock, 'elmindex',  elmindx, 'element')
                     CALL ncio_write_serial (fileblock, 'elmnpxl',   npxlall, 'element')
                     CALL ncio_write_serial (fileblock, 'elmpixels', elmpixels, &
                        'ncoor', 'pixel', compress = 1)
                  ENDIF
               ENDIF
            ENDIF

            IF (allocated (elmindx))   deallocate(elmindx)
            IF (allocated (npxlall))   deallocate(npxlall)
            IF (allocated (elmpixels)) deallocate(elmpixels)

            IF (allocated (nelm_rank)) deallocate(nelm_rank)
            IF (allocated (ndsp_rank)) deallocate(ndsp_rank)

#ifdef USEMPI
            CALL mpi_barrier (p_comm_group, p_err)
#endif
         ENDDO
      ENDDO

      IF (p_is_root) THEN

         CALL ncio_create_file (filename)

         CALL ncio_define_dimension (filename, 'xblk', gblock%nxblk)
         CALL ncio_define_dimension (filename, 'yblk', gblock%nyblk)
         CALL ncio_write_serial (filename, 'nelm_blk', nelm_blk, 'xblk', 'yblk')

         CALL ncio_define_dimension (filename, 'longitude', gridmesh%nlon)
         CALL ncio_define_dimension (filename, 'latitude' , gridmesh%nlat)

         allocate (lon (gridmesh%nlon))
         allocate (lat (gridmesh%nlat))

         DO i = 1, gridmesh%nlon
            lon(i) = (gridmesh%lon_w(i) + gridmesh%lon_e(i)) * 0.5
            IF (gridmesh%lon_w(i) > gridmesh%lon_e(i)) THEN
               lon(i) = lon(i) + 180.0
               CALL normalize_longitude (lon(i))
            ENDIF
         ENDDO
         CALL ncio_write_serial (filename, 'longitude', lon, 'longitude')

         DO i = 1, gridmesh%nlat
            lat(i) = (gridmesh%lat_s(i) + gridmesh%lat_n(i)) * 0.5
         ENDDO
         CALL ncio_write_serial (filename, 'latitude', lat, 'latitude')

#ifdef GRIDBASED
         CALL ncio_write_serial (filename, 'lat_s', gridmesh%lat_s, 'latitude' )
         CALL ncio_write_serial (filename, 'lat_n', gridmesh%lat_n, 'latitude' )
         CALL ncio_write_serial (filename, 'lon_w', gridmesh%lon_w, 'longitude')
         CALL ncio_write_serial (filename, 'lon_e', gridmesh%lon_e, 'longitude')
#endif

         deallocate (lon)
         deallocate (lat)

      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_root) write(*,*) 'SAVE land elements done.'

   END SUBROUTINE mesh_save_to_file

   !------------------------------------
   SUBROUTINE prepare_subset_eindex(subset_eindex, sorted_subset)

   USE MOD_Utils, only: quicksort
   IMPLICIT NONE

   integer*8, intent(in) :: subset_eindex(:)
   integer*8, allocatable, intent(out) :: sorted_subset(:)

   integer, allocatable :: order(:)
   integer :: i

      allocate(sorted_subset(size(subset_eindex)))
      sorted_subset = subset_eindex

      IF (size(sorted_subset) > 1) THEN
         allocate(order(size(sorted_subset)))
         order = (/ (i, i = 1, size(sorted_subset)) /)
         CALL quicksort(size(sorted_subset), sorted_subset, order)
         deallocate(order)
      ENDIF

   END SUBROUTINE prepare_subset_eindex

   !------------------------------------
   logical FUNCTION eindex_in_subset(eindex, sorted_subset) result(keep)

   USE MOD_Utils, only: find_in_sorted_list1
   IMPLICIT NONE

   integer*8, intent(in) :: eindex
   integer*8, allocatable, intent(in) :: sorted_subset(:)

      IF (.not. allocated(sorted_subset)) THEN
         keep = .true.
      ELSEIF (size(sorted_subset) < 1) THEN
         keep = .false.
      ELSE
         keep = find_in_sorted_list1(eindex, size(sorted_subset), sorted_subset) > 0
      ENDIF

   END FUNCTION eindex_in_subset

   !------------------------------------
   integer FUNCTION count_subset_eindex(eindex, sorted_subset) result(nkeep)

   IMPLICIT NONE

   integer*8, intent(in) :: eindex(:)
   integer*8, allocatable, intent(in) :: sorted_subset(:)

   integer :: i

      nkeep = 0
      DO i = 1, size(eindex)
         IF (eindex_in_subset(eindex(i), sorted_subset)) nkeep = nkeep + 1
      ENDDO

   END FUNCTION count_subset_eindex

   !------------------------------------
   SUBROUTINE mesh_load_from_file (dir_landdata, lc_year, subset_eindex)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_NetCDFSerial
   USE MOD_Mesh
   IMPLICIT NONE

   integer         , intent(in) :: lc_year
   character(len=*), intent(in) :: dir_landdata
   integer*8, optional, intent(in) :: subset_eindex(:)

   ! Local variables
   character(len=256) :: filename, fileblock, cyear
   integer :: iblkme, iblk, jblk, ie, nelm, ndsp, pdsp
   integer*8, allocatable :: elmindx(:)
   integer*8, allocatable :: subset_sorted(:)
   integer,   allocatable :: datasize(:)
   integer,   allocatable :: npxl(:), pixels(:,:), pixels2d(:,:,:)
   logical,   allocatable :: keep_elm(:)
   logical :: use_subset

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_root) THEN
         write(*,*) 'Loading land elements ...'
      ENDIF

      ! add parameter input for time year
      write(cyear,'(i4.4)') lc_year
      filename = trim(dir_landdata) // '/mesh/' // trim(cyear) // '/mesh.nc'
      CALL ncio_read_bcast_serial (filename, 'nelm_blk', nelm_blk)
      use_subset = present(subset_eindex)
      IF (use_subset) CALL prepare_subset_eindex(subset_eindex, subset_sorted)

      IF (p_is_active) THEN

         numelm = sum(nelm_blk, mask = gblock%pio == p_iam_glb)

         IF (numelm > 0) THEN

            IF (allocated(mesh)) deallocate(mesh)
            allocate (mesh (numelm))

            ndsp = 0
            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               nelm = nelm_blk(iblk,jblk)

               IF (nelm > 0) THEN

                  CALL get_filename_block (filename, iblk, jblk, fileblock)
                  CALL ncio_read_serial (fileblock, 'elmindex',  elmindx)
                  CALL ncio_read_serial (fileblock, 'elmnpxl',   npxl   )

                  IF (use_subset) THEN
                     allocate (keep_elm(nelm))
                     DO ie = 1, nelm
                        keep_elm(ie) = eindex_in_subset(elmindx(ie), subset_sorted)
                     ENDDO
                     IF (count(keep_elm) < 1) THEN
                        deallocate(keep_elm)
                        IF (allocated(elmindx)) deallocate(elmindx)
                        IF (allocated(npxl)) deallocate(npxl)
                        CYCLE
                     ENDIF
                  ENDIF

                  CALL ncio_inquire_varsize (fileblock, 'elmpixels', datasize)
                  IF (size(datasize) == 3) THEN
                     CALL ncio_read_serial (fileblock, 'elmpixels', pixels2d)
                  ELSE
                     CALL ncio_read_serial (fileblock, 'elmpixels', pixels)
                  ENDIF

                  pdsp = 0
                  DO ie = 1, nelm
                     IF (use_subset) THEN
                        IF (.not. keep_elm(ie)) THEN
                           IF (size(datasize) /= 3) pdsp = pdsp + npxl(ie)
                           CYCLE
                        ENDIF
                     ENDIF

                     ndsp = ndsp + 1
                     mesh(ndsp)%indx = elmindx(ie)
                     mesh(ndsp)%npxl = npxl(ie)
                     mesh(ndsp)%xblk = iblk
                     mesh(ndsp)%yblk = jblk

                     allocate (mesh(ndsp)%ilon (npxl(ie)))
                     allocate (mesh(ndsp)%ilat (npxl(ie)))

                     IF (size(datasize) == 3) THEN
                        mesh(ndsp)%ilon = pixels2d(1,1:npxl(ie),ie)
                        mesh(ndsp)%ilat = pixels2d(2,1:npxl(ie),ie)
                     ELSE
                        mesh(ndsp)%ilon = pixels(1,pdsp+1:pdsp+npxl(ie))
                        mesh(ndsp)%ilat = pixels(2,pdsp+1:pdsp+npxl(ie))
                        pdsp = pdsp + npxl(ie)
                     ENDIF
                  ENDDO

                  IF (allocated(keep_elm)) deallocate(keep_elm)
                  IF (allocated(elmindx)) deallocate(elmindx)
                  IF (allocated(npxl)) deallocate(npxl)
                  IF (allocated(datasize)) deallocate(datasize)
                  IF (allocated(pixels)) deallocate(pixels)
                  IF (allocated(pixels2d)) deallocate(pixels2d)
               ENDIF
            ENDDO

            numelm = ndsp
         ENDIF

         IF (use_subset .and. allocated(nelm_blk)) THEN
            nelm_blk(:,:) = 0
            DO ie = 1, numelm
               nelm_blk(mesh(ie)%xblk, mesh(ie)%yblk) = &
                  nelm_blk(mesh(ie)%xblk, mesh(ie)%yblk) + 1
            ENDDO
         ENDIF

         IF (allocated(elmindx ))  deallocate(elmindx )
         IF (allocated(npxl    ))  deallocate(npxl    )
         IF (allocated(datasize))  deallocate(datasize)
         IF (allocated(pixels  ))  deallocate(pixels  )
         IF (allocated(pixels2d))  deallocate(pixels2d)
         IF (allocated(subset_sorted)) deallocate(subset_sorted)

      ENDIF

#ifdef CoLMDEBUG
      IF (p_is_active) write(*,'(I10,A,I4)') numelm, ' elements on group ', p_iam_active
#endif

#ifdef USEMPI
#ifndef MPAS_EMBEDDED_COLM
      CALL scatter_mesh_legacy_roles
#endif
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_root) THEN
         write(*,*) 'Loading land elements done.'
      ENDIF

   END SUBROUTINE mesh_load_from_file

   !------------------------------------------------
   SUBROUTINE pixelset_save_to_file (dir_landdata, psetname, pixelset, lc_year)

   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_NetCDFVector
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*),    intent(in) :: dir_landdata
   character(len=*),    intent(in) :: psetname
   type(pixelset_type), intent(in) :: pixelset
   integer         ,    intent(in) :: lc_year

   ! Local variables
   character(len=256)   :: filename, cyear

#ifdef MPAS_EMBEDDED_COLM
      CALL CoLM_stop('MPAS embedded CoLM does not support standalone pixelset_save_to_file.')
#endif

      write(cyear,'(i4.4)') lc_year
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_root) THEN
         write(*,*) 'Saving Pixel Sets ' // trim(psetname) // ' ...'
         CALL system('mkdir -p ' // trim(dir_landdata) // '/' // trim(psetname) // '/' // trim(cyear))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      filename = trim(dir_landdata) // '/' // trim(psetname) // '/' // trim(cyear) // '/' // trim(psetname) // '.nc'

      CALL ncio_create_file_vector (filename, pixelset)
      CALL ncio_define_dimension_vector (filename, pixelset, trim(psetname))

      CALL ncio_write_vector (filename, 'eindex', trim(psetname), pixelset, pixelset%eindex, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (filename, 'ipxstt', trim(psetname), pixelset, pixelset%ipxstt, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (filename, 'ipxend', trim(psetname), pixelset, pixelset%ipxend, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (filename, 'settyp', trim(psetname), pixelset, pixelset%settyp, DEF_Srfdata_CompressLevel)

      IF (pixelset%has_shared) THEN
         CALL ncio_write_vector (filename, 'pctshared', trim(psetname), pixelset, pixelset%pctshared, DEF_Srfdata_CompressLevel)
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_root) write(*,*) 'SAVE Pixel Sets ' // trim(psetname) // ' done.'

   END SUBROUTINE pixelset_save_to_file


   !---------------------------
   SUBROUTINE pixelset_load_from_file (dir_landdata, psetname, pixelset, numset, lc_year, subset_eindex)

   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_NetCDFSerial
   USE MOD_NetCDFVector
   USE MOD_Mesh
   USE MOD_Pixelset
   IMPLICIT NONE

   integer         ,    intent(in) :: lc_year
   character(len=*),    intent(in) :: dir_landdata
   character(len=*),    intent(in) :: psetname
   type(pixelset_type), intent(inout) :: pixelset
   integer, intent(out) :: numset
   integer*8, optional, intent(in) :: subset_eindex(:)

   ! Local variables
   character(len=256) :: filename, fileblock, blockname, cyear
	   integer :: iset, nset, nset_file, ndsp, iblkme, iblk, jblk, ie, je, nave, nres, left, iproc, ipos
   integer :: nsend, nrecv
   integer*8, allocatable :: rbuff(:), sbuff(:)
   integer*8, allocatable :: subset_sorted(:)
   integer,   allocatable :: owner_rank(:)
   logical,   allocatable :: msk(:)
   logical,   allocatable :: keep_set(:)
   logical :: fexists, fexists_any
   logical :: use_subset

      write(cyear,'(i4.4)') lc_year
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      use_subset = present(subset_eindex)
      IF (use_subset) CALL prepare_subset_eindex(subset_eindex, subset_sorted)

      IF (p_is_root) THEN
         write(*,*) 'Loading Pixel Sets ' // trim(psetname) // ' ...'
      ENDIF

      filename = trim(dir_landdata) // '/' // trim(psetname) // '/' // trim(cyear) // '/' // trim(psetname) // '.nc'

      IF (p_is_active) THEN

         pixelset%nset = 0

         fexists_any = .false.

         DO iblkme = 1, gblock%nblkme
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)

#if (defined VectorInOneFileS || defined VectorInOneFileP)
            CALL get_blockname (iblk, jblk, blockname)
            CALL ncio_inquire_length_grp (filename, 'eindex', &
               trim(psetname)//'_'//trim(blockname), nset)
            IF (use_subset .and. nset > 0) THEN
               CALL ncio_read_serial_grp_int64_1d (filename, 'eindex', &
                  trim(psetname)//'_'//trim(blockname), rbuff)
               nset = count_subset_eindex(rbuff, subset_sorted)
               deallocate(rbuff)
            ENDIF
            pixelset%nset = pixelset%nset + nset
#else
            CALL get_filename_block (filename, iblk, jblk, fileblock)

            inquire (file=trim(fileblock), exist=fexists)
            IF (fexists) THEN
               IF (use_subset) THEN
                  CALL ncio_read_serial (fileblock, 'eindex', rbuff)
                  nset = count_subset_eindex(rbuff, subset_sorted)
                  deallocate(rbuff)
               ELSE
                  CALL ncio_inquire_length (fileblock, 'eindex', nset)
               ENDIF
               pixelset%nset = pixelset%nset + nset
            ENDIF

            fexists_any = fexists_any .or. fexists
#endif
         ENDDO

#if (defined VectorInOneFileS || defined VectorInOneFileP)
         inquire(file=trim(filename), exist=fexists_any)
#endif

#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, fexists_any, 1, MPI_LOGICAL, MPI_LOR, p_comm_active, p_err)
#endif
         IF (.not. fexists_any) THEN
            write(*,*) 'Warning : restart file ' //trim(filename)// ' not found.'
            CALL CoLM_stop ()
         ENDIF

         IF (pixelset%nset > 0) THEN

	            allocate (pixelset%eindex (pixelset%nset))
	            allocate (pixelset%srcpos (pixelset%nset))

            ndsp = 0
            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

#if (defined VectorInOneFileS || defined VectorInOneFileP)
               CALL get_blockname (iblk, jblk, blockname)
               CALL ncio_inquire_length_grp (filename, 'eindex', &
                  trim(psetname)//'_'//trim(blockname), nset)

               IF (nset > 0) THEN

                  CALL ncio_read_serial_grp_int64_1d (filename, 'eindex', &
                     trim(psetname)//'_'//trim(blockname), rbuff)

                  nset_file = size(rbuff)
                  IF (use_subset) THEN
                     allocate(keep_set(nset_file))
                     DO iset = 1, nset_file
                        keep_set(iset) = eindex_in_subset(rbuff(iset), subset_sorted)
                     ENDDO
                     nset = count(keep_set)
                  ELSE
                     nset = nset_file
                  ENDIF

                  IF (nset > 0) THEN
	                     IF (use_subset) THEN
	                        pixelset%eindex(ndsp+1:ndsp+nset) = pack(rbuff, keep_set)
	                        pixelset%srcpos(ndsp+1:ndsp+nset) = &
	                           pack((/ (ipos, ipos = 1, nset_file) /), keep_set)
	                     ELSE
	                        pixelset%eindex(ndsp+1:ndsp+nset) = rbuff
	                        pixelset%srcpos(ndsp+1:ndsp+nset) = (/ (ipos, ipos = 1, nset) /)
	                     ENDIF

	                     ndsp = ndsp + nset
                  ENDIF

                  IF (allocated(keep_set)) deallocate(keep_set)
                  deallocate(rbuff)
               ENDIF
#else
               CALL get_filename_block (filename, iblk, jblk, fileblock)
               inquire (file=trim(fileblock), exist=fexists)
               IF (fexists) THEN

                  CALL ncio_read_serial (fileblock, 'eindex', rbuff)

	                  nset_file = size(rbuff)
	                  IF (use_subset) THEN
	                     allocate(keep_set(nset_file))
	                     DO iset = 1, nset_file
	                        keep_set(iset) = eindex_in_subset(rbuff(iset), subset_sorted)
	                     ENDDO
	                     nset = count(keep_set)
	                  ELSE
	                     nset = nset_file
	                  ENDIF

	                  IF (nset > 0) THEN
	                     IF (use_subset) THEN
	                        pixelset%eindex(ndsp+1:ndsp+nset) = pack(rbuff, keep_set)
	                        pixelset%srcpos(ndsp+1:ndsp+nset) = &
	                           pack((/ (ipos, ipos = 1, nset_file) /), keep_set)
	                     ELSE
	                        pixelset%eindex(ndsp+1:ndsp+nset) = rbuff
	                        pixelset%srcpos(ndsp+1:ndsp+nset) = (/ (ipos, ipos = 1, nset) /)
	                     ENDIF

	                     ndsp = ndsp + nset
	                  ENDIF

	                  IF (allocated(keep_set)) deallocate(keep_set)
	                  deallocate(rbuff)
               ENDIF
#endif

            ENDDO
         ENDIF
      ENDIF


#if defined(USEMPI) && !defined(MPAS_EMBEDDED_COLM)
      IF (p_is_active) THEN
         IF (pixelset%nset > 0) THEN
            allocate (owner_rank (pixelset%nset))
            allocate (msk     (pixelset%nset))

            ie = 1
            je = 1
            iblk = mesh(ie)%xblk
            jblk = mesh(ie)%yblk
            DO iset = 1, pixelset%nset
               DO WHILE (pixelset%eindex(iset) /= mesh(ie)%indx)
                  ie = ie + 1
                  je = je + 1
                  IF ((mesh(ie)%xblk /= iblk) .or. (mesh(ie)%yblk /= jblk)) THEN
                     je = 1
                     iblk = mesh(ie)%xblk
                     jblk = mesh(ie)%yblk
                  ENDIF
               ENDDO

               nave = nelm_blk(iblk,jblk) / (p_np_group-1)
               nres = mod(nelm_blk(iblk,jblk), p_np_group-1)
               left = (nave+1) * nres
               IF (je <= left) THEN
                  owner_rank(iset) = (je-1) / (nave+1) + 1
               ELSE
                  owner_rank(iset) = (je-left-1) / nave + 1 + nres
               ENDIF
            ENDDO

            DO iproc = 1, p_np_group-1
               msk = (owner_rank == iproc)
               nsend = count(msk)
               CALL mpi_send (nsend, 1, MPI_INTEGER, iproc, mpi_tag_size, p_comm_group, p_err)

               IF (nsend > 0) THEN
                  allocate (sbuff(nsend))
                  sbuff = pack(pixelset%eindex, msk)
                  CALL mpi_send (sbuff, nsend, MPI_INTEGER8, iproc, mpi_tag_data, p_comm_group, p_err)
                  deallocate (sbuff)
               ENDIF
            ENDDO
         ELSE
            DO iproc = 1, p_np_group-1
               nsend = 0
               CALL mpi_send (nsend, 1, MPI_INTEGER, iproc, mpi_tag_size, p_comm_group, p_err)
            ENDDO
         ENDIF

      ENDIF

      IF (p_is_compute) THEN

         CALL mpi_recv (nrecv, 1, MPI_INTEGER, p_root, mpi_tag_size, p_comm_group, p_stat, p_err)

         pixelset%nset = nrecv
         IF (nrecv > 0) THEN
            allocate (pixelset%eindex (nrecv))
            CALL mpi_recv (pixelset%eindex, nrecv, MPI_INTEGER8, &
               p_root, mpi_tag_data, p_comm_group, p_stat, p_err)
         ENDIF
      ENDIF
#endif


      CALL pixelset%set_vecgs

      CALL ncio_read_vector (filename, 'ipxstt', pixelset, pixelset%ipxstt)
      CALL ncio_read_vector (filename, 'ipxend', pixelset, pixelset%ipxend)
      CALL ncio_read_vector (filename, 'settyp', pixelset, pixelset%settyp)

      IF (p_is_compute) THEN
         IF (pixelset%nset > 0) THEN

            allocate (pixelset%ielm (pixelset%nset))
            ie = 1
            DO iset = 1, pixelset%nset
               DO WHILE (pixelset%eindex(iset) /= mesh(ie)%indx)
                  ie = ie + 1
               ENDDO
               pixelset%ielm(iset) = ie
            ENDDO

         ELSE
            write(*,*) 'Warning: 0 ',trim(psetname), ' on rank :', p_iam_glb
         ENDIF
      ENDIF

      numset = pixelset%nset

      pixelset%has_shared = .false.
      IF (p_is_compute) THEN
         DO iset = 1, pixelset%nset-1
            IF ((pixelset%ielm(iset) == pixelset%ielm(iset+1)) &
               .and. (pixelset%ipxstt(iset) == pixelset%ipxstt(iset+1))) THEN
               pixelset%has_shared = .true.
               exit
            ENDIF
         ENDDO
      ENDIF

#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, pixelset%has_shared, 1, MPI_LOGICAL, &
         MPI_LOR, p_comm_glb, p_err)
#endif

      IF (pixelset%has_shared) THEN
         CALL ncio_read_vector (filename, 'pctshared', pixelset, pixelset%pctshared)
      ENDIF

#ifdef CoLMDEBUG
      IF (p_is_active)  write(*,*) numset, trim(psetname), ' on group', p_iam_active
#endif

      IF (allocated(subset_sorted)) deallocate(subset_sorted)

   END SUBROUTINE pixelset_load_from_file

END MODULE MOD_SrfdataRestart
