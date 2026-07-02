#include <define.h>

MODULE MOD_CoLM_MPAS_Interface

   USE MOD_Precision
   USE MOD_LandPatch, only: numpatch, landpatch, elm_patch
   USE MOD_Vars_Global, only: spval
   USE MOD_Vars_1DForcing, only: forc_pco2m, forc_po2m, forc_us, forc_vs, forc_t, forc_q, &
      forc_prc, forc_prl, forc_rain, forc_snow, forc_psrf, forc_pbot, forc_sols, forc_soll, &
      forc_solsd, forc_solld, forc_frl, forc_swrad, forc_hgt_u, forc_hgt_t, forc_hgt_q, &
      forc_rhoair, forc_ozone, forc_hpbl, forc_aerdep
#ifdef HYPERSPECTRAL
   USE MOD_Vars_1DForcing, only: forc_solarin
#endif
   USE MOD_Vars_1DFluxes, only: oroflag, fsena, lfevpa, fevpa, fgrnd, rnof, rsur, rsub
   USE MOD_Vars_TimeInvariants, only: patchlatr, patchlonr, patchmask
   USE MOD_Vars_TimeVariables, only: t_grnd, tref, qref, emis, z0m, alb
   USE MOD_TimeManager, only: timestamp

   IMPLICIT NONE
   PRIVATE

	   PUBLIC :: colm_mpas_initialize_from_namelist
	   PUBLIC :: colm_mpas_finalize
	   PUBLIC :: colm_mpas_ready
	   PUBLIC :: colm_mpas_find_patch
	   PUBLIC :: colm_mpas_find_element
   PUBLIC :: colm_mpas_set_forcing
   PUBLIC :: colm_mpas_set_element_forcing
   PUBLIC :: colm_mpas_step
   PUBLIC :: colm_mpas_get_surface
   PUBLIC :: colm_mpas_get_element_surface

   logical, save :: colm_mpas_initialized = .false.
   real(r8), allocatable, save :: colm_mpas_elm_lonr(:)
   real(r8), allocatable, save :: colm_mpas_elm_latr(:)
   character(len=256), save :: colm_mpas_casename = ''
   character(len=256), save :: colm_mpas_dir_restart = ''
   integer, save :: colm_mpas_lc_year = -1
   integer, save :: colm_mpas_last_idate(3) = -1
   integer, save :: colm_mpas_last_restart_idate(3) = -1
   type(timestamp), save :: colm_mpas_ptstamp
   type(timestamp), save :: colm_mpas_etstamp
   logical, save :: colm_mpas_restart_ready = .false.

   INTERFACE
      SUBROUTINE CoLMDRIVER(idate,deltim,dolai,doalb,dosst,oro)
         USE MOD_Precision
         USE MOD_LandPatch, only: numpatch
         integer, intent(in) :: idate(3)
         real(r8), intent(in) :: deltim
         logical, intent(in) :: dolai, doalb, dosst
         real(r8), intent(inout) :: oro(numpatch)
      END SUBROUTINE CoLMDRIVER
   END INTERFACE

CONTAINS

   SUBROUTINE colm_mpas_initialize_from_namelist(nlfile, ierr, mpas_comm, mpas_lon_rad, mpas_lat_rad, &
                                                mpas_cell_id, n_mpas_cells, cell_to_element)
	      USE MOD_Namelist, only: read_namelist, DEF_CASE_NAME, DEF_dir_landdata, &
		         DEF_dir_restart, DEF_LC_YEAR, DEF_simulation_time, DEF_USE_SNICAR, &
		         DEF_file_snowoptics, DEF_file_snowaging, DEF_forcing, DEF_Reservoir_Method, &
		         DEF_WRST_FREQ, DEF_HIST_FREQ, DEF_HIST_WriteBack
	      USE MOD_Vars_Global, only: Init_GlobalVars
	      USE MOD_SPMD_Task, only: spmd_init
	      USE MOD_Const_LC, only: Init_LC_Const
	      USE MOD_Const_PFT, only: Init_PFT_Const, rho_p, tau_p
      USE MOD_TimeManager, only: initimetype, monthday2julian, adj2begin, adj2end
      USE MOD_Block, only: gblock
      USE MOD_Pixel, only: pixel
      USE MOD_Mesh, only: numelm
      USE MOD_LandElm, only: landelm
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      USE MOD_LandPFT, only: landpft, numpft, map_patch_to_pft
#endif
      USE MOD_SrfdataRestart, only: mesh_load_from_file, pixelset_load_from_file
      USE MOD_Vars_TimeInvariants, only: allocate_TimeInvariants, READ_TimeInvariants
	      USE MOD_Vars_TimeVariables, only: allocate_TimeVariables, READ_TimeVariables
	      USE MOD_Vars_1DForcing, only: allocate_1D_Forcing
	      USE MOD_Vars_1DFluxes, only: allocate_1D_Fluxes
#ifdef GridRiverLakeFlow
	      USE MOD_Grid_RiverLakeNetwork, only: build_riverlake_network
	      USE MOD_Grid_Reservoir, only: reservoir_init
	      USE MOD_Grid_RiverLakeFlow, only: grid_riverlake_flow_init
#endif
#ifdef HYPERSPECTRAL
	      USE MOD_SnowSnicar_HiRes, only: SnowAge_init, SnowOptics_init
	      USE MOD_HighRes_Parameters, only: flux_frac_init, leaf_property_init, &
         get_water_optical_properties
#else
      USE MOD_SnowSnicar, only: SnowAge_init, SnowOptics_init
#endif
      character(len=*), intent(in) :: nlfile
      integer, intent(out) :: ierr
      integer, intent(in), optional :: mpas_comm
      real(r8), intent(in), optional :: mpas_lon_rad(:)
      real(r8), intent(in), optional :: mpas_lat_rad(:)
      integer, intent(in), optional :: mpas_cell_id(:)
      integer, intent(in), optional :: n_mpas_cells
      integer, intent(out), optional :: cell_to_element(:)

      character(len=256) :: casename
      character(len=256) :: dir_landdata
      character(len=256) :: dir_restart
      integer :: lc_year
      integer :: sdate(3)
      integer :: jdate(3)
      integer :: s_julian
      integer :: p_julian
      integer :: e_julian
      integer :: n_mpas
      integer :: i
      integer*8, allocatable :: mpas_cell_id_i8(:)

      ierr = 1
      IF (colm_mpas_initialized) THEN
         ierr = 0
         RETURN
      ENDIF

      IF (present(mpas_comm)) THEN
         CALL spmd_init(mpas_comm)
      ELSE
         CALL spmd_init()
      ENDIF

	      CALL read_namelist(trim(nlfile))
#ifdef MPAS_EMBEDDED_COLM
	      CALL colm_mpas_check_embedded_io(ierr)
	      IF (ierr /= 0) RETURN
#endif

	      casename = DEF_CASE_NAME
      dir_landdata = DEF_dir_landdata
      dir_restart = DEF_dir_restart
      lc_year = DEF_LC_YEAR

      CALL initimetype(DEF_simulation_time%greenwich)
      CALL monthday2julian(DEF_simulation_time%start_year, DEF_simulation_time%start_month, &
                           DEF_simulation_time%start_day, s_julian)
      CALL monthday2julian(DEF_simulation_time%spinup_year, DEF_simulation_time%spinup_month, &
                           DEF_simulation_time%spinup_day, p_julian)
      CALL monthday2julian(DEF_simulation_time%end_year, DEF_simulation_time%end_month, &
                           DEF_simulation_time%end_day, e_julian)
      sdate(1) = DEF_simulation_time%start_year
      sdate(2) = s_julian
      sdate(3) = DEF_simulation_time%start_sec

      colm_mpas_casename = casename
      colm_mpas_dir_restart = dir_restart
      colm_mpas_lc_year = lc_year
      colm_mpas_ptstamp%year = DEF_simulation_time%spinup_year
      colm_mpas_ptstamp%day = p_julian
      colm_mpas_ptstamp%sec = DEF_simulation_time%spinup_sec
      colm_mpas_etstamp%year = DEF_simulation_time%end_year
      colm_mpas_etstamp%day = e_julian
      colm_mpas_etstamp%sec = DEF_simulation_time%end_sec
      colm_mpas_last_idate(:) = -1
      colm_mpas_last_restart_idate(:) = -1
      colm_mpas_restart_ready = .true.

      CALL Init_GlobalVars
      CALL Init_LC_Const
      CALL Init_PFT_Const

      n_mpas = 0
      IF (present(n_mpas_cells)) THEN
         IF (.not. present(mpas_lon_rad)) RETURN
         IF (.not. present(mpas_lat_rad)) RETURN
         IF (.not. present(mpas_cell_id)) RETURN
         IF (.not. present(cell_to_element)) RETURN

         n_mpas = n_mpas_cells
         IF (n_mpas < 1) RETURN
         IF (size(mpas_lon_rad) < n_mpas) RETURN
         IF (size(mpas_lat_rad) < n_mpas) RETURN
         IF (size(mpas_cell_id) < n_mpas) RETURN
         IF (size(cell_to_element) < n_mpas) RETURN

         allocate(mpas_cell_id_i8(n_mpas))
         DO i = 1, n_mpas
            mpas_cell_id_i8(i) = int(mpas_cell_id(i), 8)
         ENDDO
      ENDIF

      CALL pixel%load_from_file(dir_landdata)
      CALL gblock%load_from_file(dir_landdata)

      IF (n_mpas > 0) THEN
         CALL colm_mpas_claim_owned_blocks(mpas_lon_rad, mpas_lat_rad, n_mpas, ierr)
         IF (ierr /= 0) RETURN

         CALL mesh_load_from_file(dir_landdata, lc_year, subset_eindex=mpas_cell_id_i8)
         CALL pixelset_load_from_file(dir_landdata, 'landelm', landelm, numelm, lc_year, &
                                      subset_eindex=mpas_cell_id_i8)
         CALL pixelset_load_from_file(dir_landdata, 'landpatch', landpatch, numpatch, lc_year, &
                                      subset_eindex=mpas_cell_id_i8)
      ELSE
         CALL mesh_load_from_file(dir_landdata, lc_year)
         CALL pixelset_load_from_file(dir_landdata, 'landelm', landelm, numelm, lc_year)
         CALL pixelset_load_from_file(dir_landdata, 'landpatch', landpatch, numpatch, lc_year)
      ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (n_mpas > 0) THEN
         CALL pixelset_load_from_file(dir_landdata, 'landpft', landpft, numpft, lc_year, &
                                      subset_eindex=mpas_cell_id_i8)
      ELSE
         CALL pixelset_load_from_file(dir_landdata, 'landpft', landpft, numpft, lc_year)
      ENDIF
#endif

      IF (n_mpas > 0) THEN
         CALL colm_mpas_restrict_to_mpas_cells(mpas_lon_rad, mpas_lat_rad, mpas_cell_id, &
                                               n_mpas, cell_to_element, ierr)
         IF (ierr /= 0) RETURN
      ENDIF

	      CALL elm_patch%build(landelm, landpatch, use_frac = .true.)

#ifdef GridRiverLakeFlow
#ifdef MPAS_EMBEDDED_COLM
	      CALL colm_mpas_check_embedded_riverlake(ierr)
	      IF (ierr /= 0) RETURN
#endif
	      CALL build_riverlake_network()
	      IF (DEF_Reservoir_Method > 0) CALL reservoir_init()
#endif

	      IF (allocated(colm_mpas_elm_lonr)) deallocate(colm_mpas_elm_lonr)
	      IF (allocated(colm_mpas_elm_latr)) deallocate(colm_mpas_elm_latr)
	      IF (numelm > 0) THEN
         allocate(colm_mpas_elm_lonr(numelm))
         allocate(colm_mpas_elm_latr(numelm))
         CALL landelm%get_lonlat_radian(colm_mpas_elm_lonr, colm_mpas_elm_latr)
      ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      CALL map_patch_to_pft
#endif

      CALL adj2end(sdate)
      jdate = sdate
      CALL adj2begin(jdate)

      CALL allocate_TimeInvariants()
      CALL READ_TimeInvariants(lc_year, casename, dir_restart)
      CALL allocate_TimeVariables()
      CALL READ_TimeVariables(jdate, lc_year, casename, dir_restart)

      IF (DEF_USE_SNICAR) THEN
         CALL SnowOptics_init(DEF_file_snowoptics)
         CALL SnowAge_init(DEF_file_snowaging)
      ENDIF

#ifdef HYPERSPECTRAL
      CALL flux_frac_init()
      CALL leaf_property_init(rho_p, tau_p)
      CALL get_water_optical_properties()
#endif

	      CALL allocate_1D_Forcing()
	      CALL allocate_1D_Fluxes()
	      DEF_forcing%has_missing_value = .false.

#ifdef GridRiverLakeFlow
	      CALL grid_riverlake_flow_init()
#endif

	      colm_mpas_initialized = .true.
	      ierr = 0
	   END SUBROUTINE colm_mpas_initialize_from_namelist

#if defined(GridRiverLakeFlow) && defined(MPAS_EMBEDDED_COLM)
	   SUBROUTINE colm_mpas_check_embedded_riverlake(ierr)
	      USE MOD_Namelist, only: DEF_USE_SEDIMENT
	      USE MOD_SPMD_Task, only: p_np_glb, p_is_master
	      integer, intent(out) :: ierr

	      ierr = 0
	      IF (DEF_USE_SEDIMENT) THEN
	         IF (p_is_master) THEN
	            write(*,'(A)') 'CoLM2024 MPAS embedded mode does not support GridRiverLakeSediment yet.'
	            write(*,'(A)') 'Set DEF_USE_SEDIMENT = .false. until sediment routing is migrated to MPAS-owned rank decomposition.'
	         ENDIF
	         ierr = 1
	         RETURN
	      ENDIF
	      IF (p_np_glb > 1 .and. p_is_master) THEN
	         write(*,'(A)') 'CoLM2024 MPAS embedded GridRiverLakeFlow uses MPAS communicator ranks for distributed routing.'
	         write(*,'(A)') 'CoLM standalone master/io/worker MPI pools and replicated full-river-network fallback are disabled.'
	      ENDIF
		   END SUBROUTINE colm_mpas_check_embedded_riverlake
#endif

	   SUBROUTINE colm_mpas_check_embedded_io(ierr)
	      USE MOD_Namelist, only: DEF_HIST_FREQ, DEF_HIST_WriteBack, USE_SITE_HistWriteBack
	      USE MOD_SPMD_Task, only: p_is_master
	      integer, intent(out) :: ierr

	      ierr = 0
	      USE_SITE_HistWriteBack = .false.
	      IF ((trim(adjustl(DEF_HIST_FREQ)) /= 'none' .and. trim(adjustl(DEF_HIST_FREQ)) /= 'NONE') .or. &
	          DEF_HIST_WriteBack) THEN
	         IF (p_is_master) THEN
	            write(*,'(A)') 'CoLM2024 MPAS embedded mode currently writes fluxes through MPAS streams.'
	            write(*,'(A)') 'Disable CoLM DEF_HIST_FREQ/DEF_HIST_WriteBack; CoLM restart files remain supported for patch/PFT state.'
	            write(*,'(A)') 'CoLM USE_SITE_HistWriteBack is forced off in MPAS embedded mode.'
	         ENDIF
	         ierr = 1
	      ENDIF
	   END SUBROUTINE colm_mpas_check_embedded_io

		   SUBROUTINE colm_mpas_claim_owned_blocks(mpas_lon_rad, mpas_lat_rad, n_mpas_cells, ierr)
	      USE MOD_Block, only: gblock
	      USE MOD_SPMD_Task, only: p_iam_glb
	      USE MOD_Utils, only: normalize_longitude, find_nearest_south, find_nearest_west
	      real(r8), intent(in) :: mpas_lon_rad(:)
	      real(r8), intent(in) :: mpas_lat_rad(:)
	      integer, intent(in) :: n_mpas_cells
	      integer, intent(out) :: ierr

	      logical, allocatable :: keep_block(:,:)
	      integer :: i
	      integer :: iblk
	      integer :: jblk
	      integer :: iblkme
	      real(r8) :: lon_deg
	      real(r8) :: lat_deg
	      real(r8), parameter :: rad_to_deg = 57.295779513082320876798154814105_r8

	      ierr = 1
	      IF (n_mpas_cells < 1) RETURN
	      IF (.not. allocated(gblock%pio)) RETURN
	      IF (.not. allocated(gblock%lon_w)) RETURN
	      IF (.not. allocated(gblock%lat_s)) RETURN

	      allocate(keep_block(gblock%nxblk, gblock%nyblk))
	      keep_block = .false.

	      DO i = 1, n_mpas_cells
	         lon_deg = mpas_lon_rad(i) * rad_to_deg
	         lat_deg = mpas_lat_rad(i) * rad_to_deg
	         CALL normalize_longitude(lon_deg)
	         lat_deg = max(-90._r8, min(90._r8, lat_deg))

	         iblk = find_nearest_west(lon_deg, gblock%nxblk, gblock%lon_w)
	         jblk = find_nearest_south(lat_deg, gblock%nyblk, gblock%lat_s)
	         IF (iblk >= 1 .and. iblk <= gblock%nxblk .and. jblk >= 1 .and. jblk <= gblock%nyblk) THEN
	            keep_block(iblk,jblk) = .true.
	         ENDIF
	      ENDDO

	      gblock%pio(:,:) = -1
	      WHERE (keep_block)
	         gblock%pio = p_iam_glb
	      END WHERE

	      IF (allocated(gblock%xblkme)) deallocate(gblock%xblkme)
	      IF (allocated(gblock%yblkme)) deallocate(gblock%yblkme)
	      gblock%nblkme = count(keep_block)
	      IF (gblock%nblkme < 1) THEN
	         deallocate(keep_block)
	         RETURN
	      ENDIF

	      allocate(gblock%xblkme(gblock%nblkme))
	      allocate(gblock%yblkme(gblock%nblkme))
	      iblkme = 0
	      DO iblk = 1, gblock%nxblk
	         DO jblk = 1, gblock%nyblk
	            IF (keep_block(iblk,jblk)) THEN
	               iblkme = iblkme + 1
	               gblock%xblkme(iblkme) = iblk
	               gblock%yblkme(iblkme) = jblk
	            ENDIF
	         ENDDO
	      ENDDO

	      deallocate(keep_block)
	      ierr = 0
	   END SUBROUTINE colm_mpas_claim_owned_blocks

	   SUBROUTINE colm_mpas_finalize(ierr)
#ifdef GridRiverLakeFlow
	      USE MOD_Grid_RiverLakeFlow, only: grid_riverlake_flow_final
#endif
	      USE MOD_SPMD_Task, only: spmd_exit
	      integer, intent(out) :: ierr

	      ierr = 0
	      IF (.not. colm_mpas_initialized) RETURN

	      IF (colm_mpas_last_idate(1) > 0) THEN
	         CALL colm_mpas_write_restart_if_due(colm_mpas_last_idate, 0._r8, .true., ierr)
	         IF (ierr /= 0) RETURN
	      ENDIF

#ifdef GridRiverLakeFlow
	      CALL grid_riverlake_flow_final()
#endif

	      IF (allocated(colm_mpas_elm_lonr)) deallocate(colm_mpas_elm_lonr)
	      IF (allocated(colm_mpas_elm_latr)) deallocate(colm_mpas_elm_latr)
	      CALL spmd_exit()
	      colm_mpas_initialized = .false.
	      colm_mpas_restart_ready = .false.
	   END SUBROUTINE colm_mpas_finalize

	   SUBROUTINE colm_mpas_restrict_to_mpas_cells(mpas_lon_rad, mpas_lat_rad, mpas_cell_id, &
                                               n_mpas_cells, cell_to_element, ierr)
      USE MOD_LandElm, only: landelm
      USE MOD_LandPatch, only: landpatch, numpatch
      USE MOD_Mesh, only: numelm
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      USE MOD_LandPFT, only: landpft, numpft
#endif
      real(r8), intent(in) :: mpas_lon_rad(:)
      real(r8), intent(in) :: mpas_lat_rad(:)
      integer, intent(in) :: mpas_cell_id(:)
      integer, intent(in) :: n_mpas_cells
      integer, intent(out) :: cell_to_element(:)
      integer, intent(out) :: ierr

      logical, allocatable :: keep_elm(:)
      logical, allocatable :: keep_patch(:)
      logical, allocatable :: keep_pft(:)
      integer, allocatable :: old_to_new(:)
      integer, allocatable :: old_element_for_cell(:)
      real(r8), allocatable :: all_lon(:)
      real(r8), allocatable :: all_lat(:)
      integer :: i
      integer :: old_element
      integer :: packed_count
      real(r8) :: match_dist2
      real(r8), parameter :: max_coord_match_dist2 = 1.e-10_r8

      ierr = 1
      IF (n_mpas_cells < 1) RETURN
      IF (numelm < 1) RETURN
      IF (landelm%nset < 1) RETURN

      allocate(keep_elm(numelm))
      allocate(old_to_new(numelm))
      allocate(old_element_for_cell(n_mpas_cells))
      allocate(all_lon(numelm))
      allocate(all_lat(numelm))

      keep_elm = .false.
      old_to_new = 0
      CALL landelm%get_lonlat_radian(all_lon, all_lat)

      DO i = 1, n_mpas_cells
         CALL colm_mpas_find_element_by_eindex(mpas_cell_id(i), old_element)
         IF (old_element <= 0) THEN
            match_dist2 = huge(match_dist2)
            CALL colm_mpas_find_element_in_arrays(mpas_lon_rad(i), mpas_lat_rad(i), &
                                                  all_lon, all_lat, old_element, match_dist2)
            IF (match_dist2 > max_coord_match_dist2) old_element = 0
         ENDIF
         IF (old_element <= 0 .or. old_element > numelm) RETURN

         old_element_for_cell(i) = old_element
         keep_elm(old_element) = .true.
      ENDDO

      packed_count = 0
      DO i = 1, numelm
         IF (keep_elm(i)) THEN
            packed_count = packed_count + 1
            old_to_new(i) = packed_count
         ENDIF
      ENDDO
      IF (packed_count < 1) RETURN

      CALL colm_mpas_pack_mesh(keep_elm, old_to_new)
      CALL landelm%pset_pack(keep_elm, packed_count)
      CALL colm_mpas_remap_pixelset_ielm(landelm, old_to_new, ierr)
      IF (ierr /= 0) RETURN

      IF (landpatch%nset > 0) THEN
         allocate(keep_patch(landpatch%nset))
         DO i = 1, landpatch%nset
            keep_patch(i) = landpatch%ielm(i) >= 1 .and. landpatch%ielm(i) <= size(old_to_new)
            IF (keep_patch(i)) keep_patch(i) = old_to_new(landpatch%ielm(i)) > 0
         ENDDO
         CALL landpatch%pset_pack(keep_patch, numpatch)
         CALL colm_mpas_remap_pixelset_ielm(landpatch, old_to_new, ierr)
         deallocate(keep_patch)
         IF (ierr /= 0) RETURN
      ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (landpft%nset > 0) THEN
         allocate(keep_pft(landpft%nset))
         DO i = 1, landpft%nset
            keep_pft(i) = landpft%ielm(i) >= 1 .and. landpft%ielm(i) <= size(old_to_new)
            IF (keep_pft(i)) keep_pft(i) = old_to_new(landpft%ielm(i)) > 0
         ENDDO
         CALL landpft%pset_pack(keep_pft, numpft)
         CALL colm_mpas_remap_pixelset_ielm(landpft, old_to_new, ierr)
         deallocate(keep_pft)
         IF (ierr /= 0) RETURN
      ENDIF
#endif

      DO i = 1, n_mpas_cells
         cell_to_element(i) = old_to_new(old_element_for_cell(i))
         IF (cell_to_element(i) <= 0) RETURN
      ENDDO

      ierr = 0
   END SUBROUTINE colm_mpas_restrict_to_mpas_cells

   SUBROUTINE colm_mpas_pack_mesh(keep_elm, old_to_new)
      USE MOD_Mesh, only: irregular_elm_type, mesh, numelm, copy_elm
      logical, intent(in) :: keep_elm(:)
      integer, intent(in) :: old_to_new(:)

      type(irregular_elm_type), allocatable :: packed_mesh(:)
      integer :: old_elm
      integer :: new_elm
      integer :: packed_count

      packed_count = count(keep_elm)
      allocate(packed_mesh(packed_count))

      DO old_elm = 1, size(keep_elm)
         IF (.not. keep_elm(old_elm)) CYCLE
         new_elm = old_to_new(old_elm)
         CALL copy_elm(mesh(old_elm), packed_mesh(new_elm))
      ENDDO

      IF (allocated(mesh)) deallocate(mesh)
      CALL move_alloc(packed_mesh, mesh)
      numelm = packed_count
   END SUBROUTINE colm_mpas_pack_mesh

   SUBROUTINE colm_mpas_remap_pixelset_ielm(pixelset, old_to_new, ierr)
      USE MOD_Pixelset, only: pixelset_type
      type(pixelset_type), intent(inout) :: pixelset
      integer, intent(in) :: old_to_new(:)
      integer, intent(out) :: ierr

      integer :: iset
      integer :: old_elm

      ierr = 1
      DO iset = 1, pixelset%nset
         old_elm = pixelset%ielm(iset)
         IF (old_elm < 1 .or. old_elm > size(old_to_new)) RETURN
         IF (old_to_new(old_elm) <= 0) RETURN
         pixelset%ielm(iset) = old_to_new(old_elm)
      ENDDO

      ierr = 0
   END SUBROUTINE colm_mpas_remap_pixelset_ielm

   SUBROUTINE colm_mpas_find_element_by_eindex(cell_id, element)
      USE MOD_LandElm, only: landelm
      integer, intent(in) :: cell_id
      integer, intent(out) :: element

      integer :: i
      integer*8 :: cell_id_i8

      element = 0
      cell_id_i8 = cell_id
      DO i = 1, landelm%nset
         IF (landelm%eindex(i) == cell_id_i8) THEN
            element = i
            RETURN
         ENDIF
      ENDDO
   END SUBROUTINE colm_mpas_find_element_by_eindex

   SUBROUTINE colm_mpas_find_element_in_arrays(lon_rad, lat_rad, lon_array, lat_array, element, best_dist2_out)
      real(r8), intent(in) :: lon_rad
      real(r8), intent(in) :: lat_rad
      real(r8), intent(in) :: lon_array(:)
      real(r8), intent(in) :: lat_array(:)
      integer, intent(out) :: element
      real(r8), intent(out), optional :: best_dist2_out

      integer :: i
      real(r8) :: dlon
      real(r8) :: dist2
      real(r8) :: best_dist2
      real(r8), parameter :: two_pi = 6.283185307179586476925286766559_r8

      element = 0
      IF (size(lon_array) /= size(lat_array)) RETURN

      best_dist2 = huge(best_dist2)
      DO i = 1, size(lat_array)
         dlon = abs(lon_array(i) - lon_rad)
         dlon = min(dlon, two_pi - dlon)
         dist2 = dlon * dlon + (lat_array(i) - lat_rad) * (lat_array(i) - lat_rad)
         IF (dist2 < best_dist2) THEN
            best_dist2 = dist2
            element = i
         ENDIF
      ENDDO
      IF (present(best_dist2_out)) best_dist2_out = best_dist2
   END SUBROUTINE colm_mpas_find_element_in_arrays

   SUBROUTINE colm_mpas_ready(ready, patch_count)
      logical, intent(out) :: ready
      integer, intent(out), optional :: patch_count

      ready = allocated(forc_t) .and. allocated(oroflag) .and. allocated(fsena) .and. allocated(t_grnd) &
         .and. allocated(elm_patch%substt) .and. allocated(elm_patch%subfrc) &
         .and. allocated(colm_mpas_elm_lonr) .and. allocated(colm_mpas_elm_latr)
      IF (present(patch_count)) THEN
         IF (allocated(oroflag)) THEN
            patch_count = size(oroflag)
         ELSE
            patch_count = 0
         ENDIF
      ENDIF
   END SUBROUTINE colm_mpas_ready

   SUBROUTINE colm_mpas_find_patch(lon_rad, lat_rad, patch, ierr)
      real(r8), intent(in) :: lon_rad
      real(r8), intent(in) :: lat_rad
      integer, intent(out) :: patch
      integer, intent(out) :: ierr

      integer :: i
      real(r8) :: dlon
      real(r8) :: dist2
      real(r8) :: best_dist2
      real(r8), parameter :: two_pi = 6.283185307179586476925286766559_r8

      patch = 0
      ierr = 1
      IF (.not. allocated(patchlatr)) RETURN

      best_dist2 = huge(best_dist2)
      DO i = 1, size(patchlatr)
         IF (allocated(patchmask)) THEN
            IF (.not. patchmask(i)) CYCLE
         ENDIF
         dlon = abs(patchlonr(i) - lon_rad)
         dlon = min(dlon, two_pi - dlon)
         dist2 = dlon * dlon + (patchlatr(i) - lat_rad) * (patchlatr(i) - lat_rad)
         IF (dist2 < best_dist2) THEN
            best_dist2 = dist2
            patch = i
         ENDIF
      ENDDO

      IF (patch > 0) ierr = 0
   END SUBROUTINE colm_mpas_find_patch

   SUBROUTINE colm_mpas_find_element(lon_rad, lat_rad, element, ierr)
      real(r8), intent(in) :: lon_rad
      real(r8), intent(in) :: lat_rad
      integer, intent(out) :: element
      integer, intent(out) :: ierr

      integer :: i
      real(r8) :: dlon
      real(r8) :: dist2
      real(r8) :: best_dist2
      real(r8), parameter :: two_pi = 6.283185307179586476925286766559_r8

      element = 0
      ierr = 1
      IF (.not. allocated(colm_mpas_elm_latr)) RETURN

      best_dist2 = huge(best_dist2)
      DO i = 1, size(colm_mpas_elm_latr)
         dlon = abs(colm_mpas_elm_lonr(i) - lon_rad)
         dlon = min(dlon, two_pi - dlon)
         dist2 = dlon * dlon + (colm_mpas_elm_latr(i) - lat_rad) * (colm_mpas_elm_latr(i) - lat_rad)
         IF (dist2 < best_dist2) THEN
            best_dist2 = dist2
            element = i
         ENDIF
      ENDDO

      IF (element > 0) ierr = 0
   END SUBROUTINE colm_mpas_find_element

   SUBROUTINE colm_mpas_set_forcing(patch, pco2m, po2m, us, vs, tair, qair, prc, prl, rain, snow, &
                                    psrf, pbot, sols, soll, solsd, solld, frl, hgt_u, hgt_t, hgt_q, &
                                    rhoair, hpbl, aerdep, oro, ozone, ierr)
      integer, intent(in) :: patch
      real(r8), intent(in) :: pco2m, po2m, us, vs, tair, qair, prc, prl, rain, snow
      real(r8), intent(in) :: psrf, pbot, sols, soll, solsd, solld, frl
      real(r8), intent(in) :: hgt_u, hgt_t, hgt_q, rhoair, hpbl
      real(r8), intent(in) :: aerdep(14)
      real(r8), intent(in), optional :: oro
      real(r8), intent(in), optional :: ozone
      integer, intent(out) :: ierr

      ierr = 1
      IF (.not. allocated(forc_t)) RETURN
      IF (patch < 1 .or. patch > size(forc_t)) RETURN

      forc_pco2m(patch) = pco2m
      forc_po2m (patch) = po2m
      forc_us   (patch) = us
      forc_vs   (patch) = vs
      forc_t    (patch) = tair
      forc_q    (patch) = qair
      forc_prc  (patch) = prc
      forc_prl  (patch) = prl
      forc_rain (patch) = rain
      forc_snow (patch) = snow
      forc_psrf (patch) = psrf
      forc_pbot (patch) = pbot
      forc_sols (patch) = sols
      forc_soll (patch) = soll
      forc_solsd(patch) = solsd
      forc_solld(patch) = solld
      forc_frl  (patch) = frl
      forc_swrad(patch) = sols + soll + solsd + solld
#ifdef HYPERSPECTRAL
      forc_solarin(patch) = forc_swrad(patch)
#endif
      forc_hgt_u(patch) = hgt_u
      forc_hgt_t(patch) = hgt_t
      forc_hgt_q(patch) = hgt_q
      forc_rhoair(patch) = rhoair
      forc_hpbl (patch) = hpbl
      forc_aerdep(:,patch) = aerdep(:)
      IF (present(oro) .and. allocated(oroflag)) oroflag(patch) = oro
      IF (present(ozone)) THEN
         forc_ozone(patch) = ozone
      ELSE
         forc_ozone(patch) = 0._r8
      ENDIF
      ierr = 0
   END SUBROUTINE colm_mpas_set_forcing

   SUBROUTINE colm_mpas_set_element_forcing(element, pco2m, po2m, us, vs, tair, qair, prc, prl, rain, snow, &
                                            psrf, pbot, sols, soll, solsd, solld, frl, hgt_u, hgt_t, hgt_q, &
                                            rhoair, hpbl, aerdep, oro, ozone, ierr)
      integer, intent(in) :: element
      real(r8), intent(in) :: pco2m, po2m, us, vs, tair, qair, prc, prl, rain, snow
      real(r8), intent(in) :: psrf, pbot, sols, soll, solsd, solld, frl
      real(r8), intent(in) :: hgt_u, hgt_t, hgt_q, rhoair, hpbl
      real(r8), intent(in) :: aerdep(14)
      real(r8), intent(in), optional :: oro
      real(r8), intent(in), optional :: ozone
      integer, intent(out) :: ierr

      integer :: patch
      integer :: patch_ierr
      integer :: istt
      integer :: iend
      logical :: did_set

      ierr = 1
      IF (.not. allocated(elm_patch%substt)) RETURN
      IF (.not. allocated(elm_patch%subfrc)) RETURN
      IF (element < 1 .or. element > size(elm_patch%substt)) RETURN

      istt = elm_patch%substt(element)
      iend = elm_patch%subend(element)
      IF (istt < 1 .or. iend < istt) RETURN

      did_set = .false.
      DO patch = istt, iend
         IF (patch < 1 .or. patch > numpatch) CYCLE
         IF (allocated(patchmask)) THEN
            IF (.not. patchmask(patch)) CYCLE
         ENDIF
         IF (present(oro) .and. present(ozone)) THEN
            CALL colm_mpas_set_forcing(patch, pco2m, po2m, us, vs, tair, qair, prc, prl, rain, snow, &
                                       psrf, pbot, sols, soll, solsd, solld, frl, hgt_u, hgt_t, hgt_q, &
                                       rhoair, hpbl, aerdep, oro=oro, ozone=ozone, ierr=patch_ierr)
         ELSEIF (present(oro)) THEN
            CALL colm_mpas_set_forcing(patch, pco2m, po2m, us, vs, tair, qair, prc, prl, rain, snow, &
                                       psrf, pbot, sols, soll, solsd, solld, frl, hgt_u, hgt_t, hgt_q, &
                                       rhoair, hpbl, aerdep, oro=oro, ierr=patch_ierr)
         ELSEIF (present(ozone)) THEN
            CALL colm_mpas_set_forcing(patch, pco2m, po2m, us, vs, tair, qair, prc, prl, rain, snow, &
                                       psrf, pbot, sols, soll, solsd, solld, frl, hgt_u, hgt_t, hgt_q, &
                                       rhoair, hpbl, aerdep, ozone=ozone, ierr=patch_ierr)
         ELSE
            CALL colm_mpas_set_forcing(patch, pco2m, po2m, us, vs, tair, qair, prc, prl, rain, snow, &
                                       psrf, pbot, sols, soll, solsd, solld, frl, hgt_u, hgt_t, hgt_q, &
                                       rhoair, hpbl, aerdep, ierr=patch_ierr)
         ENDIF
         IF (patch_ierr /= 0) RETURN
         did_set = .true.
      ENDDO

      IF (did_set) ierr = 0
   END SUBROUTINE colm_mpas_set_element_forcing

	   SUBROUTINE colm_mpas_step(idate, deltim, dolai, doalb, dosst, ierr)
#ifdef GridRiverLakeFlow
	      USE MOD_Grid_RiverLakeFlow, only: grid_riverlake_flow
#endif
	      integer, intent(in) :: idate(3)
	      real(r8), intent(in) :: deltim
	      logical, intent(in) :: dolai, doalb, dosst
      integer, intent(out) :: ierr

      logical :: ready

      CALL colm_mpas_ready(ready)
      ierr = 1
      IF (.not. ready) RETURN

	      CALL CoLMDRIVER(idate, deltim, dolai, doalb, dosst, oroflag)
#ifdef GridRiverLakeFlow
	      CALL grid_riverlake_flow(idate(1), deltim)
#endif
	      colm_mpas_last_idate(:) = idate(:)
	      CALL colm_mpas_write_restart_if_due(idate, deltim, .false., ierr)
	      IF (ierr /= 0) RETURN
	      ierr = 0
	   END SUBROUTINE colm_mpas_step

	   SUBROUTINE colm_mpas_write_restart_if_due(idate, deltim, force, ierr)
	      USE MOD_Namelist, only: DEF_WRST_FREQ
	      USE MOD_TimeManager, only: adj2begin
	      USE MOD_Vars_TimeVariables, only: save_to_restart, WRITE_TimeVariables
	      integer, intent(in) :: idate(3)
	      real(r8), intent(in) :: deltim
	      logical, intent(in) :: force
	      integer, intent(out) :: ierr

	      type(timestamp) :: itstamp
	      integer :: write_idate(3)
	      integer :: write_lc_year
	      logical :: should_write
	      character(len=256) :: wrst_freq

	      ierr = 0
	      IF (.not. colm_mpas_restart_ready) RETURN

	      write_idate(:) = idate(:)
	      CALL adj2begin(write_idate)

	      IF (all(write_idate == colm_mpas_last_restart_idate)) RETURN

	      itstamp%year = idate(1)
	      itstamp%day = idate(2)
	      itstamp%sec = idate(3)

	      wrst_freq = trim(adjustl(DEF_WRST_FREQ))
	      IF (wrst_freq == '' .or. wrst_freq == 'none' .or. wrst_freq == 'NONE') THEN
	         should_write = force .or. colm_mpas_timestamp_reached(itstamp, colm_mpas_etstamp)
	      ELSEIF (force) THEN
	         should_write = .true.
	      ELSE
	         should_write = save_to_restart(idate, deltim, itstamp, colm_mpas_ptstamp, colm_mpas_etstamp)
	      ENDIF
	      IF (.not. should_write) RETURN

#ifdef LULCC
	      IF (write_idate(1) >= 2000) THEN
	         write_lc_year = write_idate(1)
	      ELSE
	         write_lc_year = (write_idate(1) / 5) * 5
	      ENDIF
#else
	      write_lc_year = colm_mpas_lc_year
#endif

	      CALL WRITE_TimeVariables(write_idate, write_lc_year, colm_mpas_casename, colm_mpas_dir_restart)
	      colm_mpas_last_restart_idate(:) = write_idate(:)
	   END SUBROUTINE colm_mpas_write_restart_if_due

	   LOGICAL FUNCTION colm_mpas_timestamp_reached(tstamp, target)
	      type(timestamp), intent(in) :: tstamp
	      type(timestamp), intent(in) :: target

	      colm_mpas_timestamp_reached = .false.
	      IF (tstamp%year > target%year) THEN
	         colm_mpas_timestamp_reached = .true.
	      ELSEIF (tstamp%year == target%year .and. tstamp%day > target%day) THEN
	         colm_mpas_timestamp_reached = .true.
	      ELSEIF (tstamp%year == target%year .and. tstamp%day == target%day .and. tstamp%sec >= target%sec) THEN
	         colm_mpas_timestamp_reached = .true.
	      ENDIF
	   END FUNCTION colm_mpas_timestamp_reached

   SUBROUTINE colm_mpas_get_surface(patch, sensible, latent, evaporation, ground_heat, runoff, &
                                    surface_runoff, subsurface_runoff, skin_temp, t2m, q2m, &
                                    emissivity, roughness, albedo, ierr)
      integer, intent(in) :: patch
      real(r8), intent(out) :: sensible, latent, evaporation, ground_heat, runoff
      real(r8), intent(out) :: surface_runoff, subsurface_runoff, skin_temp, t2m, q2m
      real(r8), intent(out) :: emissivity, roughness, albedo
      integer, intent(out) :: ierr

      ierr = 1
      sensible = spval
      latent = spval
      evaporation = spval
      ground_heat = spval
      runoff = spval
      surface_runoff = spval
      subsurface_runoff = spval
      skin_temp = spval
      t2m = spval
      q2m = spval
      emissivity = spval
      roughness = spval
      albedo = spval

      IF (.not. allocated(fsena)) RETURN
      IF (patch < 1 .or. patch > size(fsena)) RETURN

      sensible = fsena(patch)
      latent = lfevpa(patch)
      evaporation = fevpa(patch)
      ground_heat = fgrnd(patch)
      runoff = rnof(patch)
      surface_runoff = rsur(patch)
      subsurface_runoff = rsub(patch)
      skin_temp = t_grnd(patch)
      t2m = tref(patch)
      q2m = qref(patch)
      emissivity = emis(patch)
      roughness = z0m(patch)
      IF (allocated(alb)) THEN
         ! ponytail: bulk albedo average; replace with band/beam mapping when MPAS consumes CoLM spectral albedo.
         albedo = sum(alb(:,:,patch)) / real(size(alb(:,:,patch)), r8)
      ENDIF
      ierr = 0
   END SUBROUTINE colm_mpas_get_surface

   SUBROUTINE colm_mpas_get_element_surface(element, sensible, latent, evaporation, ground_heat, runoff, &
                                            surface_runoff, subsurface_runoff, skin_temp, t2m, q2m, &
                                            emissivity, roughness, albedo, ierr)
      integer, intent(in) :: element
      real(r8), intent(out) :: sensible, latent, evaporation, ground_heat, runoff
      real(r8), intent(out) :: surface_runoff, subsurface_runoff, skin_temp, t2m, q2m
      real(r8), intent(out) :: emissivity, roughness, albedo
      integer, intent(out) :: ierr

      integer :: patch
      integer :: patch_ierr
      integer :: istt
      integer :: iend
      real(r8) :: wt
      real(r8) :: sumwt
      real(r8) :: albedo_sum
      real(r8) :: albedo_wt
      real(r8) :: patch_sensible, patch_latent, patch_evaporation, patch_ground_heat
      real(r8) :: patch_runoff, patch_surface_runoff, patch_subsurface_runoff
      real(r8) :: patch_skin_temp, patch_t2m, patch_q2m
      real(r8) :: patch_emissivity, patch_roughness, patch_albedo

      ierr = 1
      sensible = spval
      latent = spval
      evaporation = spval
      ground_heat = spval
      runoff = spval
      surface_runoff = spval
      subsurface_runoff = spval
      skin_temp = spval
      t2m = spval
      q2m = spval
      emissivity = spval
      roughness = spval
      albedo = spval

      IF (.not. allocated(elm_patch%substt)) RETURN
      IF (.not. allocated(elm_patch%subfrc)) RETURN
      IF (element < 1 .or. element > size(elm_patch%substt)) RETURN

      istt = elm_patch%substt(element)
      iend = elm_patch%subend(element)
      IF (istt < 1 .or. iend < istt) RETURN

      sensible = 0._r8
      latent = 0._r8
      evaporation = 0._r8
      ground_heat = 0._r8
      runoff = 0._r8
      surface_runoff = 0._r8
      subsurface_runoff = 0._r8
      skin_temp = 0._r8
      t2m = 0._r8
      q2m = 0._r8
      emissivity = 0._r8
      roughness = 0._r8
      albedo_sum = 0._r8
      albedo_wt = 0._r8
      sumwt = 0._r8

      DO patch = istt, iend
         IF (patch < 1 .or. patch > numpatch) CYCLE
         IF (allocated(patchmask)) THEN
            IF (.not. patchmask(patch)) CYCLE
         ENDIF
         wt = elm_patch%subfrc(patch)
         IF (wt <= 0._r8) CYCLE
         CALL colm_mpas_get_surface(patch, patch_sensible, patch_latent, patch_evaporation, patch_ground_heat, &
                                    patch_runoff, patch_surface_runoff, patch_subsurface_runoff, &
                                    patch_skin_temp, patch_t2m, patch_q2m, patch_emissivity, &
                                    patch_roughness, patch_albedo, patch_ierr)
         IF (patch_ierr /= 0) CYCLE

         sumwt = sumwt + wt
         sensible = sensible + wt * patch_sensible
         latent = latent + wt * patch_latent
         evaporation = evaporation + wt * patch_evaporation
         ground_heat = ground_heat + wt * patch_ground_heat
         runoff = runoff + wt * patch_runoff
         surface_runoff = surface_runoff + wt * patch_surface_runoff
         subsurface_runoff = subsurface_runoff + wt * patch_subsurface_runoff
         skin_temp = skin_temp + wt * patch_skin_temp
         t2m = t2m + wt * patch_t2m
         q2m = q2m + wt * patch_q2m
         emissivity = emissivity + wt * patch_emissivity
         roughness = roughness + wt * patch_roughness
         IF (patch_albedo > 0._r8 .and. patch_albedo < 1._r8) THEN
            albedo_sum = albedo_sum + wt * patch_albedo
            albedo_wt = albedo_wt + wt
         ENDIF
      ENDDO

      IF (sumwt <= 0._r8) THEN
         sensible = spval
         latent = spval
         evaporation = spval
         ground_heat = spval
         runoff = spval
         surface_runoff = spval
         subsurface_runoff = spval
         skin_temp = spval
         t2m = spval
         q2m = spval
         emissivity = spval
         roughness = spval
         RETURN
      ENDIF

      sensible = sensible / sumwt
      latent = latent / sumwt
      evaporation = evaporation / sumwt
      ground_heat = ground_heat / sumwt
      runoff = runoff / sumwt
      surface_runoff = surface_runoff / sumwt
      subsurface_runoff = subsurface_runoff / sumwt
      skin_temp = skin_temp / sumwt
      t2m = t2m / sumwt
      q2m = q2m / sumwt
      emissivity = emissivity / sumwt
      roughness = roughness / sumwt
      IF (albedo_wt > 0._r8) albedo = albedo_sum / albedo_wt
      ierr = 0
   END SUBROUTINE colm_mpas_get_element_surface

END MODULE MOD_CoLM_MPAS_Interface
