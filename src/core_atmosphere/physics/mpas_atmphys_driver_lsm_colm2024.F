! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_lsm_colm2024
 use mpas_kind_types
 use mpas_derived_types, only: dm_info
 use mpas_pool_routines

 use mpas_atmphys_constants, only: R_d, R_v
 use mpas_atmphys_manager, only: year, curr_julday
 use mpas_atmphys_utilities, only: physics_error_fatal
 use mpas_atmphys_vars, only: dt_pbl, xice_threshold

#ifdef MPAS_COLM2024
	 use MOD_Precision, only: colm_r8 => r8
	 use MOD_CoLM_MPAS_Interface, only: colm_mpas_initialize_from_namelist, colm_mpas_ready, &
	                                    colm_mpas_set_element_forcing, &
	                                    colm_mpas_step, colm_mpas_get_element_surface, &
	                                    colm_mpas_finalize
#endif

 implicit none
 private
 public :: init_lsm_colm2024, driver_lsm_colm2024, finalize_lsm_colm2024

 logical, save :: colm_first_step = .true.
#ifdef MPAS_COLM2024
 integer, allocatable, save :: colm_element_for_cell(:)
#endif

 contains

!=================================================================================================================
 subroutine init_lsm_colm2024(dminfo,configs,mesh,diag_physics,sfc_input)
!=================================================================================================================

 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: sfc_input
#ifdef MPAS_COLM2024
 logical :: colm_ready
 integer :: colm_ierr
 integer :: colm_mpas_comm
	 integer :: colm_patch_count
 character(len=StrKIND),pointer:: config_colm2024_namelist
 integer,pointer:: nCellsSolve
 integer,dimension(:),pointer:: indexToCellID
 real(kind=RKIND),dimension(:),pointer:: latCell, lonCell
#endif

#ifdef MPAS_COLM2024
 call mpas_pool_get_config(configs,'config_colm2024_namelist',config_colm2024_namelist)
 if(trim(config_colm2024_namelist) == 'none' .or. trim(config_colm2024_namelist) == '') then
    call physics_error_fatal("config_lsm_scheme='sf_colm2024' requires config_colm2024_namelist.")
 endif

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)
 call mpas_pool_get_array(mesh,'indexToCellID',indexToCellID)
	 if(allocated(colm_element_for_cell)) deallocate(colm_element_for_cell)
	 allocate(colm_element_for_cell(size(latCell)))
	 colm_element_for_cell(:) = 0

#ifdef MPAS_USE_MPI_F08
 colm_mpas_comm = dminfo % comm % mpi_val
#else
 colm_mpas_comm = dminfo % comm
#endif

 call colm_mpas_initialize_from_namelist(trim(config_colm2024_namelist), colm_ierr, &
                                        mpas_comm=colm_mpas_comm, &
                                        mpas_lon_rad=real(lonCell(1:nCellsSolve), colm_r8), &
                                        mpas_lat_rad=real(latCell(1:nCellsSolve), colm_r8), &
                                        mpas_cell_id=indexToCellID(1:nCellsSolve), &
                                        n_mpas_cells=nCellsSolve, &
                                        cell_to_element=colm_element_for_cell(1:nCellsSolve))
 if(colm_ierr /= 0) call physics_error_fatal("CoLM2024 initialization from config_colm2024_namelist failed.")

 call colm_mpas_ready(colm_ready, colm_patch_count)
 if(.not. colm_ready) call colm2024_adapter_not_ready

 colm_first_step = .true.
#else
 call colm2024_adapter_not_built
#endif

 end subroutine init_lsm_colm2024

!=================================================================================================================
 subroutine driver_lsm_colm2024(configs,mesh,state,time_lev,diag,diag_physics,sfc_input,itimestep,its,ite)
!=================================================================================================================

 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: sfc_input
 integer,intent(in):: time_lev
 integer,intent(in):: itimestep
 integer,intent(in):: its,ite
#ifdef MPAS_COLM2024
 logical :: colm_ready
 integer :: colm_ierr
 integer :: colm_patch_count
 integer :: colm_idate(3)
#endif

#ifdef MPAS_COLM2024
 call colm_mpas_ready(colm_ready, colm_patch_count)
 if(.not. colm_ready) call colm2024_adapter_not_ready
	 if(.not. allocated(colm_element_for_cell)) call physics_error_fatal("CoLM2024 MPAS cell-to-element map is not initialized.")
	 if(size(colm_element_for_cell) < ite) call physics_error_fatal("CoLM2024 MPAS cell-to-element map is smaller " // &
	                                                                "than the MPAS cell range.")

	 call colm2024_push_mpas_forcing(configs, mesh, state, time_lev, diag, diag_physics, sfc_input, itimestep, its, ite)

 colm_idate(1) = year
 colm_idate(2) = int(curr_julday) + 1
 colm_idate(3) = nint((curr_julday - real(int(curr_julday), RKIND)) * 86400._RKIND)
 call colm_mpas_step(colm_idate, real(dt_pbl, colm_r8), colm_first_step, .true., .false., colm_ierr)
 if(colm_ierr /= 0) call physics_error_fatal("CoLM2024 driver returned before completing a land-surface step.")

 call colm2024_pull_mpas_surface(diag_physics, sfc_input, real(dt_pbl, colm_r8), its, ite)
 colm_first_step = .false.
#else
 call colm2024_adapter_not_built
#endif

 end subroutine driver_lsm_colm2024

!=================================================================================================================
	 subroutine finalize_lsm_colm2024()
	!=================================================================================================================

#ifdef MPAS_COLM2024
	 integer :: colm_ierr

	 call colm_mpas_finalize(colm_ierr)
	 if(colm_ierr /= 0) call physics_error_fatal("CoLM2024 finalization failed.")
#endif

	 end subroutine finalize_lsm_colm2024

#ifdef MPAS_COLM2024
!=================================================================================================================
 subroutine colm2024_push_mpas_forcing(configs,mesh,state,time_lev,diag,diag_physics,sfc_input,itimestep,its,ite)
!=================================================================================================================

	 type(mpas_pool_type),intent(in):: configs
	 type(mpas_pool_type),intent(in):: mesh
	 type(mpas_pool_type),intent(in):: state
	 type(mpas_pool_type),intent(in):: diag
	 type(mpas_pool_type),intent(in):: diag_physics
	 type(mpas_pool_type),intent(in):: sfc_input
	 integer,intent(in):: time_lev
	 integer,intent(in):: itimestep
	 integer,intent(in):: its,ite

	 character(len=StrKIND),pointer:: config_convection_scheme, config_microp_scheme
	 integer,pointer:: index_qv
	 integer:: i, k
	 integer:: ierr
	 real(kind=colm_r8):: aerdep(14)
	 real(kind=colm_r8):: dt, grid_snow, hgt, oro, pbot, psrf, prc, prl, qair, rain, rhoair, snow, &
	                      solld, soll, solsd, sols, spectral_sw, swdir, swdif, tair, uair, vair
	 real(kind=RKIND):: qv_cell
	 real(kind=RKIND),dimension(:),pointer:: glw, swddir, swddif, swvisdir, swvisdif, swnirdir, swnirdif
	 real(kind=RKIND),dimension(:),pointer:: surface_pressure, sr, xice, xland
	 real(kind=RKIND),dimension(:),pointer:: raincv, rainncv, snowncv, graupelncv
	 real(kind=RKIND),dimension(:,:),pointer:: exner, pressure_b, pressure_p, theta_m, u, v, zgrid
	 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 call mpas_pool_get_config(configs,'config_convection_scheme',config_convection_scheme)
	 call mpas_pool_get_config(configs,'config_microp_scheme'    ,config_microp_scheme    )

	 call mpas_pool_get_dimension(state,'index_qv',index_qv)
	 call mpas_pool_get_array(mesh,'zgrid',zgrid)
	 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
	 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
	 call mpas_pool_get_array(diag,'exner'                 ,exner     )
	 call mpas_pool_get_array(diag,'pressure_base'         ,pressure_b)
	 call mpas_pool_get_array(diag,'pressure_p'            ,pressure_p)
	 call mpas_pool_get_array(diag,'uReconstructZonal'     ,u         )
	 call mpas_pool_get_array(diag,'uReconstructMeridional',v         )
	 call mpas_pool_get_array(diag,'surface_pressure',surface_pressure)
	 call mpas_pool_get_array(diag_physics,'glw'   ,glw   )
	 call mpas_pool_get_array(diag_physics,'swddir',swddir)
	 call mpas_pool_get_array(diag_physics,'swddif',swddif)
	 call mpas_pool_get_array(diag_physics,'swvisdir',swvisdir)
	 call mpas_pool_get_array(diag_physics,'swvisdif',swvisdif)
	 call mpas_pool_get_array(diag_physics,'swnirdir',swnirdir)
	 call mpas_pool_get_array(diag_physics,'swnirdif',swnirdif)
	 call mpas_pool_get_array(diag_physics,'raincv'    ,raincv    )
	 call mpas_pool_get_array(diag_physics,'rainncv'   ,rainncv   )
	 call mpas_pool_get_array(diag_physics,'snowncv'   ,snowncv   )
	 call mpas_pool_get_array(diag_physics,'graupelncv',graupelncv)
	 call mpas_pool_get_array(diag_physics,'sr'        ,sr        )
	 call mpas_pool_get_array(sfc_input,'xland',xland)
	 call mpas_pool_get_array(sfc_input,'xice' ,xice )

 dt = max(real(dt_pbl, colm_r8), 1._colm_r8)
	 aerdep(:) = 0._colm_r8

	 k = lbound(theta_m, 1)
	 do i = its,ite
	    if(colm_element_for_cell(i) <= 0) call physics_error_fatal("CoLM2024 MPAS cell-to-element map has an invalid element index.")

	    psrf = real(surface_pressure(i), colm_r8)
	    pbot = real(pressure_b(k,i) + pressure_p(k,i), colm_r8)
	    qv_cell = max(0._RKIND, scalars(index_qv,k,i))
	    tair = real((theta_m(k,i) / (1._RKIND + R_v / R_d * qv_cell)) * exner(k,i), colm_r8)
	    qair = real(qv_cell / (1._RKIND + qv_cell), colm_r8)
	    uair = real(u(k,i), colm_r8)
	    vair = real(v(k,i), colm_r8)
	    hgt = 2._colm_r8
	    if(ubound(zgrid, 1) >= k + 1) hgt = max(2._colm_r8, real(0.5_RKIND * (zgrid(k+1,i) - zgrid(k,i)), colm_r8))
	    rhoair = max(0.1_colm_r8, pbot / (R_d * tair * (1._colm_r8 + 0.61_colm_r8 * qair)))
	    prc = 0._colm_r8
	    prl = 0._colm_r8
	    rain = 0._colm_r8
		    snow = 0._colm_r8

		    if(config_convection_scheme .ne. 'off') then
		       prc = max(0._colm_r8, real(raincv(i), colm_r8) / dt)
		       rain = rain + prc
		    endif
		    if(config_microp_scheme .ne. 'off') then
		       prl = max(0._colm_r8, real(rainncv(i), colm_r8) / dt)
		       grid_snow = max(0._colm_r8, real(snowncv(i) + graupelncv(i), colm_r8) / dt)
		       if(grid_snow <= 0._colm_r8 .and. sr(i) > 0._RKIND) then
		          grid_snow = prl * min(1._colm_r8, real(sr(i), colm_r8))
		       endif
		       snow = min(prl, grid_snow)
		       rain = rain + max(0._colm_r8, prl - snow)
		    endif

	    swdir = max(0._colm_r8, real(swddir(i), colm_r8))
	    swdif = max(0._colm_r8, real(swddif(i), colm_r8))
	    sols = max(0._colm_r8, real(swvisdir(i), colm_r8))
	    solsd = max(0._colm_r8, real(swvisdif(i), colm_r8))
	    soll = max(0._colm_r8, real(swnirdir(i), colm_r8))
	    solld = max(0._colm_r8, real(swnirdif(i), colm_r8))
	    spectral_sw = sols + solsd + soll + solld
	    if(spectral_sw <= 0._colm_r8 .and. swdir + swdif > 0._colm_r8) then
	       sols = 0.5_colm_r8 * swdir
	       soll = 0.5_colm_r8 * swdir
	       solsd = 0.5_colm_r8 * swdif
	       solld = 0.5_colm_r8 * swdif
	    endif
	    oro = 1._colm_r8
	    if(xland(i) >= 1.5_RKIND) oro = 0._colm_r8
	    if(xice(i) >= xice_threshold) oro = 2._colm_r8

	    call colm_mpas_set_element_forcing(colm_element_for_cell(i), psrf * 420.e-6_colm_r8, &
	                                       psrf * 0.209_colm_r8, uair, vair, tair, qair, &
	                                       prc, prl, rain, snow, psrf, pbot, &
	                                       sols, soll, solsd, solld, &
	                                       real(glw(i), colm_r8), hgt, hgt, hgt, rhoair, 1000._colm_r8, &
	                                       aerdep, oro=oro, ierr=ierr)
	    if(ierr /= 0) call physics_error_fatal("Failed to pass MPAS forcing into CoLM2024.")
	 enddo

 end subroutine colm2024_push_mpas_forcing

!=================================================================================================================
 subroutine colm2024_pull_mpas_surface(diag_physics,sfc_input,dt,its,ite)
!=================================================================================================================

 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: sfc_input
 real(kind=colm_r8),intent(in):: dt
 integer,intent(in):: its,ite

 integer:: i
 integer:: ierr
 real(kind=colm_r8):: sensible, latent, evaporation, ground_heat, runoff
 real(kind=colm_r8):: surface_runoff, subsurface_runoff, skin_temp, t2m_colm, q2m_colm
 real(kind=colm_r8):: emissivity, roughness, albedo
 real(kind=RKIND),dimension(:),pointer:: hfx, lh, qfx, grdflx, sfcrunoff, udrunoff
 real(kind=RKIND),dimension(:),pointer:: sfc_albedo, sfc_emiss, z0, znt, qsfc, t2m, q2
 real(kind=RKIND),dimension(:),pointer:: skintemp

 call mpas_pool_get_array(diag_physics,'hfx'       ,hfx       )
 call mpas_pool_get_array(diag_physics,'lh'        ,lh        )
 call mpas_pool_get_array(diag_physics,'qfx'       ,qfx       )
 call mpas_pool_get_array(diag_physics,'grdflx'    ,grdflx    )
 call mpas_pool_get_array(diag_physics,'sfcrunoff' ,sfcrunoff )
 call mpas_pool_get_array(diag_physics,'udrunoff'  ,udrunoff  )
 call mpas_pool_get_array(diag_physics,'sfc_albedo',sfc_albedo)
 call mpas_pool_get_array(diag_physics,'sfc_emiss' ,sfc_emiss )
 call mpas_pool_get_array(diag_physics,'z0'        ,z0        )
 call mpas_pool_get_array(diag_physics,'znt'       ,znt       )
 call mpas_pool_get_array(diag_physics,'qsfc'      ,qsfc      )
 call mpas_pool_get_array(diag_physics,'t2m'       ,t2m       )
 call mpas_pool_get_array(diag_physics,'q2'        ,q2        )
 call mpas_pool_get_array(sfc_input,'skintemp',skintemp)

	 do i = its,ite
	    if(colm_element_for_cell(i) <= 0) call physics_error_fatal("CoLM2024 MPAS cell-to-element map has an invalid element index.")

	    call colm_mpas_get_element_surface(colm_element_for_cell(i), sensible, latent, evaporation, ground_heat, runoff, &
	                                       surface_runoff, subsurface_runoff, skin_temp, t2m_colm, q2m_colm, &
	                                       emissivity, roughness, albedo, ierr)
    if(ierr /= 0) call physics_error_fatal("Failed to retrieve CoLM2024 surface fluxes.")

		    hfx(i) = real(sensible, RKIND)
		    lh(i) = real(latent, RKIND)
		    qfx(i) = real(evaporation, RKIND)
		    grdflx(i) = real(ground_heat, RKIND)
		    sfcrunoff(i) = real(surface_runoff * dt, RKIND)
		    udrunoff(i) = real(subsurface_runoff * dt, RKIND)
		    skintemp(i) = real(skin_temp, RKIND)
		    t2m(i) = real(t2m_colm, RKIND)
		    qsfc(i) = real(q2m_colm, RKIND)
	    if(q2m_colm > 0._colm_r8 .and. q2m_colm < 1._colm_r8) then
	       q2(i) = real(q2m_colm / (1._colm_r8 - q2m_colm), RKIND)
	    else
	       q2(i) = real(q2m_colm, RKIND)
	    endif
		    sfc_emiss(i) = real(emissivity, RKIND)
		    z0(i) = real(roughness, RKIND)
		    znt(i) = real(roughness, RKIND)
		    if(albedo > 0._colm_r8 .and. albedo < 1._colm_r8) sfc_albedo(i) = real(albedo, RKIND)
 enddo

 end subroutine colm2024_pull_mpas_surface
#endif

!=================================================================================================================
 subroutine colm2024_adapter_not_ready()
!=================================================================================================================

 call physics_error_fatal("CoLM2024 is linked, but CoLM landpatch/restart state has not been initialized " // &
                          "for this MPAS run.")

 end subroutine colm2024_adapter_not_ready

!=================================================================================================================
 subroutine colm2024_adapter_not_built()
!=================================================================================================================

 call physics_error_fatal("CoLM2024 requires building MPAS with COLM2024=true.")

 end subroutine colm2024_adapter_not_built

!=================================================================================================================
 end module mpas_atmphys_driver_lsm_colm2024
!=================================================================================================================
