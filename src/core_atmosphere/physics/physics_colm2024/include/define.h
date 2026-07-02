! 1. Spatial structure:
!    Select one of the following options.
#define GRIDBASED
#undef CATCHMENT
#undef UNSTRUCTURED
#undef SinglePoint

! 2. Land subgrid type classification:
!    Select one of the following options.
#undef LULC_USGS
#undef LULC_IGBP
#undef LULC_IGBP_PFT
#define LULC_IGBP_PC

! 2.1 3D Urban model (put it temporarily here):
#undef URBAN_MODEL
!    Dependence:  only LULC_IGBP subgrid type for
!    single point URBAN_MODEL right now.
#if (defined URBAN_MODEL && defined SinglePoint)
#define LULC_IGBP
#undef LULC_USGS
#undef LULC_IGBP_PFT
#undef LULC_IGBP_PC
#endif

! 3. If defined, debug information is output.
#define CoLMDEBUG
! 3.1 If defined, range of variables is checked.
#define RangeCheck
! 3.1 If defined, surface data in vector is mapped to gridded data for checking.
#undef SrfdataDiag

! 4. CoLM uses MPI collectives through MPAS-owned communicators.  The legacy
!    CoLM process-pool decomposition is disabled in MOD_SPMD_Task.
#define USEMPI

! 5. Hydrological process options.
! 5.1 Two soil hydraulic models can be used.
#undef   Campbell_SOIL_MODEL
#define  vanGenuchten_Mualem_SOIL_MODEL
! 5.2 If defined, lateral flow is modeled.
#define CatchLateralFlow
!    Conflicts :
#ifndef CATCHMENT
#undef CatchLateralFlow
#endif

! 6. CaMa-Flood is not built in the MPAS-embedded CoLM package.
#undef CaMa_Flood
#if (defined SinglePoint)
#undef CaMa_Flood
#endif
#ifndef USEMPI
#undef CaMa_Flood
#endif

#define GridRiverLakeFlow
!    Conflicts :
#if (defined CATCHMENT || defined SinglePoint)
#undef GridRiverLakeFlow
#endif

#undef GridRiverLakeSediment
#if (!defined GridRiverLakeFlow)
#undef GridRiverLakeSediment
#endif

! 7. If defined, BGC model is used.
#undef BGC

!    Conflicts :  only used when LULC_IGBP_PFT is defined.
#ifndef LULC_IGBP_PFT
#ifndef LULC_IGBP_PC
#undef BGC
#endif
#endif
! 7.1 If defined, CROP model is used
#undef CROP
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef CROP
#endif

! 8. If defined, open Land use and land cover change mode.
#undef LULCC

! 9. If defined, data assimilation is used.
#undef DataAssimilation
#if (defined DataAssimilation)
#define LULC_IGBP
#undef LULC_USGS
#undef LULC_IGBP_PFT
#undef LULC_IGBP_PC
#endif

! 10. Interface to AI model.
#undef USESplitAI

! 11. External lake models.
#undef EXTERNAL_LAKE

! 12. Hyperspectral scheme.
#define HYPERSPECTRAL

! MPAS embeds CoLM as a local land-surface physics package on each MPAS rank.
! Reusing the legacy CoLM process split inside MPAS leaves some
! ranks without patch forcing/flux arrays, so the embedded library keeps CoLM
! land state on MPAS-owned cells. MPAS also supplies broadband shortwave forcing.
#define MPAS_EMBEDDED_COLM
#ifdef MPAS_EMBEDDED_COLM
#ifndef USEMPI
#define USEMPI
#endif
#undef HYPERSPECTRAL
#undef CoLMDEBUG
! CoLM range diagnostics still assume standalone IO/worker roles.
#undef RangeCheck
#endif

#undef COLM_PARALLEL
#if defined(USEMPI)
#define COLM_PARALLEL
#endif
