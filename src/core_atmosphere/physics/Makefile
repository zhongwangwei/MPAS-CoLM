.SUFFIXES: .F .o

ifeq ($(CORE),atmosphere)
COREDEF = -Dmpas
endif

ifneq (,$(findstring -DMPAS_COLM2024,$(CPPFLAGS)))
COLM2024_ENABLED = true
COLM2024_CORE = core_physics_colm2024_embedded
COLM2024_BUILD_DIR = physics_colm2024/.mpas_bld
COLM2024_HEADER = physics_colm2024/include/define.h
COLM2024_STAMP = $(COLM2024_BUILD_DIR)/.built
COLM2024_MOD_CMD ?= -J
NF_CONFIG ?= nf-config
COLM2024_NETCDF_INC ?= $(shell $(NF_CONFIG) --includedir 2>/dev/null || nc-config --includedir 2>/dev/null)
COLM2024_INCLUDES = -I./physics_colm2024/include -I./$(COLM2024_BUILD_DIR)
ifneq ($(strip $(COLM2024_NETCDF_INC)),)
COLM2024_INCLUDES += -I$(COLM2024_NETCDF_INC)
endif
ifneq ($(strip $(PNETCDF)),)
COLM2024_INCLUDES += -I$(PNETCDF)/include
endif
ifneq ($(strip $(PNETCDF_INC)),)
COLM2024_INCLUDES += -I$(PNETCDF_INC)
endif
COLM2024_FOPTS ?= $(CPPFLAGS) $(COREDEF) $(FFLAGS) -fdefault-real-8 -ffree-form -cpp -ffree-line-length-0 -fallow-argument-mismatch
COLM2024_SOURCE_DIRS = \
	physics_colm2024/share \
	physics_colm2024/main \
	physics_colm2024/main/HYDRO \
	physics_colm2024/main/BGC \
	physics_colm2024/main/URBAN \
	physics_colm2024/main/DA \
	physics_colm2024/main/ParaOpt \
	physics_colm2024/drivers/mpas
COLM2024_SOURCES = $(foreach dir,$(COLM2024_SOURCE_DIRS),$(wildcard $(dir)/*.F90))
endif

all:
	./../tools/manage_externals/checkout_externals --externals ./../Externals.cfg
	$(RM) libphys.a
	$(MAKE) lookup_tables $(COLM2024_CORE) core_physics_init core_physics_mmm core_UGWP core_physics_wrf core_physics_noahmp core_physics

dummy:
	echo "****** compiling physics ******"

OBJS_init = \
	ccpp_kind_types.o        \
	mpas_atmphys_constants.o \
	mpas_atmphys_date_time.o \
	mpas_atmphys_functions.o \
	mpas_atmphys_utilities.o

OBJS = \
	mpas_atmphys_camrad_init.o         \
	mpas_atmphys_control.o             \
	mpas_atmphys_driver.o              \
	mpas_atmphys_driver_cloudiness.o   \
	mpas_atmphys_driver_convection.o   \
	mpas_atmphys_driver_gwdo.o         \
	mpas_atmphys_driver_lsm.o          \
	mpas_atmphys_driver_lsm_noahmp.o   \
	mpas_atmphys_driver_lsm_colm2024.o \
	mpas_atmphys_driver_microphysics.o \
	mpas_atmphys_driver_oml.o          \
	mpas_atmphys_driver_pbl.o          \
	mpas_atmphys_driver_radiation_lw.o \
	mpas_atmphys_driver_radiation_sw.o \
	mpas_atmphys_driver_seaice.o       \
	mpas_atmphys_driver_sfclayer.o     \
	mpas_atmphys_finalize.o            \
	mpas_atmphys_init.o                \
	mpas_atmphys_init_microphysics.o   \
	mpas_atmphys_interface.o           \
	mpas_atmphys_landuse.o             \
	mpas_atmphys_lsm_noahinit.o        \
	mpas_atmphys_lsm_noahmpinit.o      \
	mpas_atmphys_lsm_noahmpfinalize.o  \
	mpas_atmphys_lsm_shared.o          \
	mpas_atmphys_manager.o             \
	mpas_atmphys_o3climatology.o       \
	mpas_atmphys_packages.o            \
	mpas_atmphys_rrtmg_lwinit.o        \
	mpas_atmphys_rrtmg_swinit.o        \
	mpas_atmphys_sfc_diagnostics.o     \
	mpas_atmphys_todynamics.o          \
	mpas_atmphys_update_surface.o      \
	mpas_atmphys_update.o              \
	mpas_atmphys_vars.o

lookup_tables:
	./checkout_data_files.sh

core_physics_mmm: core_physics_init
	(cd physics_mmm; $(MAKE) -f Makefile.mpas all)

core_UGWP: core_physics_init
	(cd physics_noaa/UGWP; $(MAKE) all)

core_physics_wrf: core_physics_init core_physics_mmm core_UGWP
	(cd physics_wrf; $(MAKE) all COREDEF="$(COREDEF)")

core_physics_noahmp:
	(cd physics_noahmp/utility; $(MAKE) all COREDEF="$(COREDEF)")
	(cd physics_noahmp/src; $(MAKE) all COREDEF="$(COREDEF)")
	(cd physics_noahmp/drivers/mpas; $(MAKE) all COREDEF="$(COREDEF)")

ifeq "$(COLM2024_ENABLED)" "true"
COLM2024_OBJS_SHARED_NAMES = \
	MOD_Precision.o \
	MOD_SPMD_Task.o \
	MOD_Namelist.o \
	MOD_Vars_Global.o \
	MOD_Const_Physical.o \
	MOD_Const_LC.o \
	MOD_Utils.o \
	MOD_IncompleteGamma.o \
	MOD_UserDefFun.o \
	MOD_TimeManager.o \
	MOD_Const_PFT.o \
	MOD_NetCDFSerial.o \
	MOD_Block.o \
	MOD_Grid.o \
	MOD_Pixel.o \
	MOD_DataType.o \
	MOD_NetCDFPoint.o \
	MOD_NetCDFBlock.o \
	MOD_CatchmentDataReadin.o \
	MOD_5x5DataReadin.o \
	MOD_Mesh.o \
	MOD_Pixelset.o \
	MOD_NetCDFVector.o \
	MOD_RangeCheck.o \
	MOD_SpatialMapping.o \
	MOD_ComputePushData.o \
	MOD_AggregationRequestData.o \
	MOD_PixelsetShared.o \
	MOD_LandElm.o \
	MOD_LandHRU.o \
	MOD_LandPatch.o \
	MOD_Land2mWMO.o \
	MOD_LandCrop.o \
	MOD_LandPFT.o \
	MOD_LandUrban.o \
	MOD_Urban_Const_LCZ.o \
	MOD_SingleSrfdata.o \
	MOD_SrfdataRestart.o \
	MOD_ElmVector.o \
	MOD_HRUVector.o \
	MOD_MeshFilter.o

COLM2024_OBJS_BASIC_NAMES = \
	MOD_Vector_ReadWrite.o \
	MOD_dataSpec_PDB.o \
	MOD_tav_abs.o \
	MOD_prospect_DB.o \
	MOD_Catch_BasinNetwork.o \
	MOD_Catch_Vars_TimeVariables.o \
	MOD_Catch_Vars_1DFluxes.o \
	MOD_Grid_RiverLakeNetwork.o \
	MOD_Grid_Reservoir.o \
	MOD_Grid_RiverLakeTimeVars.o \
	MOD_BGC_Vars_1DFluxes.o \
	MOD_BGC_Vars_1DPFTFluxes.o \
	MOD_BGC_Vars_PFTimeVariables.o \
	MOD_BGC_Vars_TimeInvariants.o \
	MOD_BGC_Vars_TimeVariables.o \
	MOD_Urban_Vars_1DFluxes.o \
	MOD_Urban_Vars_TimeVariables.o \
	MOD_Urban_Vars_TimeInvariants.o \
	MOD_DA_Vars_1DFluxes.o \
	MOD_Vars_TimeInvariants.o \
	MOD_DA_Vars_TimeVariables.o \
	MOD_Vars_TimeVariables.o \
	MOD_Vars_1DPFTFluxes.o \
	MOD_Vars_1DFluxes.o \
	MOD_Vars_1DForcing.o \
	MOD_Hydro_SoilFunction.o \
	MOD_Hydro_SoilWater.o \
	MOD_Eroot.o \
	MOD_Qsadv.o \
	MOD_LAIEmpirical.o \
	MOD_LAIReadin.o \
	MOD_CropReadin.o \
	MOD_NitrifData.o \
	MOD_NdepData.o \
	MOD_FireData.o \
	MOD_OrbCoszen.o \
	MOD_OrbCosazi.o \
	MOD_HighRes_Parameters.o \
	MOD_3DCanopyRadiation.o \
	MOD_Aerosol.o \
	MOD_SnowSnicar.o \
	MOD_Albedo.o \
	MOD_SnowSnicar_HiRes.o \
	MOD_Albedo_HiRes.o \
	MOD_SnowFraction.o \
	MOD_Urban_LAIReadin.o \
	MOD_Urban_Shortwave.o \
	MOD_Urban_Albedo.o \
	MOD_MonthlyinSituCO2MaunaLoa.o \
	MOD_BGC_CNSummary.o \
	MOD_ElementNeighbour.o \
	MOD_Catch_HillslopeNetwork.o \
	MOD_Catch_RiverLakeNetwork.o \
	MOD_Catch_Reservoir.o \
	MOD_VicParaReadin.o

COLM2024_OBJS_MAIN_NAMES = \
	MOD_Catch_HillslopeFlow.o \
	MOD_Catch_SubsurfaceFlow.o \
	MOD_Catch_RiverLakeFlow.o \
	MOD_Catch_Hist.o \
	MOD_Catch_WriteParameters.o \
	MOD_BGC_CNCStateUpdate1.o \
	MOD_BGC_CNCStateUpdate2.o \
	MOD_BGC_CNCStateUpdate3.o \
	MOD_BGC_CNNStateUpdate1.o \
	MOD_BGC_CNNStateUpdate2.o \
	MOD_BGC_CNNStateUpdate3.o \
	MOD_BGC_Soil_BiogeochemNStateUpdate1.o \
	MOD_BGC_Soil_BiogeochemNitrifDenitrif.o \
	MOD_BGC_Soil_BiogeochemCompetition.o \
	MOD_BGC_Soil_BiogeochemDecompCascadeBGC.o \
	MOD_BGC_Soil_BiogeochemDecomp.o \
	MOD_BGC_Soil_BiogeochemLittVertTransp.o \
	MOD_BGC_Soil_BiogeochemNLeaching.o \
	MOD_BGC_Soil_BiogeochemPotential.o \
	MOD_BGC_Soil_BiogeochemVerticalProfile.o \
	MOD_BGC_Veg_CNGapMortality.o \
	MOD_BGC_Veg_CNGResp.o \
	MOD_BGC_Veg_CNMResp.o \
	MOD_BGC_Daylength.o \
	MOD_BGC_Veg_CNPhenology.o \
	MOD_BGC_Veg_NutrientCompetition.o \
	MOD_BGC_Veg_CNVegStructUpdate.o \
	MOD_BGC_CNAnnualUpdate.o \
	MOD_BGC_CNZeroFluxes.o \
	MOD_BGC_CNBalanceCheck.o \
	MOD_BGC_CNSASU.o \
	MOD_BGC_Veg_CNNDynamics.o \
	MOD_BGC_Veg_CNFireBase.o \
	MOD_BGC_Veg_CNFireLi2016.o \
	MOD_Vars_2DForcing.o \
	MOD_UserSpecifiedForcing.o \
	MOD_ForcingDownscaling.o \
	MOD_Forcing.o \
	MOD_DA_TWS.o \
	MOD_DA_Const.o \
	MOD_DA_RTM.o \
	MOD_DA_EnKF.o \
	MOD_DA_SM.o \
	MOD_DA_Ensemble.o \
	MOD_DA_Main.o \
	MOD_Opt_Baseflow.o \
	MOD_ParameterOptimization.o \
	MOD_AssimStomataConductance.o \
	MOD_PlantHydraulic.o \
	MOD_FrictionVelocity.o \
	MOD_TurbulenceLEddy.o \
	MOD_Ozone.o \
	MOD_CanopyLayerProfile.o \
	MOD_LeafTemperature.o \
	MOD_LeafTemperaturePC.o \
	MOD_SoilThermalParameters.o \
	MOD_Hydro_VIC_Variables.o \
	MOD_Hydro_VIC.o \
	MOD_Runoff.o \
	MOD_SoilSnowHydrology.o \
	MOD_SnowLayersCombineDivide.o \
	MOD_PhaseChange.o \
	MOD_Glacier.o \
	MOD_Lake.o \
	MOD_SimpleOcean.o \
	MOD_GroundFluxes.o \
	MOD_GroundTemperature.o \
	MOD_LeafInterception.o \
	MOD_NetSolar.o \
	MOD_NetSolar_Hyper.o \
	MOD_WetBulb.o \
	MOD_RainSnowTemp.o \
	MOD_SoilSurfaceResistance.o \
	MOD_NewSnow.o \
	MOD_Thermal.o \
	MOD_Vars_1DAccFluxes.o \
	MOD_Irrigation.o \
	MOD_BGC_driver.o \
	MOD_HistWriteBack.o \
	MOD_HistGridded.o \
	MOD_HistVector.o \
	MOD_HistSingle.o \
	MOD_Grid_RiverLakeHist.o \
	MOD_Hist.o \
	MOD_CheckEquilibrium.o \
	MOD_LightningData.o \
	MOD_Catch_LateralFlow.o \
	MOD_Grid_RiverLakeFlow.o \
	MOD_Urban_Longwave.o \
	MOD_Urban_NetSolar.o \
	MOD_Urban_Flux.o \
	MOD_Urban_GroundFlux.o \
	MOD_Urban_RoofFlux.o \
	MOD_Urban_RoofTemperature.o \
	MOD_Urban_WallTemperature.o \
	MOD_Urban_PerviousTemperature.o \
	MOD_Urban_ImperviousTemperature.o \
	MOD_Urban_Hydrology.o \
	MOD_Urban_BEM.o \
	MOD_Urban_LUCY.o \
	MOD_Urban_Thermal.o \
	CoLMMAIN_Urban.o \
	CoLMDRIVER.o \
	MOD_CoLM_MPAS_Interface.o \
	CoLMMAIN.o

COLM2024_OBJ_NAMES = $(COLM2024_OBJS_SHARED_NAMES) $(COLM2024_OBJS_BASIC_NAMES) $(COLM2024_OBJS_MAIN_NAMES)
COLM2024_OBJS = $(addprefix $(COLM2024_BUILD_DIR)/,$(COLM2024_OBJ_NAMES))

.PHONY: core_physics_colm2024_embedded
core_physics_colm2024_embedded: $(COLM2024_STAMP)

$(COLM2024_BUILD_DIR):
	mkdir -p $(COLM2024_BUILD_DIR)

$(COLM2024_STAMP): $(COLM2024_HEADER) $(COLM2024_SOURCES) | $(COLM2024_BUILD_DIR)
	@echo 'compiling embedded CoLM2024 physics >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	@set -e; \
	for obj in $(COLM2024_OBJ_NAMES); do \
		src=""; \
		stem="$${obj%.o}"; \
		for dir in $(COLM2024_SOURCE_DIRS); do \
			if [ -f "$$dir/$$stem.F90" ]; then src="$$dir/$$stem.F90"; break; fi; \
		done; \
		if [ -z "$$src" ]; then \
			echo "CoLM2024 source not found for $$obj"; \
			exit 1; \
		fi; \
		echo "  $$src"; \
		$(FC) -c $(COLM2024_FOPTS) $(COLM2024_INCLUDES) -o "$(COLM2024_BUILD_DIR)/$$obj" "$$src" $(COLM2024_MOD_CMD)$(COLM2024_BUILD_DIR); \
	done; \
	touch $@
endif

core_physics_init: $(OBJS_init)

core_physics: core_physics_wrf core_physics_noahmp $(COLM2024_CORE)
	($(MAKE) phys_interface COREDEF="$(COREDEF)")
	ar -ru libphys.a $(OBJS_init) $(OBJS)
	($(MAKE) -C ./physics_mmm -f Makefile.mpas physics_mmm_lib)
	($(MAKE) -C ./physics_wrf physics_wrf_lib)
	($(MAKE) -C ./physics_noahmp/drivers/mpas driver_lib)
	($(MAKE) -C ./physics_noahmp/src src_lib)
	($(MAKE) -C ./physics_noahmp/utility utility_lib)
ifeq "$(COLM2024_ENABLED)" "true"
	ar -ru libphys.a $(COLM2024_OBJS)
endif

phys_interface: $(OBJS)

# DEPENDENCIES:
mpas_atmphys_camrad_init.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_control.o: \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver.o: \
	mpas_atmphys_driver_cloudiness.o \
	mpas_atmphys_driver_convection.o \
	mpas_atmphys_driver_gwdo.o \
	mpas_atmphys_driver_lsm.o \
	mpas_atmphys_driver_lsm_noahmp.o \
	mpas_atmphys_driver_lsm_colm2024.o \
	mpas_atmphys_driver_pbl.o \
	mpas_atmphys_driver_radiation_lw.o \
	mpas_atmphys_driver_radiation_sw.o \
	mpas_atmphys_driver_seaice.o \
	mpas_atmphys_driver_sfclayer.o \
	mpas_atmphys_driver_oml.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_interface.o \
	mpas_atmphys_sfc_diagnostics.o \
	mpas_atmphys_update.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_cloudiness.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_convection.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_gwdo.o: \
	mpas_atmphys_vars.o \
	mpas_atmphys_manager.o

mpas_atmphys_driver_lsm.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_landuse.o \
	mpas_atmphys_lsm_noahinit.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_lsm_noahmp.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_manager.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_microphysics.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_init_microphysics.o \
	mpas_atmphys_interface.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_oml.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_landuse.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_pbl.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_radiation_lw.o: \
	mpas_atmphys_camrad_init.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_driver_radiation_sw.o \
	mpas_atmphys_manager.o \
	mpas_atmphys_rrtmg_lwinit.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_radiation_sw.o: \
	mpas_atmphys_camrad_init.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_manager.o \
	mpas_atmphys_rrtmg_swinit.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_seaice.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_lsm_shared.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_sfclayer.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_finalize.o: \
	mpas_atmphys_lsm_noahmpfinalize.o \
	mpas_atmphys_driver_lsm_colm2024.o

mpas_atmphys_init.o: \
	mpas_atmphys_driver_convection.o \
	mpas_atmphys_driver_lsm.o \
	mpas_atmphys_driver_microphysics.o \
	mpas_atmphys_driver_pbl.o \
	mpas_atmphys_driver_radiation_lw.o \
	mpas_atmphys_driver_radiation_sw.o \
	mpas_atmphys_driver_sfclayer.o \
	mpas_atmphys_lsm_noahmpinit.o \
	mpas_atmphys_driver_lsm_colm2024.o \
	mpas_atmphys_landuse.o \
	mpas_atmphys_o3climatology.o \
	mpas_atmphys_vars.o

mpas_atmphys_interface.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_landuse.o: \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

mpas_atmphys_lsm_noahinit.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_lsm_noahmpinit.o: \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

mpas_atmphys_lsm_noahmpfinalize.o : \
	mpas_atmphys_vars.o

mpas_atmphys_driver_lsm_colm2024.o : mpas_atmphys_utilities.o $(COLM2024_STAMP)

mpas_atmphys_manager.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_o3climatology.o \
	mpas_atmphys_update_surface.o \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

mpas_atmphys_o3climatology.o: \
	mpas_atmphys_date_time.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_rrtmg_lwinit.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_rrtmg_swinit.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_sfc_diagnostics.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_todynamics.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_update_surface.o: \
	mpas_atmphys_date_time.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_update.o: \
	mpas_atmphys_driver_convection.o \
	mpas_atmphys_vars.o

clean:
	$(RM) *.o *.mod *.f90 libphys.a
	( cd physics_wrf; $(MAKE) clean )
	( if [ -d physics_mmm ]; then cd physics_mmm; $(MAKE) -f Makefile.mpas clean; fi; )
	( cd physics_noahmp/drivers/mpas; $(MAKE) clean )
	( cd physics_noahmp/src; $(MAKE) clean )
	( cd physics_noahmp/utility; $(MAKE) clean )
	$(RM) -r physics_colm2024/.mpas_bld
	( if [ -d physics_noaa/UGWP ]; then cd physics_noaa/UGWP; $(MAKE) clean; fi )
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(CPPINCLUDES) -I../../framework $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) $(COLM2024_INCLUDES) -I./physics_mmm -I./physics_wrf -I./physics_noahmp -I./physics_noahmp/utility -I./physics_noahmp/drivers/mpas -I./physics_noahmp/src -I./physics_noaa/UGWP -I.. -I../../framework $(MPAS_ESMF_INC)
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(HYDROSATIC) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) $(COLM2024_INCLUDES) -I./physics_mmm -I./physics_wrf -I./physics_noahmp -I./physics_noahmp/utility -I./physics_noahmp/drivers/mpas -I./physics_noahmp/src -I./physics_noaa/UGWP -I.. -I../../framework $(MPAS_ESMF_INC)
endif
