# Makefile for CoLM program

include include/Makeoptions
HEADER = include/define.h

INCLUDE_DIR = -Iinclude -I.bld/ -I${NETCDF_INC}
ifneq ($(strip $(PNETCDF)),)
INCLUDE_DIR += -I$(PNETCDF)/include
endif
ifneq ($(strip $(PNETCDF_INC)),)
INCLUDE_DIR += -I$(PNETCDF_INC)
endif
VPATH = include : share : mksrfdata \
	: main : main/HYDRO : main/BGC : main/URBAN : main/DA \
	: main/ParaOpt : drivers/mpas : .bld

.PHONY: all
all : mpas_objects

.NOTPARALLEL:

.PHONY: mkdir_build
mkdir_build :
	mkdir -p .bld

OBJS_SHARED =    \
				  MOD_Precision.o              \
				  MOD_SPMD_Task.o              \
				  MOD_Namelist.o               \
				  MOD_Vars_Global.o            \
				  MOD_Const_Physical.o         \
				  MOD_Const_LC.o               \
				  MOD_Utils.o                  \
				  MOD_IncompleteGamma.o        \
				  MOD_UserDefFun.o             \
				  MOD_TimeManager.o            \
				  MOD_Const_PFT.o              \
				  MOD_NetCDFSerial.o           \
				  MOD_Block.o                  \
				  MOD_Grid.o                   \
				  MOD_Pixel.o                  \
				  MOD_DataType.o               \
				  MOD_NetCDFPoint.o            \
				  MOD_NetCDFBlock.o            \
				  MOD_CatchmentDataReadin.o    \
				  MOD_5x5DataReadin.o          \
				  MOD_Mesh.o                   \
				  MOD_Pixelset.o               \
				  MOD_NetCDFVector.o           \
				  MOD_RangeCheck.o             \
				  MOD_SpatialMapping.o         \
				  MOD_WorkerPushData.o         \
				  MOD_AggregationRequestData.o \
				  MOD_PixelsetShared.o         \
				  MOD_LandElm.o                \
				  MOD_LandHRU.o                \
				  MOD_LandPatch.o              \
				  MOD_Land2mWMO.o              \
				  MOD_LandCrop.o               \
				  MOD_LandPFT.o                \
				  MOD_LandUrban.o              \
				  MOD_Urban_Const_LCZ.o        \
				  MOD_SingleSrfdata.o          \
				  MOD_SrfdataRestart.o         \
				  MOD_ElmVector.o              \
				  MOD_HRUVector.o              \
				  MOD_MeshFilter.o

${OBJS_SHARED} : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< ${MOD_CMD}.bld

OBJS_SHARED_T = $(addprefix .bld/,${OBJS_SHARED})

OBJS_BASIC =    \
				 MOD_Vector_ReadWrite.o         \
				 MOD_dataSpec_PDB.o             \
				 MOD_tav_abs.o                  \
				 MOD_prospect_DB.o              \
				 MOD_Catch_BasinNetwork.o       \
				 MOD_Catch_Vars_TimeVariables.o \
				 MOD_Catch_Vars_1DFluxes.o      \
				 MOD_Grid_RiverLakeNetwork.o    \
				 MOD_Grid_Reservoir.o           \
				 MOD_Grid_RiverLakeTimeVars.o   \
				 MOD_BGC_Vars_1DFluxes.o        \
				 MOD_BGC_Vars_1DPFTFluxes.o     \
				 MOD_BGC_Vars_PFTimeVariables.o \
				 MOD_BGC_Vars_TimeInvariants.o  \
				 MOD_BGC_Vars_TimeVariables.o   \
				 MOD_Urban_Vars_1DFluxes.o      \
				 MOD_Urban_Vars_TimeVariables.o \
				 MOD_Urban_Vars_TimeInvariants.o\
				 MOD_DA_Vars_1DFluxes.o         \
				 MOD_Vars_TimeInvariants.o      \
				 MOD_DA_Vars_TimeVariables.o    \
				 MOD_Vars_TimeVariables.o       \
				 MOD_Vars_1DPFTFluxes.o         \
				 MOD_Vars_1DFluxes.o            \
				 MOD_Vars_1DForcing.o           \
				 MOD_Hydro_SoilFunction.o       \
				 MOD_Hydro_SoilWater.o          \
				 MOD_Eroot.o                    \
				 MOD_Qsadv.o                    \
				 MOD_LAIEmpirical.o             \
				 MOD_LAIReadin.o                \
				 MOD_CropReadin.o               \
				 MOD_NitrifData.o               \
				 MOD_NdepData.o                 \
				 MOD_FireData.o                 \
				 MOD_OrbCoszen.o                \
				 MOD_OrbCosazi.o                \
				 MOD_HighRes_Parameters.o		\
				 MOD_3DCanopyRadiation.o        \
				 MOD_Aerosol.o                  \
				 MOD_SnowSnicar.o               \
				 MOD_Albedo.o                   \
				 MOD_SnowSnicar_HiRes.o			\
				 MOD_Albedo_HiRes.o				\
				 MOD_SnowFraction.o             \
				 MOD_Urban_LAIReadin.o          \
				 MOD_Urban_Shortwave.o          \
				 MOD_Urban_Albedo.o             \
				 MOD_MonthlyinSituCO2MaunaLoa.o \
				 MOD_BGC_CNSummary.o            \
				 MOD_ElementNeighbour.o         \
				 MOD_Catch_HillslopeNetwork.o   \
				 MOD_Catch_RiverLakeNetwork.o   \
				 MOD_Catch_Reservoir.o          \
				 MOD_VicParaReadin.o


$(OBJS_BASIC) : %.o : %.F90 ${HEADER} ${OBJS_SHARED}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< ${MOD_CMD}.bld

OBJS_BASIC_T = $(addprefix .bld/,${OBJS_BASIC})

OBJS_MAIN = \
				MOD_Catch_HillslopeFlow.o                 \
				MOD_Catch_SubsurfaceFlow.o                \
				MOD_Catch_RiverLakeFlow.o                 \
				MOD_Catch_Hist.o                          \
				MOD_Catch_WriteParameters.o               \
				MOD_BGC_CNCStateUpdate1.o                 \
				MOD_BGC_CNCStateUpdate2.o                 \
				MOD_BGC_CNCStateUpdate3.o                 \
				MOD_BGC_CNNStateUpdate1.o                 \
				MOD_BGC_CNNStateUpdate2.o                 \
				MOD_BGC_CNNStateUpdate3.o                 \
				MOD_BGC_Soil_BiogeochemNStateUpdate1.o    \
				MOD_BGC_Soil_BiogeochemNitrifDenitrif.o   \
				MOD_BGC_Soil_BiogeochemCompetition.o      \
				MOD_BGC_Soil_BiogeochemDecompCascadeBGC.o \
				MOD_BGC_Soil_BiogeochemDecomp.o           \
				MOD_BGC_Soil_BiogeochemLittVertTransp.o   \
				MOD_BGC_Soil_BiogeochemNLeaching.o        \
				MOD_BGC_Soil_BiogeochemPotential.o        \
				MOD_BGC_Soil_BiogeochemVerticalProfile.o  \
				MOD_BGC_Veg_CNGapMortality.o              \
				MOD_BGC_Veg_CNGResp.o                     \
				MOD_BGC_Veg_CNMResp.o                     \
				MOD_BGC_Daylength.o                       \
				MOD_BGC_Veg_CNPhenology.o                 \
				MOD_BGC_Veg_NutrientCompetition.o         \
				MOD_BGC_Veg_CNVegStructUpdate.o           \
				MOD_BGC_CNAnnualUpdate.o                  \
				MOD_BGC_CNZeroFluxes.o                    \
				MOD_BGC_CNBalanceCheck.o                  \
				MOD_BGC_CNSASU.o                          \
				MOD_BGC_Veg_CNNDynamics.o                 \
				MOD_BGC_Veg_CNFireBase.o                  \
				MOD_BGC_Veg_CNFireLi2016.o                \
				MOD_Vars_2DForcing.o                      \
				MOD_UserSpecifiedForcing.o                \
				MOD_ForcingDownscaling.o                  \
				MOD_Forcing.o                             \
				MOD_DA_TWS.o                              \
				MOD_DA_Const.o                            \
				MOD_DA_RTM.o                              \
				MOD_DA_EnKF.o                             \
				MOD_DA_SM.o                               \
				MOD_DA_Ensemble.o                         \
				MOD_DA_Main.o                             \
				MOD_Opt_Baseflow.o                        \
				MOD_ParameterOptimization.o               \
				MOD_AssimStomataConductance.o             \
				MOD_PlantHydraulic.o                      \
				MOD_FrictionVelocity.o                    \
				MOD_TurbulenceLEddy.o                     \
				MOD_Ozone.o                               \
				MOD_CanopyLayerProfile.o                  \
				MOD_LeafTemperature.o                     \
				MOD_LeafTemperaturePC.o                   \
				MOD_SoilThermalParameters.o               \
				MOD_Hydro_VIC_Variables.o                 \
				MOD_Hydro_VIC.o                           \
				MOD_Runoff.o                              \
				MOD_SoilSnowHydrology.o                   \
				MOD_SnowLayersCombineDivide.o             \
				MOD_PhaseChange.o                         \
				MOD_Glacier.o                             \
				MOD_Lake.o                                \
				MOD_SimpleOcean.o                         \
				MOD_GroundFluxes.o                        \
				MOD_GroundTemperature.o                   \
				MOD_LeafInterception.o                    \
				MOD_NetSolar.o                            \
				MOD_NetSolar_Hyper.o                      \
				MOD_WetBulb.o                             \
				MOD_RainSnowTemp.o                        \
				MOD_SoilSurfaceResistance.o               \
				MOD_NewSnow.o                             \
				MOD_Thermal.o                             \
				MOD_Vars_1DAccFluxes.o                    \
				MOD_Irrigation.o                          \
				MOD_BGC_driver.o                          \
				MOD_HistWriteBack.o                       \
				MOD_HistGridded.o                         \
				MOD_HistVector.o                          \
				MOD_HistSingle.o                          \
				MOD_Grid_RiverLakeHist.o                  \
				MOD_Hist.o                                \
				MOD_CheckEquilibrium.o                    \
				MOD_LightningData.o                       \
				MOD_Catch_LateralFlow.o                   \
				MOD_Grid_RiverLakeFlow.o                  \
				MOD_Urban_Longwave.o                      \
				MOD_Urban_NetSolar.o                      \
				MOD_Urban_Flux.o                          \
				MOD_Urban_GroundFlux.o                    \
				MOD_Urban_RoofFlux.o                      \
				MOD_Urban_RoofTemperature.o               \
				MOD_Urban_WallTemperature.o               \
				MOD_Urban_PerviousTemperature.o           \
				MOD_Urban_ImperviousTemperature.o         \
				MOD_Urban_Hydrology.o                     \
				MOD_Urban_BEM.o                           \
				MOD_Urban_LUCY.o                          \
				MOD_Urban_Thermal.o                       \
				CoLMMAIN_Urban.o                          \
					CoLMDRIVER.o                              \
					MOD_CoLM_MPAS_Interface.o                 \
					CoLMMAIN.o

$(OBJS_MAIN) : %.o : %.F90 ${HEADER} ${OBJS_SHARED} ${OBJS_BASIC}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< ${MOD_CMD}.bld

MOD_Urban_Thermal.o: MOD_Urban_Flux.o
MOD_Grid_RiverLakeTimeVars.o: MOD_Grid_RiverLakeNetwork.o MOD_Grid_Reservoir.o
MOD_Grid_RiverLakeFlow.o: MOD_Grid_RiverLakeHist.o

OBJS_MPAS_T = $(addprefix .bld/,${OBJS_SHARED} ${OBJS_BASIC} ${OBJS_MAIN})

.PHONY: mpas_objects
mpas_objects : mkdir_build ${HEADER} ${OBJS_SHARED} ${OBJS_BASIC} ${OBJS_MAIN}
	@echo ''
	@echo 'staging CoLM objects for MPAS physics library >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	rm -rf .mpas
	mkdir -p .mpas
	cd .mpas && for obj in ${OBJS_MPAS_T}; do ln -sf ../$$obj ./; done

.PHONY: clean
clean :
	rm -rf .bld
	rm -rf .mpas
