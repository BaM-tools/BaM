#========================================================================
# Makefile to build BaM for LINUX/WINDOWS/MACOS based on MiniDMSL
#========================================================================
#

# Define the Fortran Compiler and its options
FC  = gfortran
FLAGS = -fcheck=all -Wall -O3 -static # DEBUG: -fcheck=all -Wall -g -O0

# Define directories
DMSL_DIR = ../../miniDMSL/src/
BMSL_DIR = ../../BMSL/src/
Models_DIR = ../src/Models/
BaM_DIR = ../src/

# Define DMSL files
DMSL = \
	$(DMSL_DIR)kinds_dmsl_kit.f90              \
	$(DMSL_DIR)types_dmsl_kit.f90              \
	$(DMSL_DIR)utilities_dmsl_kit.f90          \
	$(DMSL_DIR)mt19937-64.f95                  \
	$(DMSL_DIR)uniran1_minidmsl_mod.f90        \
	$(DMSL_DIR)linalg_dmsl_kit.f90             \
	$(DMSL_DIR)numerix_dmsl_kit.f90            \
	$(DMSL_DIR)dynamicModelLibrary.f90

# Define BMSL files
BMSL = \
	$(BMSL_DIR)Geodesy_tools.f90                 \
	$(BMSL_DIR)Dates_tools.f90                   \
	$(BMSL_DIR)Aggregation_tools.f90             \
	$(BMSL_DIR)EmpiricalStats_tools.f90          \
	$(BMSL_DIR)TimeSeries_tools.f90              \
	$(BMSL_DIR)DataRW_tools.f90                  \
	$(BMSL_DIR)Distribution_tools.f90            \
	$(BMSL_DIR)MCMC_tools.f90                    \
	$(BMSL_DIR)MCMCStrategy_tools.f90            \
	$(BMSL_DIR)BayesianEstimation_tools.f90

# Define Models files
Models = \
	$(Models_DIR)Linear_model.f90                      \
	$(Models_DIR)Mixture_model.f90                     \
	$(Models_DIR)BaRatin_model.f90                     \
	$(Models_DIR)BaRatinBAC_model.f90                  \
	$(Models_DIR)StageFallDischarge_model.f90          \
	$(Models_DIR)SFDTidal_model.f90                    \
	$(Models_DIR)SFDTidal2_model.f90                   \
	$(Models_DIR)SFDTidal4_model.f90                   \
	$(Models_DIR)SFDTidalJones_model.f90               \
	$(Models_DIR)StageGradientDischarge_model.f90      \
	$(Models_DIR)Vegetation_model.f90                  \
	$(Models_DIR)SWOT_model.f90                        \
	$(Models_DIR)SedimentTransport_model.f90           \
	$(Models_DIR)SuspendedLoad_model.f90               \
	$(Models_DIR)Orthorectification_model.f90          \
	$(Models_DIR)GR4J_model.f90                        \
	$(Models_DIR)Recession_model.f90                   \
	$(Models_DIR)Segmentation_model.f90                \
	$(Models_DIR)TextFile_model.f90                    \
	$(Models_DIR)AlgaeBiomass_model.f90                \
	$(Models_DIR)DynamicVegetation_model.f90           \
	$(Models_DIR)TidalODE_model.f90                    \
	$(Models_DIR)TidalRemenieras_model.f90             \
	$(Models_DIR)SFDTidal_Sw_correction_model.f90      \
	$(Models_DIR)SFDTidal_Qmec0_model.f90              \
	$(Models_DIR)SFDTidal_Qmec_model.f90               \
	$(Models_DIR)SFDTidal_Qmec2_model.f90              \
	$(Models_DIR)SMASH/SMASH_model.f90                 \
	$(Models_DIR)MAGE_model.f90	                       \
	$(Models_DIR)MAGE_TEMP_model.f90	               \
	$(Models_DIR)MAGE_ZQV_model.f90	               \
	$(Models_DIR)HydraulicControl_section_model.f90

# Define BaM files
BaM = \
	$(BaM_DIR)ModelLibrary_tools.f90       \
	$(BaM_DIR)BaM_tools.f90                \
	$(BaM_DIR)BaM_main.f90

all:build_exe

build_exe: compile_bam
	$(FC)  $(FLAGS) $(DMSL) $(BMSL) $(Models) $(BaM) -o BaM

compile_bam: compile_models
	$(FC) $(FLAGS) -c $(BaM)

compile_models: compile_bmsl
	$(FC) $(FLAGS) -c $(Models)

compile_bmsl: compile_dmsl
	$(FC) $(FLAGS) -c $(BMSL)

compile_dmsl:
	$(FC) $(FLAGS) -c $(DMSL)





