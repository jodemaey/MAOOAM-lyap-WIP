
ifndef COMPILER
  COMPILER=ifort
endif
ifdef INTMODE
  # Define integrator module
  INTEGRATOR = $(INTMODE)_integrator
  TL_AD_INTEGRATOR = $(INTMODE)_tl_ad_integrator
else
  # Define integrator module
  INTEGRATOR = rk2_integrator
  TL_AD_INTEGRATOR = rk2_tl_ad_integrator
endif


# Define relevant compiler options for debug case and normal case for each compiler
# Normal Case:
IFORT_COMPILER_FLAGS = -assume byterecl -O3
IFORT_LD_FLAGS = -mkl
GCC_COMPILER_FLAGS = -frecord-marker=8 -O3 -funroll-loops
GCC_LD_FLAGS = -llapack -lblas
# Debug Case:
IFORT_COMPILER_FLAGS_DBG = -O0 -check all -traceback -fpe0 -check bounds -debug all -check uninit -ftrapuv -assume byterecl
IFORT_LD_FLAGS_DBG = -llapack -lblas
GCC_COMPILER_FLAGS_DBG = -frecord-marker=8 -O3 
GCC_LD_FLAGS_DBG = -llapack -lblas


#Adapt to host system
ifeq (,$(findstring Windows,$(OS)))
  RM = rm -f
else
  RM = del
endif

#Derive gfortran version from gcc version and check if >=4.6.0.

ifeq ($(COMPILER),gfortran)
  GVERSION := $(shell gcc -dumpversion)
  GE40600 := $(shell expr `gcc -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40600)
  
  ifeq ($(COMPILER),gfortran)
    ifeq ($(GE40600), 0)
      $(error GCC/gfortran version $(GVERSION) detected. This program uses features of gfortran >= 4.6.0. Please update your toolchain)
    endif
  else
    $(warning This code was optimized for gfortran >= 4.6.0 and has not been tested with other compilers)
  endif
  COMPILER_FLAGS = $(GCC_COMPILER_FLAGS)
  LD_FLAGS = $(GCC_LD_FLAGS)
  COMPILER_FLAGS_DBG = $(GCC_COMPILER_FLAGS_DBG)
  LD_FLAGS_DBG = $(GCC_LD_FLAGS_DBG)
  IFPORT = ifport.mod

endif
ifeq ($(COMPILER),ifort)
  COMPILER_FLAGS = $(IFORT_COMPILER_FLAGS)
  LD_FLAGS = $(IFORT_LD_FLAGS)
  COMPILER_FLAGS_DBG = $(IFORT_COMPILER_FLAGS_DBG)
  LD_FLAGS_DBG = $(IFORT_LD_FLAGS_DBG)
endif


.SUFFIXES:
.SUFFIXES: .f90 .o .mod .out .test

MODULES = util.mod params.mod inprod_analytic.mod aotensor_def.mod ic_def.mod ic_def_ext.mod $(INTEGRATOR).mod stat.mod tensor.mod tl_ad_tensor.mod $(TL_AD_INTEGRATOR).mod lyap_vectors.mod lyap_stat.mod $(IFPORT) energy_stat.mod
MODULE_OBJECTS = $(MODULES:.mod=.o)

PROGRAMS = maooam maooam_lyapvectors maooam_ene

all: $(PROGRAMS)
basic: maooam 
lyap: maooam_lyapvectors 
energy: maooam_ene
debug: COMPILER_FLAGS = $(COMPILER_FLAGS_DBG) 
debug: LD_FLAGS = $(LD_FLAGS_DBG)
debug: all

%.mod: %.f90  
	$(COMPILER) $(LD_FLAGS) $(COMPILER_FLAGS) -c $<

%.o: %.f90 $(MODULES)
	$(COMPILER) $(LD_FLAGS) $(COMPILER_FLAGS) -c $<

%: %.o $(MODULE_OBJECTS)
	$(COMPILER) $(COMPILER_FLAGS) $^ -o $@ $(LD_FLAGS)

inprod_analytic.mod: params.mod
aotensor_def.mod: tensor.mod inprod_analytic.mod
$(INTEGRATOR).mod: aotensor_def.mod
ic_def.mod: util.mod aotensor_def.mod
ic_def_ext.mod: util.mod aotensor_def.mod
stat.mod: params.mod
tl_ad_tensor.mod: aotensor_def.mod
$(TL_AD_INTEGRATOR).mod: tl_ad_tensor.mod
lyap_vectors.mod: params.mod
util.mod: $(IFPORT)
energy_stat.mod: params.mod inprod_analytic.mod aotensor_def.mod

clean:
	$(RM) *.o *.mod $(PROGRAMS) $(TEST_PROGRAMS) tests/*.out

.PHONY: clean all test %.test
