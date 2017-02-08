
outdir          = ./bin
app             = $(outdir)/fmm
SUPERGLUE_DIR   = /pica/h1/afshin/sg/superglue/include
SUPERGLUE_FLAGS = -pthread -I$(SUPERGLUE_DIR)  #-pedantic -Wno-long-long -Wno-format

DUCTTEIP_INC    = /pica/h1/afshin/Damavand/D4/DuctTeip/ResStudio/Prototype/include
DUCTTEIP_LIB    = /pica/h1/afshin/Damavand/D3/ductteip/bin/Debug
DT_COMP_FLAGS   = -I$(DUCTTEIP_INC) 
DT_LINK_FLAGS   = -L$(DUCTTEIP_LIB) -lductteip

ACML_DIR        = /pica/h1/afshin/acml/gnu4.9/gfortran64
ACML_LIB        = $(ACML_DIR)/lib/libacml.a
ACML_FLAGS      = -I$(ACML_DIR)/include
ACML_LINK_FLAGS = $(ACML_LIB)
SOURCE_DIR      = ./src
HEADER_DIR      = ./include
CPP             = mpic++
LINKER          = mpic++
FMM_FLAGS       = -I$(HEADER_DIR)


#-----------------------GCC Compiler set---------------------
	GCOV_FLAGS=-fprofile-arcs -ftest-coverage
	OPTIM_FLAGS=-fopenmp -mavx -march=bdver1 -mfma4 -Ofast  -Wwrite-strings
	SPECIAL_FLAGS=$(OPTIM_FLAGS) -DOMP_TASKS     
#	LINK_FLAGS=-lm -lrt -lpthread -Wl,--allow-multiple-definition -lgfortran -fopenmp $(ACML_LINK_FLAGS)
#	COMP_FLAGS= $(SUPERGLUE_FLAGS) $(ACML_FLAGS) $(FMM_FLAGS) -std=c++11 -DOMP_TASKS -g $(SPECIAL_FLAGS) 
	LINK_FLAGS=-lm -lrt -lpthread -Wl,--allow-multiple-definition -lgfortran -fopenmp $(ACML_LINK_FLAGS) $(DT_LINK_FLAGS)
	COMP_FLAGS= $(SUPERGLUE_FLAGS) $(DT_COMP_FLAGS) $(ACML_FLAGS) $(FMM_FLAGS) -std=c++11 -g $(SPECIAL_FLAGS) 


#########################################################
sources:=$(notdir $(shell ls -Sr $(SOURCE_DIR)/*.cpp))
objnames:=$(sources:%.cpp=%.o)
objects:=$(addprefix $(outdir)/,$(objnames))	
all: $(app)

$(objects):  $(outdir)/%.o:  $(SOURCE_DIR)/%.cpp
	$(info compile $(notdir $<) )
	$(CPP) -c -o $@ $< $(COMP_FLAGS)

$(app): $(objects) 
	$(LINKER) -o $(app) $(objects) $(LINK_FLAGS)

clean: 
	rm -f $(outdir)/*.o $(app)
