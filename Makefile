
outdir= ./bin
app= $(outdir)/fmm
SUPERGLUE_DIR=/pica/h1/afshin/sg/superglue/include
SUPERGLUE_FLAGS=-pthread -I$(SUPERGLUE_DIR)  #-pedantic -Wno-long-long -Wno-format
SOURCE_DIR=./src
HEADER_DIR=./include
CPP=g++



#-----------------------GCC Compiler set---------------------
	GCOV_FLAGS=-fprofile-arcs -ftest-coverage
	OPTIM_FLAGS=-mavx -march=bdver1 -mfma4 -Ofast  -Wwrite-strings
	SPECIAL_FLAGS=$(OPTIM_FLAGS)     
	LINK_FLAGS=-lm -lrt -lpthread
	COMP_FLAGS= $(SUPERGLUE_FLAGS) $(DUCTTEIP_FLAGS) -std=c++11 $(SPECIAL_FLAGS) 

#########################################################
headers:=$(notdir $(shell ls -Sr $(SOURCE_DIR)/*.cpp))
objnames:=$(headers:%.cpp=%.o)
objects:=$(addprefix $(outdir)/,$(objnames))	
all: $(app)


$(objects):  $(outdir)/%.o:  $(SOURCE_DIR)/%.cpp #$(HEADER_DIR)/%.hpp
	$(info compile $(notdir $<) )
	@$(CPP) -c -o $@ $< $(COMP_FLAGS)
$(app): $(objects)
	@$(CPP) -o $(app) $(objects) 


clean: 
	rm -f $(outdir)/*.o $(app)
