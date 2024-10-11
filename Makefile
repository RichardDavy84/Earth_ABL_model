# Include local defenitions
include make.macro

# Executable names
ABL = ABL.x
IOTEST = iotest.x

# Put objects and libraries in seperate directories
OBJ_DIR = .obj

LDFLAGS += -ldatetime -lncio -lnetcdff -lnetcdf -I $(OBJ_DIR)

# The source files for modules
SRC = mod_io.f90 mod_physics.f90
S77 = Initialization_module_NeXtSIM.for Command_module_for_NeXtSIM.for Physics.for Solver.for

# Resulting object files
OBJ = $(SRC:%.f90=$(OBJ_DIR)/%.o)
O77 = $(S77:%.for=$(OBJ_DIR)/%.o)

$(ABL): ABL.f90 $(OBJ) $(O77)
	$(FC) $(FFLAGS) $^ $(LDFLAGS) -o $@

$(OBJ): $(OBJ_DIR)/%.o : %.f90
	@mkdir -p $(OBJ_DIR)
	$(FC) -c $(F90FLAGS) $< -o $@ $(MODFLAG) $(OBJ_DIR)

$(O77): $(OBJ_DIR)/%.o : %.for
	@mkdir -p $(OBJ_DIR)
	$(FC) -c $(FFLAGS) $< -o $@ $(MODFLAG) $(OBJ_DIR)

$(IOTEST): iotest.for $(OBJ)
	$(FC) $(FFLAGS) $^ $(LDFLAGS) -o $@

test: $(IOTEST)
	./$(IOTEST)

clean:
	rm -rf $(OBJ_DIR) $(ABL) $(IOTEST)
