# Set the fortran compiler to nvfortran
FC = nvfortran

# Set device compute capability and cuda version
CCXY = cc86
CUDAXY = cuda11.7

# Set the flags for the nnvfortran compiler
FFLAGS = -fast -Minfo=all -gpu=$(CCXY),$(CUDAXY),managed,lineinfo -cuda -acc -mp=multicore

# Set the source files
SRC = fem1d.f90

# Set the list of objects
OBJ = $(SRC:.f90=.o)

# Set the name of the executable
EXE = FEM1D.x

# Set the all target
all: $(EXE)

# Compile the object
$(OBJ): $(SRC)
	$(FC) $(FFLAGS) -c $<

# Link the executable
$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^
	rm -rf *.o *.mod

# Clean the directory
clean:
	rm -f *.o *.mod *.x
