# Set the fortran compiler to nvfortran
FC = nvfortran
#FC = gfortran

# Set device compute capability and cuda version
CCXY = cc61
CUDAXY = cuda11.8

# Options
USE_OMP = 1
USE_ACC = 1

# Set the flags for the nnvfortran compiler
ifeq ($(FC),nvfortran)
	FFLAGS = -cpp -fast -Minfo=all -gpu=$(CCXY),$(CUDAXY),managed,lineinfo
	ifeq ($(USE_OMP),1)
		FFLAGS += -DUSE_OMP -mp=gpu
	endif
	ifeq ($(USE_ACC),1)
		FFLAGS += -DUSE_ACC -acc
	endif
else ifeq ($(FC),gfortran)
	FFLAGS = -cpp -O3 -march=native
	ifeq ($(USE_OMP),1)
		FFLAGS += -DUSE_OMP -fopenmp
	endif
endif

# Set the source files
SRC = conjGrad.f90

# Set the list of objects
OBJ = $(SRC:.f90=.o)

# Set the name of the executable
EXE = conjGrad.x

# Set the all target
all: $(EXE)

# Compile the object
$(OBJ): $(SRC)
	$(FC) $(FFLAGS) -c $<

# Link the executable
$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^
	rm -rf *.o *.mod

# Rule for running the executable with nsys (set OMP_NUM_THREADS to a number between 1 and 80 before running)
run:
	nsys profile --stats=true --force-overwrite=true --trace=cuda,nvtx,osrt,openacc,openmp --cuda-memory-usage=true --output=conjGrad_report ./$(EXE)

# Clean the directory
clean:
	rm -f *.o *.mod *.x *.nsys-rep *.sqlite *.qdrep
