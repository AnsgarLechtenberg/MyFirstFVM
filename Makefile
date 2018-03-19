FC=gfortran
#FC=ifort
CFLAGS +=-Wall
CFLAGS += -g
CFLAGS += -ffpe-trap=invalid,zero,overflow
#CFLAGS += -O3

SOURCE_DIR = src
OBJECTS_DIR = obj
EXECUTABLE_DIR = bin

VPATH := :$(SOURCE_DIR)
LIB = lib/liblapack.a lib/librefblas.a

_EXE = fv_explizit \
       fv_implizit \
       fv_preCon_explizit \
       fv_preCon_implizit \
       test

EXE = $(patsubst %,$(EXECUTABLE_DIR)/%,$(_EXE))

all: FORCE $(EXE)

FORCE:
	@mkdir -p $(EXECUTABLE_DIR)
	@mkdir -p $(OBJECTS_DIR)


.PHONY: FORCE clean \
	all 

.SILENT: clean parallel FORCE

clean:
	rm -rvf $(EXECUTABLE_DIR)
	rm -rvf $(OBJECTS_DIR)
	rm -vrf *.mod *.o
	rm -vrf *~
	rm -vrf $(SOURCE_DIR)/*.mod
	rm -vrf fort.*

$(OBJECTS_DIR)/modules.o: modules.F90 $(LIB)
	$(FC) $(CFLAGS) $< -J $(SOURCE_DIR) -c -o $@

$(EXECUTABLE_DIR)/% : %.F90 $(OBJECTS_DIR)/modules.o $(LIB)
	$(FC) $(CFLAGS) $^ -I $(SOURCE_DIR) -o $@

