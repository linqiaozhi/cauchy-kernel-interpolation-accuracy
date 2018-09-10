PROJECT = int2

#HOST = macosx
#HOST = linux-gfortran
HOST = linux-gfortran-openmp

ifeq ($(HOST),macosx)
FC = gfortran -c -w
FFLAGS = -O3
FLINK = gfortran -w -o $(PROJECT) -framework accelerate
endif

ifeq ($(HOST),linux-gfortran)
FC = gfortran -c -w
FFLAGS = -O3
FLINK = gfortran -w -o $(PROJECT) 
endif

ifeq ($(HOST),linux-gfortran-openmp)
FC = gfortran 
FFLAGS =  -c -w -O3 --openmp
FLINK = gfortran -w -o $(PROJECT) --openmp
endif




.PHONY: all clean list


SOURCES =  fitsnemt.f \
  durdec.f \
  fourt2d.f \
  prini.f \
  hkrand.f \
  dlaran.f \


OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))
#

%.o : %.f
	$(FC) $(FFLAGS) $< -o $@


all: $(OBJECTS)
	rm -f $(PROJECT)
	$(FLINK) $(OBJECTS) $(FEND)
	./$(PROJECT)

clean:
	rm -f $(OBJECTS)
	rm -f $(PROJECT)

list: $(SOURCES)
	$(warning Requires:  $^)





  
