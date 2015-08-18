# for gfortran
#
BINDIR="../"

FLAGSALWAYS =  #-ffree-form
FLAGSOPT= -O3 #-DMPI
#FLAGSDEBUG= -g

#F90=gfortran -DGFORTRAN
F90=mpif90 -DMPI
#LIBS=-L/usr/local/openmpi-1.6.3-mx-intel12/lib/
 
# choose debugging or optimization
FLAGS= ${FLAGSMACHINE} ${FLAGSALWAYS} ${LIBS} ${FLAGSOPT}  #  change the last to ${FLAGSOPT} or ${FLAGSDEBUG}


MY_DIR=`basename ${PWD}`
#------------------------------------

%.o: %.f90
	$(F90) $(FLAGS) -c $< -o $@

%.o: %.F90
	$(F90) $(FLAGS) -c $< -o $@

all:	2pcf 3pcf add_fiber

2pcf: 2pcf.o kdtree2.o extension.o
	${F90} ${FLAGS} -o 2pcf 2pcf.o kdtree2.o extension.o ; \
	mv $@ $(BINDIR)

2pcf.o: kdtree2.o extension.o

3pcf: 3pcf.o kdtree2.o extension.o
	${F90} ${FLAGS} -o 3pcf 3pcf.o kdtree2.o extension.o ; \
	mv $@ $(BINDIR)

3pcf.o: kdtree2.o

add_fiber: add_fiber.o kdtree2.o extension.o
	${F90} ${FLAGS} -o add_fiber add_fiber.o kdtree2.o extension.o ; \
	mv $@ $(BINDIR)

add_fiber.o: kdtree2.o

clean:
	rm -f *.o *.mod
