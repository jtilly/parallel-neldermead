The files contain two versions of the NelderMead Simplex method.

The first is a standard serial implementation contained in NelderMead.cpp and 
NelderMead.hpp. An example test program is contained in NelderMead_Driver.cpp.

The test program may be built using `make serial`. The correct usage is as 
follows.

./NelderTest <problem size>

For example, to run a problem size of 5 points, the command would be as follows.

./NelderTest 5

The second is our parallel distributed memory implementation contained in 
DistParNelderMead.cpp and DistParNelderMead.hpp. An example test program is 
contained in DistParNelderMead_Driver.cpp.

The test program may be built using `make parallel`, though note that the 
required MPI libraries must be installed. The correct usage is as 
follows.

./DistParNelderTest <problem size> <points per iter> <max iterations (optional)>

For example, to run on the SDSC Triton using two processors, on a problem size 
of 8, with two points updated on each processor per iteration, the command would
be as follows.

mpirun -machinefile $PBS_NODEFILE -np 2 ./DistParNelderTest 8 2

For additional customizations please consult the header files, as we have done
our best to comment each method. If something is not clear, feel free to contact
kyleklein@cs.ucsb.edu or neira@econ.ucsb.edu.