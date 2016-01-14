CFLAGS = $(USERFLAGS)

serial: NelderMead.o
	g++ $(CFLAGS) -o NelderTest.out NelderMead_Driver.cpp NelderMead.o

parallel: DistParNelderMead.o
	mpicxx $(CFLAGS) -o DistParNelderTest.out DistParNelderMead_Driver.cpp DistParNelderMead.o

NelderMead.o : NelderMead.cpp NelderMead.hpp
	g++ $(CFLAGS) -c NelderMead.cpp

DistParNelderMead.o : DistParNelderMead.cpp DistParNelderMead.hpp
	mpicxx $(CFLAGS) -c DistParNelderMead.cpp


clean:
	rm -f *.o *.out
