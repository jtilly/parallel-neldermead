CFLAGS = $(USERFLAGS)

all: serial parallel tests leewiswall

run: 
	./LeeWiswall.out 10 
	./NelderTest.out 10 
	./DistParNelderTest.out 10 1

serial: NelderMead.o
	g++ $(CFLAGS) -o NelderTest.out NelderMead_Driver.cpp NelderMead.o

tests: NelderMead.o
	g++ $(CFLAGS) -o tests.out tests.cpp NelderMead.o

parallel: DistParNelderMead.o
	mpicxx $(CFLAGS) -o DistParNelderTest.out DistParNelderMead_Driver.cpp DistParNelderMead.o

leewiswall: LeeWiswall.o
	mpicxx $(CFLAGS) -o LeeWiswall.out LeeWiswall_Driver.cpp LeeWiswall.o

NelderMead.o: NelderMead.cpp NelderMead.hpp
	g++ $(CFLAGS) -c NelderMead.cpp

DistParNelderMead.o: DistParNelderMead.cpp DistParNelderMead.hpp
	mpicxx $(CFLAGS) -c DistParNelderMead.cpp

LeeWiswall.o: LeeWiswall.cpp LeeWiswall.hpp
	mpicxx $(CFLAGS) -c LeeWiswall.cpp


clean:
	rm -f *.o *.out
