#!/bin/bash

make clean
make all

echo "============================================"
echo "Serial Nelder Mead with 20 Parameters:"
echo "============================================"
./NelderTest.out 20 
./LeeWiswall.out 20 
./DistParNelderTest.out 20 1
echo "\n>>>>> The three implementations above produce identical results!"

echo "\n\n============================================"
echo "Replace the 2 worst points on the simplex"
echo "============================================"
mpirun -np 1 ./DistParNelderTest.out 20 2
mpirun -np 2 ./LeeWiswall.out 20 
echo "\n>>>>> The two implementations above produce identical results!"

echo "\n\n============================================"
echo "Replace the 3 worst points on the simplex"
echo "============================================"
mpirun -np 1 ./DistParNelderTest.out 20 3
mpirun -np 3 ./LeeWiswall.out 20 
echo "\n>>>>> The two implementations above produce identical results!"

echo "\n\n============================================"
echo "Parallel implementations with two workers"
echo "============================================"
mpirun -np 2 ./DistParNelderTest.out 20 1
mpirun -np 2 ./LeeWiswall.out 20 
echo "\n>>>>> This takes a while."

echo "\n\n============================================"
echo "Parallel implementations with four workers"
echo "============================================"
mpirun -np 3 ./DistParNelderTest.out 20 1
mpirun -np 3 ./LeeWiswall.out 20 
echo "\n>>>>> This takes a while."