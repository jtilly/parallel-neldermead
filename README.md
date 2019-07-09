# Parallel Nelder-Mead
[![Build Status](https://travis-ci.org/jtilly/parallel-neldermead.svg?branch=master)](https://travis-ci.org/jtilly/parallel-neldermead)

This repository contains different implementations of the Nelder-Mead optimization algorithm. The implementations are based on
* Jeff Borggaard's [[Matlab code]](http://people.sc.fsu.edu/~jburkardt/m_src/nelder_mead/nelder_mead.html)
* Lee and Wiswall (2007) [[Paper]](http://www.econ.nyu.edu/user/wiswall/research/lee_wiswall_parallel_simplex_edit_2_8_2007.pdf)
* Klein and Neira (2014) [[Paper]](https://www.dropbox.com/s/5zw70tgueot3s4r/Klein_Neira.pdf) [[C++ Code]](  https://www.dropbox.com/s/z8gpwkszy0u5wz6/Klein_Neira_code.zip?dl=1)

Build and run the serial algorithm:
```{shell}
# build and run a simple example
make serial 
./NelderTest.out <problem_size>
# try a range of different objective functions
make tests 
./tests.out
```

Build and run the parallel Lee and Wiswall implementation:
```{shell}
# build and run a simple example
make leewiswall  
mpirun -np <num_proc> ./LeeWiswall.out <problem_size>
```

Build and run the parallel Klein and Neira implementation (this doesn't really work):
```{shell}
# build and run a simple example
make parallel 
mpirun -np <num_proc> ./DistParNelderTest.out <problem_size> <points_per_iter>
```

