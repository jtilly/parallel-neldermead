# Parallel Nelder-Mead

This repository contains three implementations of the Nelder-Mead optimization algorithm. The implementations are based on
* Jeff Borggaard's [[Matlab code]](http://people.sc.fsu.edu/~jburkardt/m_src/nelder_mead/nelder_mead.html)
* Lee and Wiswall (2007) [[Paper]](http://www.econ.nyu.edu/user/wiswall/research/lee_wiswall_parallel_simplex_edit_2_8_2007.pdf)
* Klein and Neira (2014) [[Paper]](http://www.cs.ucsb.edu/~kyleklein/publications/neldermead.pdf) [[C++ Code]](https://dl.dropboxusercontent.com/u/17629709/Klein_Neira_code.zip)

Build and run the serial algorithm:
```{shell}
# build and run a simple example
make serial 
./NelderTest.out <problem_size>
# try a range of different objective functions
make tests 
./tests.out
```

Build and run the Lee and Wiswall implementation:
```{shell}
# build and run a simple example
make leewiswall  
mpirun -np <num_proc> ./LeeWiswall.out <problem_size>
```

Build and run the Klein and Neira implementation:
```{shell}
# build and run a simple example
make parallel 
mpirun -np <num_proc> ./DistParNelderTest.out <problem_size> <points_per_iter>
```

