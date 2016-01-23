# Parallel Nelder Mead

This repository contains three implementations of the Nelder-Mead algorithm. The implementations are based on
* Jeff Borggaard's [Matlab code](http://people.sc.fsu.edu/~jburkardt/m_src/nelder_mead/nelder_mead.html)
* Lee and Wiswall (2007) [Paper](http://www.econ.nyu.edu/user/wiswall/research/lee_wiswall_parallel_simplex_edit_2_8_2007.pdf)
* Klein and Neira (2014) [Paper](http://www.cs.ucsb.edu/~kyleklein/publications/neldermead.pdf) [C++ Code](https://dl.dropboxusercontent.com/u/17629709/Klein_Neira_code.zip)

Build and run the serial algorithm:
```{shell}
make serial
./NelderTest <problem_size>
```

Build and run the Klein and Neira implementation:
```{shell}
make parallel
runmpi -np <num_proc> ./NelderTest <problem_size> <points_per_iter>
```

Build and run the Lee and Wiswall implementation:
```{shell}
make leewiswall
runmpi -np <num_proc> ./NelderTest <problem_size>
```

