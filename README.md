# Parallel Traveling Salesman Problem

This is an implementation of a parallel algorithm for the Traveling Salesman Problem (TSP). The TSP is a classic optimization problem in which we are given a set of cities, and we want to find the shortest tour that visits each city exactly once and returns to the starting city.

The parallel algorithm used in this implementation is based on the Simulated Annealing algorithm. Simulated Annealing is a stochastic optimization algorithm that is particularly well-suited to problems with a large search space and many local optima. The algorithm works by iteratively modifying a candidate solution by making small random changes, and accepting the new solution if it improves the objective function (i.e., reduces the length of the tour). The acceptance probability is determined by a temperature parameter that gradually decreases over time, allowing the algorithm to escape from local optima.

In the parallel version of the algorithm, multiple instances of the algorithm are run in parallel, each starting from a different initial solution. The solutions are periodically exchanged between the instances to allow for exploration of different regions of the search space. This helps to avoid getting stuck in local optima.

## Implementation

The implementation is in C++ using the Message Passing Interface (MPI) library for parallelization. The program reads in a list of cities from a file, and the TSP is solved using the Simulated Annealing algorithm with parallel tempering. The algorithm is run for a specified number of iterations, and the best solution found across all instances is output at the end.

## Usage

To use the program, you need to have MPI installed on your system. You can compile the program using a C++ compiler and the MPI library, as follows:

```
mpic++ tsp.cpp -o tsp
```

You can then run the program using the `mpiexec` command, specifying the number of instances to run in parallel and the input file containing the cities:

```
mpiexec -n <num_instances> ./tsp <input_file>
```

For example, to run the program with 4 instances and an input file called `cities.txt`, you would use:

```
mpiexec -n 4 ./tsp cities.txt
```

The program outputs the best solution found and its length.

## Input Format

The input file should contain one line per city, with each line specifying the x and y coordinates of the city. The coordinates should be separated by a space. For example:

```
0.779073 0.287351
0.064111 0.570261
0.90577 0.655971
0.511472 0.238187
0.575265 0.42304
0.794072 0.846565
0.862257 0.343014
0.21567 0.629088
0.22838 0.0832957
0.230298 0.276451
```

## Future Work

There are several ways in which this implementation could be improved or extended. For example:

- Implement other optimization algorithms for the TSP, such as genetic algorithms or ant colony optimization.
- Add more advanced parallelization techniques, such as hybrid parallelization using both MPI and OpenMP.
- Explore different parallelization strategies, such as parallel divide-and-conquer or parallel branch-and-bound.
- Develop better methods for exchanging solutions between instances, such as using more advanced communication patterns or adapting the exchange frequency based on performance metrics.
- Optimize the code for performance, for example by using more efficient data structures or implementing specialized functions for common operations.

## References

1. [Student Projects in Parallel Computing](https://www.ccs.neu.edu/home/gene/projects.html#simulated_annealing)
2. [Parallelizing simulated annealing algorithms based on high-performance computer](https://www.cs.uoi.gr/~lagaris/GRAD_GLOPT/projects/parallel_SA.pdf)
3. [Implenting a Parallel Simulated Annealing
Algorithm](https://sci-hub.hkvisa.net/10.1007/978-3-642-14390-8_16)
4. [Fast parallel simulated annealing for
traveling salesman problem on SIMD
machines with linear interconnections ](https://sci-hub.hkvisa.net/10.1016/s0167-8191(05)80107-3)

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.