# vortex


# Modelling Von Karman Vortices in Parallel by the Lattice Boltzmann Method
## Oliver Chapman
### 31/03/2019

#### Abstract

Modelling fluids is a difficult task. It is computationally intensive, mathematically complex and expresses rich behaviours which are often hard to examine computationally. This code is designed to model Von Karman Vortices - an observable phenomenon which results in rhythmic oscilaltions from a randomly perturbed inital state.

<img src="mpi_vel.gif">

Figure 1: This gif shows the velocity in m/s of a fluid passing a circular boundary. The initial state is comprised of random fluctuations around 0.1 m/s which eventually lead to the Von Karman Vortices observed thereafter.

This code also seeks to demonstrate the Message Passing Interface (MPI) which allows tasks to be distributed amongst several computational units to be run in parallel. With the aid of the highest level of optimisation provided by gcc, this parallelised C code is capable of running 1,000 to 10,000 faster than similar python code. Parallelisation is therefore a very powerful tool for physics simulations.
