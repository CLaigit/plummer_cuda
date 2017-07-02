# Plummer model simulation

This project uses plummer model to simulate star clusters. This model describes a star cluster consists of 10240 stars. 

## Theoretical part

1. The inital positions and velocities of all stars is given by Plummer 3D density profile.
2. Leapfrog integration is used to integrate this hamiltonian system. 

## Programming part

This project is written in C. Leapfrog integration part is sped up with CUDA.

## Simulation result

1. stable star clusters. If the velocity distribution is set by Plummer model, the star clusters will maintain in a stable status.
2. unstable star clusters. If the velocity distribution is set randomly, the star clusters will spread out after a while.
