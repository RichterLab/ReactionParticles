# Reaction Particles
Code to simulate particle's movement through random walk and annihilation through proximity and random chance with other particles.

## Installation

```
$ module load cmake/3.6.3
$ module load gcc/5.2.0
$ module load cuda/8.0
$ git clone https://github.com/RichterLab/ReactionParticles.git
$ cd ReactionParticles
$ mkdir build
$ cd build
$ cmake -DBUILD_CUDA=ON ..
$ make
```

## Branches
To switch between GPU and CPU versions of the code, you can select the branch with ```git checkout``` before the ```cmake ..``` step.
The GPU branch is ```gpu_cuda``` and the CPU branch is ```master```

## Usage

### Generate Input Matrices
To generate an input velocity, you must execute generate_fields within matlab.
The resulting matrix must be saved as an ascii file for the software to read it. 
The code, velocity_function.m, already saves the file in this format.

### Execute Simulation
```
$ ./reaction-particle -v vel1.mat -p 5000
```

#### Execution Parameters
The executable provides a help option to list command line arguments. 
```
$ ./reaction-particle --help
```

Most of these parameters have a relation to the generation of the velocity fields:
* width - This must match the value of Lx in generate_fields.m
* height - This must match the value of Ly in generate_fields.m
* grid - This must match the value specified for the multiplier for numgridy and numgridx in permeability.m
