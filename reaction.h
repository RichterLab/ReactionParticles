#ifndef BUILD_REACTION_H_
#define BUILD_REACTION_H_

#include <random>
#include <iostream>
#include <vector>
#include <cmath>

#ifndef BUILD_CUDA
#define DEVICE
#define HOST
#define GLOBAL
#define SHARED
#define CONSTANT
#else
#include <curand.h>
#include <curand_kernel.h>

#define CUDA_BLOCK_THREADS 128
#define DEVICE __device__
#define HOST __host__
#define GLOBAL __global__
#define CONSTANT __constant__
#define SHARED __shared__
#define gpuErrchk(ans) \
	{ gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
	if(code != cudaSuccess) {
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if(abort) exit(code);
	}
}
#endif

struct Index {
    int x, y;
    HOST DEVICE Index(int x, int y) : x(x), y(y) {}
};

struct Particle{
    bool Alive;
    double x, y, u, v;

    HOST DEVICE float PeriodicDistance(const Particle& b, const double xLength, const double yLength) {
        const float xDiff = b.x - x;
        const float yDiff = b.y - y;

        // Center
        float result = std::sqrtf(std::powf(xDiff, 2.0) + std::powf(yDiff, 2.0));

        // Bottom Left
        result = min(result, std::sqrtf(std::powf(xDiff-xLength, 2.0) + std::powf(yDiff-yLength, 2.0)));

        // Bottom
        result = min(result, std::sqrtf(std::powf(xDiff, 2.0) + std::powf(yDiff-yLength, 2.0)));

        // Bottom Right
        result = min(result, std::sqrtf(std::powf(xDiff+xLength, 2.0) + std::powf(yDiff-yLength, 2.0)));

        // Left
        result = min(result, std::sqrtf(std::powf(xDiff-xLength, 2.0) + std::powf(yDiff, 2.0)));

        // Right
        result = min(result, std::sqrtf(std::powf(xDiff+xLength, 2.0) + std::powf(yDiff, 2.0)));

        // Top Left
        result = min(result, std::sqrtf(std::powf(xDiff-xLength, 2.0) + std::powf(yDiff+yLength, 2.0)));

        // Top
        result = min(result, std::sqrtf(std::powf(xDiff, 2.0) + std::powf(yDiff+yLength, 2.0)));

        // Top Right
        result = min(result, std::sqrtf(std::powf(xDiff+xLength, 2.0) + std::powf(yDiff+yLength, 2.0)));

        return result;
    }
};

template<typename T>
HOST DEVICE const T LinearAccess(T* array, const size_t x, const size_t y, const size_t Width) {
    return array[x + y * Width];
}

template<typename T>
HOST DEVICE void LinearSet(T* array, const size_t x, const size_t y, const size_t Width, const T value) {
    array[x + y * Width] = value;
}

GLOBAL void UpdateConcentration(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double dx, const double dy, const size_t ConcentrationWidth, const size_t ConcentrationHeight, unsigned int *CountA, unsigned int *CountB);
GLOBAL void Interpolate(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const size_t VelocityWidth, const size_t VelocityHeight, const double VelocityDX, const double VelocityDY, double *U, double *V);

#ifdef BUILD_CUDA
GLOBAL void InitializeRandom(unsigned int seed, curandState_t* states);
#endif

#ifdef BUILD_CUDA
GLOBAL void UpdateParticles(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double TimeStep, const double Diffusion, const unsigned int FieldWidth, const unsigned int FieldHeight, curandState_t* states );
GLOBAL void UpdateReactions(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double TimeStep, const double Diffusion, const double ReactionProbability, const unsigned int FieldWidth, const unsigned int FieldHeight, curandState_t* states);
#else
void UpdateParticles(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double TimeStep, const double Diffusion, const unsigned int FieldWidth, const unsigned int FieldHeight, std::uniform_real_distribution<double> &mRandom, std::mt19937_64 &gen );
void UpdateReactions(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double TimeStep, const double Diffusion, const double ReactionProbability, const unsigned int FieldWidth, const unsigned int FieldHeight, std::uniform_real_distribution<double> &mRandom, std::mt19937_64 &gen);
#endif

#endif // BUILD_REACTION_H_
