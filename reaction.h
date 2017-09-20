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

    double PeriodicDistance(const Particle& b, const double xLength, const double yLength) {
        const double xDiff = b.x - x;
        const double yDiff = b.y - y;

        std::vector<double> results(9);

        // Center
        results[0] = std::sqrt(std::pow(xDiff, 2.0) + std::pow(yDiff, 2.0));

        // Bottom Left
        results[1] = std::sqrt(std::pow(xDiff-xLength, 2.0) + std::pow(yDiff-yLength, 2.0));

        // Bottom
        results[2] = std::sqrt(std::pow(xDiff, 2.0) + std::pow(yDiff-yLength, 2.0));

        // Bottom Right
        results[3] = std::sqrt(std::pow(xDiff+xLength, 2.0) + std::pow(yDiff-yLength, 2.0));

        // Left
        results[4] = std::sqrt(std::pow(xDiff-xLength, 2.0) + std::pow(yDiff, 2.0));

        // Right
        results[5] = std::sqrt(std::pow(xDiff+xLength, 2.0) + std::pow(yDiff, 2.0));

        // Top Left
        results[6] = std::sqrt(std::pow(xDiff-xLength, 2.0) + std::pow(yDiff+yLength, 2.0));

        // Top
        results[7] = std::sqrt(std::pow(xDiff, 2.0) + std::pow(yDiff+yLength, 2.0));

        // Top Right
        results[8] = std::sqrt(std::pow(xDiff+xLength, 2.0) + std::pow(yDiff+yLength, 2.0));

        return *std::min_element( results.begin(), results.end() );
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
void UpdateParticles(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double TimeStep, const double Diffusion, const unsigned int FieldWidth, const unsigned int FieldHeight, std::uniform_real_distribution<double> &mRandom, std::mt19937_64 &gen );

#endif // BUILD_REACTION_H_
