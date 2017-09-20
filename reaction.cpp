#include "reaction.h"

#include <iostream>

GLOBAL void UpdateConcentration(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double dx, const double dy, const unsigned int ConcentrationHeight, unsigned int *CountA, unsigned int *CountB, size_t CountWidth){
    int index_start = 0, index_stride = 1;
    #ifdef BUILD_CUDA
        index_start = blockIdx.x * blockDim.x + threadIdx.x;
        index_stride = blockDim.x * gridDim.x;
    #endif

    for(int particle = 0; particle < Particles; particle += index_stride) {
        const unsigned int ax = std::floor(mParticleA[particle].x / dx);
        const unsigned int ay = (ConcentrationHeight-1) - std::floor(mParticleA[particle].y / dy);
        LinearSet(CountA, ax, ay, CountWidth, LinearAccess(CountA, ax, ay, CountWidth) + 1);

        const unsigned int bx = std::floor(mParticleB[particle].x / dx);
        const unsigned int by = (ConcentrationHeight-1) - std::floor(mParticleB[particle].y / dy);
        LinearSet(CountB, bx, by, CountWidth, LinearAccess(CountB, bx, by, CountWidth) + 1);
    }
}

GLOBAL void Interpolate(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const size_t VelocityWidth, const size_t VelocityHeight, const double VelocityDX, const double VelocityDY, double *U, double *V) {
    int index_start = 0, index_stride = 1;
    #ifdef BUILD_CUDA
        index_start = blockIdx.x * blockDim.x + threadIdx.x;
        index_stride = blockDim.x * gridDim.x;
    #endif

    for(int particle = 0; particle < Particles; particle += index_stride) {
        // Interpolate Particle A
        Particle& A = mParticleA[particle];
        const Index posA( std::floor(A.x / VelocityDX), (VelocityHeight-1) - std::floor(A.y / VelocityDY));
        Index posA2(posA.x+1, posA.y-1);

        if (posA2.x == VelocityWidth) posA2.x = 0;
        if (posA2.y == -1) posA2.y = VelocityHeight-1;

        const double aa = posA2.x * VelocityDX;
        const double ab = (VelocityHeight - posA.y - 1) * VelocityDY;
        const double ac = (VelocityHeight - posA2.y - 1) * VelocityDY;
        const double ad = posA.x * VelocityDX;

        A.u = ((aa-A.x)*(A.y-ab)*LinearAccess(U, posA2.y, posA.x, VelocityHeight)+(aa-A.x)*(ac-A.y)*LinearAccess(U, posA.y, posA.x, VelocityHeight)+(A.x-ad)*(ac-A.y)*LinearAccess(U, posA.y, posA2.x, VelocityHeight)+(A.x-ad)*(A.y-ab)*LinearAccess(U, posA2.y, posA2.x, VelocityHeight))/((aa-ad)*(ac-ab));
        A.v = ((aa-A.x)*(A.y-ab)*LinearAccess(V, posA2.y, posA.x, VelocityHeight)+(aa-A.x)*(ac-A.y)*LinearAccess(V, posA.y, posA.x, VelocityHeight)+(A.x-ad)*(ac-A.y)*LinearAccess(V, posA.y, posA2.x, VelocityHeight)+(A.x-ad)*(A.y-ab)*LinearAccess(V, posA2.y, posA2.x, VelocityHeight))/((aa-ad)*(ac-ab));

        // Interpolate Particle B
        Particle& B = mParticleB[particle];
        const Index posB( std::floor(B.x / VelocityDX), (VelocityHeight-1) - std::floor(B.y / VelocityDY));
        Index posB2(posB.x+1, posB.y-1);

        if (posB2.x == VelocityWidth) posB2.x = 0;
        if (posB2.y == -1) posB2.y = VelocityHeight-1;

        const double ba = posB2.x * VelocityDX;
        const double bb = (VelocityHeight - posB.y - 1) * VelocityDY;
        const double bc = (VelocityHeight - posB2.y - 1) * VelocityDY;
        const double bd = posB.x * VelocityDX;

        B.u = ((ba-B.x)*(B.y-bb)*LinearAccess(U, posB2.y, posB.x, VelocityHeight)+(ba-B.x)*(bc-B.y)*LinearAccess(U, posB.y, posB.x, VelocityHeight)+(B.x-bd)*(bc-B.y)*LinearAccess(U, posB.y, posB2.x, VelocityHeight)+(B.x-bd)*(B.y-bb)*LinearAccess(U, posB2.y, posB2.x, VelocityHeight))/((ba-bd)*(bc-bb));
        B.v = ((ba-B.x)*(B.y-bb)*LinearAccess(V, posB2.y, posB.x, VelocityHeight)+(ba-B.x)*(bc-B.y)*LinearAccess(V, posB.y, posB.x, VelocityHeight)+(B.x-bd)*(bc-B.y)*LinearAccess(V, posB.y, posB2.x, VelocityHeight)+(B.x-bd)*(B.y-bb)*LinearAccess(V, posB2.y, posB2.x, VelocityHeight))/((ba-bd)*(bc-bb));
    }
}

GLOBAL void UpdateParticles(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double TimeStep, const double Diffusion, const unsigned int FieldWidth, const unsigned int FieldHeight, std::uniform_real_distribution<double> &mRandom, std::mt19937_64 &gen ){
    int index_start = 0, index_stride = 1;
    #ifdef BUILD_CUDA
        index_start = blockIdx.x * blockDim.x + threadIdx.x;
        index_stride = blockDim.x * gridDim.x;
    #endif

    for(int particle = 0; particle < Particles; particle += index_stride) {
        // Particle A
        mParticleA[particle].x += mParticleA[particle].u * TimeStep + std::sqrt(2 * Diffusion * TimeStep) * mRandom(gen);
        mParticleA[particle].y += mParticleA[particle].v * TimeStep + std::sqrt(2 * Diffusion * TimeStep) * mRandom(gen);

        mParticleA[particle].x = std::fmod(mParticleA[particle].x, FieldWidth);
        mParticleA[particle].y = std::fmod(mParticleA[particle].y, FieldHeight);

        // Particle B
        mParticleB[particle].x += mParticleB[particle].u * TimeStep + std::sqrt(2 * Diffusion * TimeStep) * mRandom(gen);
        mParticleB[particle].y += mParticleB[particle].v * TimeStep + std::sqrt(2 * Diffusion * TimeStep) * mRandom(gen);

        mParticleB[particle].x = std::fmod(mParticleB[particle].x, FieldWidth);
        mParticleB[particle].y = std::fmod(mParticleB[particle].y, FieldHeight);
    }
}