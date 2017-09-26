#include "reaction.h"

#include <iostream>
#include <cfloat>

GLOBAL void UpdateConcentration(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double dx, const double dy, const size_t ConcentrationWidth, const size_t ConcentrationHeight, unsigned int *CountA, unsigned int *CountB){
    int index_start = 0, index_stride = 1;
    #ifdef BUILD_CUDA
        index_start = blockIdx.x * blockDim.x + threadIdx.x;
        index_stride = blockDim.x * gridDim.x;
    #endif

#ifndef BUILD_CUDA
    memset(CountA, 0, sizeof(unsigned int) * ConcentrationHeight * ConcentrationWidth);
    memset(CountB, 0, sizeof(unsigned int) * ConcentrationHeight * ConcentrationWidth);
#endif

    for(int idx = index_start; idx < ConcentrationHeight * ConcentrationWidth; idx += index_stride) {
        for(int particle = 0; particle < Particles; particle++) {
            const unsigned int ax = std::floor(mParticleA[particle].x / dx);
            const unsigned int ay = (ConcentrationHeight-1) - std::floor(mParticleA[particle].y / dy);
            if( ax + ay * ConcentrationWidth == idx ){
                LinearSet(CountA, ax, ay, ConcentrationWidth, LinearAccess(CountA, ax, ay, ConcentrationWidth) + 1);
            }

            const unsigned int bx = std::floor(mParticleB[particle].x / dx);
            const unsigned int by = (ConcentrationHeight-1) - std::floor(mParticleB[particle].y / dy);
            if( bx + by * ConcentrationWidth == idx ){
                LinearSet(CountB, bx, by, ConcentrationWidth, LinearAccess(CountB, bx, by, ConcentrationWidth) + 1);
            }
        }
    }
}

GLOBAL void Interpolate(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const size_t VelocityWidth, const size_t VelocityHeight, const double VelocityDX, const double VelocityDY, double *U, double *V) {
    int index_start = 0, index_stride = 1;
    #ifdef BUILD_CUDA
        index_start = blockIdx.x * blockDim.x + threadIdx.x;
        index_stride = blockDim.x * gridDim.x;
    #endif

    for(int particle = index_start; particle < Particles; particle += index_stride) {
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

#ifdef BUILD_CUDA
GLOBAL void InitializeRandom(unsigned int seed, curandState_t* states) {
      curand_init(seed, /* the seed can be the same for each core, here we pass the time in from the CPU */
                  blockIdx.x, /* the sequence number should be different for each core (unless you want all
                                 cores to get the same sequence of numbers for some reason - use thread id! */
                  0, /* the offset is how much extra we advance in the sequence for each call, can be 0 */
                  &states[blockIdx.x]);
}
#endif

#ifdef BUILD_CUDA
GLOBAL void UpdateParticles(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double TimeStep, const double Diffusion, const unsigned int FieldWidth, const unsigned int FieldHeight, curandState_t* states ){
#else
void UpdateParticles(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double TimeStep, const double Diffusion, const unsigned int FieldWidth, const unsigned int FieldHeight, std::uniform_real_distribution<double> &mRandom, std::mt19937_64 &gen ){
#endif
    int index_start = 0, index_stride = 1;
    #ifdef BUILD_CUDA
        index_start = blockIdx.x * blockDim.x + threadIdx.x;
        index_stride = blockDim.x * gridDim.x;
    #endif

    #ifdef BUILD_CUDA
    #define RANDOM curand_uniform(&states[blockIdx.x])
    #else
    #define RANDOM mRandom(gen)
    #endif

    for(int particle = index_start; particle < Particles; particle += index_stride) {
        // Particle A
        mParticleA[particle].x += mParticleA[particle].u * TimeStep + std::sqrt(2 * Diffusion * TimeStep) * RANDOM;
        mParticleA[particle].y += mParticleA[particle].v * TimeStep + std::sqrt(2 * Diffusion * TimeStep) * RANDOM;

        mParticleA[particle].x = std::fmod((float)mParticleA[particle].x, (float)FieldWidth);
        mParticleA[particle].y = std::fmod((float)mParticleA[particle].y, (float)FieldHeight);

        // Particle B
        mParticleB[particle].x += mParticleB[particle].u * TimeStep + std::sqrt(2 * Diffusion * TimeStep) * RANDOM;
        mParticleB[particle].y += mParticleB[particle].v * TimeStep + std::sqrt(2 * Diffusion * TimeStep) * RANDOM;

        mParticleB[particle].x = std::fmod((float)mParticleB[particle].x, (float)FieldWidth);
        mParticleB[particle].y = std::fmod((float)mParticleB[particle].y, (float)FieldHeight);
    }
}

#ifdef BUILD_CUDA
GLOBAL void UpdateReactions(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double TimeStep, const double Diffusion, const double ReactionProbability, const unsigned int FieldWidth, const unsigned int FieldHeight, curandState_t* states){
#else
void UpdateReactions(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double TimeStep, const double Diffusion, const double ReactionProbability, const unsigned int FieldWidth, const unsigned int FieldHeight, std::uniform_real_distribution<double> &mRandom, std::mt19937_64 &gen){
#endif
    int index_start = 0, index_stride = 1;
    #ifdef BUILD_CUDA
        index_start = blockIdx.x * blockDim.x + threadIdx.x;
        index_stride = blockDim.x * gridDim.x;
    #endif

    #ifdef BUILD_CUDA
    #define RANDOM curand_uniform(&states[blockIdx.x])
    #else
    #define RANDOM mRandom(gen)
    #endif

    const double P = 0.000001;
    const double Cutoff = std::sqrt( -8.0 * Diffusion * TimeStep * std::log( 8.0 * 3.14159 * Diffusion * TimeStep * P / ReactionProbability));

    for(int a = index_start; a < Particles; a += index_stride) {
        size_t index = 0; double maximum = -DBL_MAX;
        for( size_t b = 0; b < Particles; b++ ){
            const double distance = mParticleA[a].PeriodicDistance(mParticleB[b], FieldWidth, FieldHeight);
            if( distance < Cutoff ){
                const double probability = ReactionProbability * 1.0 / (4.0 * 3.14159 * (2.0 * Diffusion) * TimeStep) * std::exp(-std::pow(distance, 2.0) / (4.0 * (2.0 * Diffusion) * TimeStep));
                const double random = probability - RANDOM;

                if( random > maximum ){
                    index = b;
                    maximum = random;
                }
            }
        }

        if( maximum > 0 ) {
            mParticleA[a].Alive = false;
            mParticleB[index].Alive = false;
        }
    }
}