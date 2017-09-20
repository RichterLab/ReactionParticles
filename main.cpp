#include "reaction.h"

#include <iostream>
#include <string>
#include <fstream>
#include <random>

#include <cxxopts.hpp>

int main( int argc, char* argv[] ) {
    std::string sVelocity;
    unsigned int FieldWidth, FieldHeight, GridScale, Particles, TimeSteps;
    double TimeStep = 0.01;

    cxxopts::Options options(argv[0]);
    options.add_options()
        ("v,velocity", "Velocity file", cxxopts::value<std::string>(sVelocity)
            ->default_value("vel1.mat")->implicit_value("vel1.mat"))
        ("w,width", "Velocity field width", cxxopts::value<unsigned int>(FieldWidth)
            ->default_value("10")->implicit_value("10"))
        ("h,height", "Velocity field height", cxxopts::value<unsigned int>(FieldHeight)
            ->default_value("1")->implicit_value("1"))
        ("g,grid", "Grid Scale", cxxopts::value<unsigned int>(GridScale)
            ->default_value("5")->implicit_value("5"))
        ("p,particles", "Initial number of particles", cxxopts::value<unsigned int>(Particles)
            ->default_value("10")->implicit_value("10"))
        ("s,steps", "Total number of steps to simulate", cxxopts::value<unsigned int>(TimeSteps)
            ->default_value("5000")->implicit_value("5000"))
        ("help", "Print help");
    options.parse(argc, argv);

    if (options.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // Setup Constants
    const size_t VelocityWidth = FieldWidth * GridScale;
    const size_t VelocityHeight = FieldHeight * GridScale;
    const double VelocityDX = FieldWidth / (double)(VelocityWidth-1);
    const double VelocityDY = FieldHeight / (double)(VelocityHeight-1);

    const size_t ConcentrationWidth = 2 * FieldWidth;
    const size_t ConcentrationHeight = 2 * FieldHeight;
    const double ConcentrationDX = FieldWidth / (double)(ConcentrationWidth-1);
    const double ConcentrationDY = FieldHeight / (double)(ConcentrationHeight-1);

    const double Diffusion = 1e-5;
    const double InitialConcentration = 1.0;
    const double ReactionRate = 10.0;
    const double TimeStepEpsilon = 1.025;
    const double TimeStepMax = 0.01;
    const double ParticleMass = InitialConcentration * (FieldWidth * FieldHeight) / Particles;
    double ReactionProbability = ReactionRate * ParticleMass * TimeStep;
    const double ReactionProbabilityMax = ReactionRate * ParticleMass / 8.0 / 3.14159 / Diffusion;

    // Initialize Random Number Generator
    std::random_device rd;
    std::mt19937_64 gen(rd());
    gen.seed(1024);

    // Create Random Number Distribution
    std::uniform_real_distribution<double> mRandom;

    // Setup Velocity Fields
    double *U = (double*)malloc( sizeof(double) * VelocityWidth * VelocityHeight );
    double *V = (double*)malloc( sizeof(double) * VelocityWidth * VelocityHeight );

    // Parse Velocity File
    std::ifstream fVelocity(sVelocity, std::ifstream::in);
    for( size_t i = 0; i < VelocityHeight; i++ ){
        for( size_t j = 0; j < VelocityWidth; j++ ){
            double value;
            fVelocity >> value;
            LinearSet(U, i, j, VelocityHeight, value);
        }
    }

    for( size_t i = 0; i < VelocityHeight; i++ ){
        for( size_t j = 0; j < VelocityWidth; j++ ){
            double value;
            fVelocity >> value;
            LinearSet(V, i, j, VelocityHeight, value);
        }
    }

    fVelocity.close();

    // Setup Particles
    Particle *mParticleA = (Particle*)malloc( sizeof(Particle) * Particles );
    Particle *mParticleB = (Particle*)malloc( sizeof(Particle) * Particles );

    for( size_t i = 0; i < Particles; i++ ){
        mParticleA[i].Alive = true;
        mParticleA[i].x = FieldWidth * mRandom(gen);
        mParticleA[i].y = FieldHeight * mRandom(gen);

        mParticleB[i].Alive = true;
        mParticleB[i].x = FieldWidth * mRandom(gen);
        mParticleB[i].y = FieldHeight * mRandom(gen);
    }

    // Create Concentration Count Grid
    unsigned int *CountA = (unsigned int*) malloc(sizeof(unsigned int) * ConcentrationHeight * ConcentrationWidth);
    unsigned int *CountB = (unsigned int*) malloc(sizeof(unsigned int) * ConcentrationHeight * ConcentrationWidth);

    // Create Statistics Grids
    std::vector<double> MeanU2(TimeSteps);
    std::vector<double> MeanCA(TimeSteps);
    std::vector<double> StepConcentration(TimeSteps);
    std::vector<double> StepTime(TimeSteps);

    // Setip Temporary Arrays
    std::vector<double> ReactionChance(Particles);

    // Calculate Steps
    for( size_t step = 0; step < TimeSteps; step++ ){
        std::cout << "Step " << step << std::endl;

        TimeStep = std::min(TimeStep * TimeStepEpsilon, TimeStepMax);
        ReactionProbability = ReactionRate * ParticleMass * TimeStep;

        const double P = 0.000001;

        UpdateConcentration(Particles, mParticleA, mParticleB, ConcentrationDX, ConcentrationDY, ConcentrationWidth, ConcentrationHeight, CountA, CountB);
        for( size_t i = 0; i < ConcentrationHeight; i++ ){
            for( size_t j = 0; j < ConcentrationWidth; j++ ){
                std::cout << LinearAccess(CountA, i, j, ConcentrationHeight) << ", ";
            }
            std::cout << std::endl;
        }

        Interpolate( Particles, mParticleA, mParticleB, VelocityWidth, VelocityHeight, VelocityDX, VelocityDY, U, V );

        double U2Mean = 0.0, CAMean = 0.0;
        for( size_t i = 0; i < (ConcentrationHeight-1); i++ ){
            double PartMeanU2 = 0.0, PartMeanCA = 0.0;
            for( size_t j = 0; j < (ConcentrationWidth-1); j++ ){
                const double cas = (LinearAccess(CountA, j, i+1, ConcentrationWidth) * ParticleMass) / (ConcentrationDX * ConcentrationDY);
                const double cbs = (LinearAccess(CountB, j, i+1, ConcentrationWidth) * ParticleMass) / (ConcentrationDX * ConcentrationDY);

                PartMeanU2 += std::pow((cas - cbs), 2.0);
                PartMeanCA += cas;
            }
            U2Mean += PartMeanU2 / (ConcentrationWidth-1);
            CAMean += PartMeanCA / (ConcentrationWidth-1);
        }
        MeanU2[step] = U2Mean / (ConcentrationHeight-1);
        MeanCA[step] = CAMean / (ConcentrationHeight-1);

        std::cout << "MeanU2 " << MeanU2[step] << std::endl;
        std::cout << "MeanCA " << MeanCA[step] << std::endl;

        // Update Particle Positions
        UpdateParticles(Particles, mParticleA, mParticleB, TimeStep, Diffusion, FieldWidth, FieldHeight, mRandom, gen);
        /*for( size_t particle = 0; particle < Particles; particle++ ){
            std::cout << mParticleA[particle].x << " " << mParticleA[particle].y << std::endl;
        }*/

        // Check if all particles have reacted
        bool isComplete = true;
        for( size_t particle = 0; particle < Particles; particle++ ){
            if( mParticleA[particle].Alive ) {
                isComplete = false;
                break;
            }
        }

        if( isComplete ){
            break;
        }

        // Calculate Reactions
        for( size_t a = 0; a < Particles; a++ ){
            size_t index = 0; double maximum = -std::numeric_limits<double>::max();
            for( size_t b = 0; b < Particles; b++ ){
                const double distance = mParticleA[a].PeriodicDistance(mParticleB[b], FieldWidth, FieldHeight);
                const double probability = ReactionProbability * 1.0 / (4.0 * 3.14159 * (2.0 * Diffusion) * TimeStep) * std::exp(-std::pow(distance, 2.0) / (4.0 * (2.0 * Diffusion) * TimeStep));
                const double random = probability - mRandom(gen);

                if( random > maximum ){
                    index = b;
                    maximum = random;
                }
            }

            if( ReactionChance[index] > 0 ) {
                mParticleA[a].Alive = false;
                mParticleB[index].Alive = false;
            }
        }

        // Update Concentration Statistics
        unsigned int alive = 0;
        for( size_t a = 0; a < Particles; a++ ){
            if( mParticleA[a].Alive ) {
                alive++;
            }
        }

        StepConcentration[step+1] = ((float)alive * ParticleMass / (FieldWidth * FieldHeight));
        StepTime[step+1] = StepTime[step] + TimeStep;

        std::cout << "Step: " << step << std::endl;
        std::cout << "\tConcentration: " << StepConcentration[step+1] << std::endl;
        std::cout << "\tTotal Time: " << StepTime[step+1] << std::endl;
    }
}