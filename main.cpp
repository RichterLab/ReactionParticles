#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <random>

#include <cxxopts.hpp>

struct Index {
    unsigned int x, y;

    Index(unsigned int x, unsigned int y) : x(x), y(y) {

    }
};

struct Particle{
    double x, y, u, v, mass;
};

struct Field {
    unsigned int Width, Height;
    double dx, dy;
    std::vector<std::vector<double>> data;

    Field(unsigned int width, unsigned int height, double xLength, double yLength) : Width(width), Height(height) {
        data.resize(height);
        for( size_t i = 0; i < height; i++ ){
            data[i].resize(width);
        }

        dx = xLength/(width-1);
        dy = yLength/(height-1);
    }

    const Index GetIndex(Particle &p) {
        return Index( std::floor(p.x/dx), (Height-1) - std::floor(p.y/dy) );
    }
};

int main( int argc, char* argv[] ) {
    std::string sVelocity;
    unsigned int FieldWidth, FieldHeight, GridScale, Particles, TimeSteps = 1;
    std::vector<std::vector<double>> U, V;
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
        ("help", "Print help");
    options.parse(argc, argv);

    if (options.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // Initialize Random Number Generator
    std::random_device rd;
    std::mt19937_64 gen(rd());

    // Resize Velocity Fields
    U.resize(FieldHeight * GridScale);
    for( size_t i = 0; i < FieldHeight * GridScale; i++ ){
        U[i].resize(FieldWidth * GridScale);
    }

    V.resize(FieldHeight * GridScale);
    for( size_t i = 0; i < FieldHeight * GridScale; i++ ){
        V[i].resize(FieldWidth * GridScale);
    }

    // Parse Velocity File
    std::ifstream fVelocity(sVelocity, std::ifstream::in);
    for( size_t i = 0; i < FieldHeight * GridScale; i++ ){
        for( size_t j = 0; j < FieldWidth * GridScale; j++ ){
            fVelocity >> U[i][j];
        }
    }

    for( size_t i = 0; i < FieldHeight * GridScale; i++ ){
        for( size_t j = 0; j < FieldWidth * GridScale; j++ ){
            fVelocity >> V[i][j];
        }
    }

    fVelocity.close();

    // Setup Constants
    const double Diffusion = 1e-5;
    const double InitialConcentration = 1.0;
    const double ReactionRate = 10.0;
    const double TimeStepEpsilon = 1.025;
    const double TimeStepMax = 0.01;
    const double ParticleMass = InitialConcentration * (FieldWidth * FieldHeight) / Particles;
    double ReactionProbability = ReactionRate * ParticleMass * TimeStep;
    const double ReactionProbabilityMax = ReactionRate * ParticleMass / 8.0 / 3.14159 / Diffusion;

    // Create Random Number Distribution
    std::uniform_real_distribution<double> mRandom;

    // Setup Particles
    std::vector<Particle> mParticleA(Particles), mParticleB(Particles);

    for( size_t i = 0; i < Particles; i++ ){
        mParticleA[i].x = FieldWidth * mRandom(gen);
        mParticleA[i].y = FieldHeight * mRandom(gen);
        mParticleA[i].mass = ParticleMass;

        mParticleB[i].x = FieldWidth * mRandom(gen);
        mParticleB[i].y = FieldHeight * mRandom(gen);
        mParticleB[i].mass = ParticleMass;
    }

    // Create Concentration Grid
    Field Concentration(2 * FieldWidth, 2 * FieldHeight, FieldWidth, FieldHeight);

    // Create Concentration Count Grid
    std::vector<std::vector<unsigned int>> CountA(2 * FieldHeight), CountB(2 * FieldHeight);
    for( size_t i = 0; i < 2 * FieldHeight; i++ ){
        CountA[i].resize(2 * FieldWidth);
        CountB[i].resize(2 * FieldWidth);
    }

    for( size_t step = 0; step < TimeSteps; step++ ){
        std::cout << "Step " << step << std::endl;

        TimeStep = std::min(TimeStep * TimeStepEpsilon, TimeStepMax);
        ReactionProbability = ReactionRate * ParticleMass * TimeStep;

        const double P = 0.000001;

        // Update Particle Concentration Count
        for( size_t particle = 0; particle < Particles; particle++ ){
            const Index posA = Concentration.GetIndex(mParticleA[particle]);
            CountA[posA.y][posA.x] += 1;

            const Index posB = Concentration.GetIndex(mParticleB[particle]);
            CountB[posB.y][posB.x] += 1;
        }

        // Interpolate Velocities
        for( size_t particle = 0; particle < Particles; particle++ ){
        }
    }
}