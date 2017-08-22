#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <random>

#include <cxxopts.hpp>

struct Index {
    int x, y;

    Index(int x, int y) : x(x), y(y) {

    }
};

struct Particle{
    double x, y, u, v, mass;
};

struct Field {
    unsigned int Width, Height;
    double dx, dy;
    std::vector<std::vector<double>> data;
    std::vector<std::vector<std::pair<double,double>>> steps;

    Field(unsigned int width, unsigned int height, double xLength, double yLength) : Width(width), Height(height) {
        data.resize(height);
        for( size_t i = 0; i < height; i++ ){
            data[i].resize(width);
        }

        dx = xLength/(width-1);
        dy = yLength/(height-1);

        steps.resize(height);
        for( size_t i = 0; i < height; i++ ){
            steps[i].resize(width);
            for( size_t j = 0; j < width; j++ ){
                steps[i][j] = std::make_pair(j * dx, (height-i-1) * dy);
            }
        }
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

    mParticleA[0].x = 6.9234;
    mParticleA[1].x = 4.3310;
    mParticleA[2].x = 4.2018;
    mParticleA[3].x = 1.4481;
    mParticleA[4].x = 5.1285;
    mParticleA[5].x = 8.0339;
    mParticleA[6].x = 9.3188;
    mParticleA[7].x = 6.4324;
    mParticleA[8].x = 5.5442;
    mParticleA[9].x = 1.1989;

    mParticleA[0].y = 0.1467;
    mParticleA[1].y = 0.1462;
    mParticleA[2].y = 0.9986;
    mParticleA[3].y = 0.0473;
    mParticleA[4].y = 0.6285;
    mParticleA[5].y = 0.6246;
    mParticleA[6].y = 0.2868;
    mParticleA[7].y = 0.2002;
    mParticleA[8].y = 0.4078;
    mParticleA[9].y = 0.9763;

    /*for( size_t i = 0; i < Particles; i++ ){
        mParticleA[i].x = FieldWidth * mRandom(gen);
        mParticleA[i].y = FieldHeight * mRandom(gen);
        mParticleA[i].mass = ParticleMass;

        mParticleB[i].x = FieldWidth * mRandom(gen);
        mParticleB[i].y = FieldHeight * mRandom(gen);
        mParticleB[i].mass = ParticleMass;
    }*/

    // Create Concentration Grid
    Field Concentration(2 * FieldWidth, 2 * FieldHeight, FieldWidth, FieldHeight);

    // Create Concentration Count Grid
    std::vector<std::vector<unsigned int>> CountA(2 * FieldHeight), CountB(2 * FieldHeight);
    for( size_t i = 0; i < 2 * FieldHeight; i++ ){
        CountA[i].resize(2 * FieldWidth);
        CountB[i].resize(2 * FieldWidth);
    }

    std::vector<std::vector<double>> CountAS((2 * FieldHeight)-1), CountBS((2 * FieldHeight)-1);
    for( size_t i = 0; i < (2 * FieldHeight)-1; i++ ){
        CountAS[i].resize((2 * FieldWidth)-1);
        CountBS[i].resize((2 * FieldWidth)-1);
    }

    // Create Velocity Grid
    Field Velocity(U[0].size(), U.size(), FieldWidth, FieldHeight);

    // Calculate Steps
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
        const std::vector<std::vector<std::pair<double,double>>>& Grid = Velocity.steps;

        for( size_t particle = 0; particle < Particles; particle++ ){
            Particle& A = mParticleA[particle];
            const Index posA = Velocity.GetIndex(A);
            Index posA2(posA.x+1, posA.y-1);

            if (posA2.x == Velocity.Width) posA2.x = 1;
            if (posA2.y == -1) posA2.y = Velocity.Height-1;

            A.u = ((Grid[1][posA2.x].first - A.x) * (A.y - Grid[posA.y][1].second) * U[posA2.y][posA.x] + (Grid[1][posA2.x].first - A.x) *(Grid[posA2.y][1].second - A.y) * U[posA.y][posA.x] + (A.x - Grid[1][posA.x].first) * (Grid[posA2.y][1].second - A.y) *U[posA.y][posA2.x] + (A.x - Grid[1][posA.x].first) * (A.y - Grid[posA.y][1].second) * U[posA2.y][posA2.x]) / ((Grid[1][posA2.x].first - Grid[1][posA.x].first) * (Grid[posA2.y][1].second - Grid[posA.y][1].second));

            A.v = ((Grid[1][posA2.x].first-A.x)*(A.y-Grid[posA.y][1].second)*V[posA2.y][posA.x]+(Grid[1][posA2.x].first-A.x)*(Grid[posA2.y][1].second-A.y)*V[posA.y][posA.x]+(A.x-Grid[1][posA.x].first)*(Grid[posA2.y][1].second-A.y)*V[posA.y][posA2.x]+(A.x-Grid[1][posA.x].first)*(A.y-Grid[posA.y][1].second)*V[posA2.y][posA2.x])/((Grid[1][posA2.x].first-Grid[1][posA.x].first)*(Grid[posA2.y][1].second-Grid[posA.y][1].second));

            Particle& B = mParticleB[particle];
            const Index posB = Velocity.GetIndex(B);
            Index posB2(posB.x+1, posB.y-1);

            if (posB2.x == Velocity.Width) posB2.x = 1;
            if (posB2.y == -1) posB2.y = Velocity.Height-1;

            B.u = ((Grid[1][posB2.x].first-B.x)*(B.y-Grid[posB.y][1].second)*U[posB2.y][posB.x]+(Grid[1][posB2.x].first-B.x)*(Grid[posB2.y][1].second-B.y)*U[posB.y][posB.x]+(B.x-Grid[1][posB.x].first)*(Grid[posB2.y][1].second-B.y)*U[posB.y][posB2.x]+(B.x-Grid[1][posB.x].first)*(B.y-Grid[posB.y][1].second)*U[posB2.y][posB2.x])/((Grid[1][posB2.x].first-Grid[1][posB.x].first)*(Grid[posB2.y][1].second-Grid[posB.y][1].second));

            B.v = ((Grid[1][posB2.x].first-B.x)*(B.y-Grid[posB.y][1].second)*V[posB2.y][posB.x]+(Grid[1][posB2.x].first-B.x)*(Grid[posB2.y][1].second-B.y)*V[posB.y][posB.x]+(B.x-Grid[1][posB.x].first)*(Grid[posB2.y][1].second-B.y)*V[posB.y][posB2.x]+(B.x-Grid[1][posB.x].first)*(B.y-Grid[posB.y][1].second)*V[posB2.y][posB2.x])/((Grid[1][posB2.x].first-Grid[1][posB.x].first)*(Grid[posB2.y][1].second-Grid[posB.y][1].second));
        }
    }
}