#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <random>

#include <cxxopts.hpp>

typedef std::vector<std::vector<double>> DoubleGrid;

DoubleGrid CreateDoubleGrid(size_t width, size_t height) {
    DoubleGrid retVal(height);
    for( size_t i = 0; i < height; i++ ){
        retVal[i].resize(width);
    }
    return retVal;
}

struct Index {
    int x, y;
    Index(int x, int y) : x(x), y(y) {}
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

struct Field {
    unsigned int Width, Height;
    double dx, dy;
    DoubleGrid data;
    std::vector<std::vector<std::pair<double,double>>> steps;

    Field(unsigned int width, unsigned int height, double xLength, double yLength) : Width(width), Height(height) {
        data = CreateDoubleGrid(width, height);

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

    // Initialize Random Number Generator
    std::random_device rd;
    std::mt19937_64 gen(rd());

    // Setup Velocity Fields
    DoubleGrid U = CreateDoubleGrid(FieldWidth * GridScale, FieldHeight * GridScale);
    DoubleGrid V = CreateDoubleGrid(FieldWidth * GridScale, FieldHeight * GridScale);

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
        mParticleA[i].Alive = true;
        mParticleB[i].Alive = true;
        mParticleA[i].x = FieldWidth * mRandom(gen);
        mParticleA[i].y = FieldHeight * mRandom(gen);

        mParticleB[i].x = FieldWidth * mRandom(gen);
        mParticleB[i].y = FieldHeight * mRandom(gen);
    }

    // Create Concentration Grid
    Field Concentration(2 * FieldWidth, 2 * FieldHeight, FieldWidth, FieldHeight);

    // Create Concentration Count Grid
    std::vector<std::vector<unsigned int>> CountA(2 * FieldHeight), CountB(2 * FieldHeight);
    for( size_t i = 0; i < 2 * FieldHeight; i++ ){
        CountA[i].resize(2 * FieldWidth);
        CountB[i].resize(2 * FieldWidth);
    }

    DoubleGrid CountAS = CreateDoubleGrid((2 * FieldWidth)-1, (2 * FieldHeight)-1);
    DoubleGrid CountBS = CreateDoubleGrid((2 * FieldWidth)-1, (2 * FieldHeight)-1);

    // Create Velocity Grid
    Field Velocity(U[0].size(), U.size(), FieldWidth, FieldHeight);

    // Create Statistics Grids
    DoubleGrid DiffSquared = CreateDoubleGrid((2 * FieldWidth)-1, (2 * FieldHeight)-1);

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
            // Interpolate Particle A
            Particle& A = mParticleA[particle];
            const Index posA = Velocity.GetIndex(A);
            Index posA2(posA.x+1, posA.y-1);

            if (posA2.x == Velocity.Width) posA2.x = 0;
            if (posA2.y == -1) posA2.y = Velocity.Height-1;

            A.u = ((Grid[1][posA2.x].first-A.x)*(A.y-Grid[posA.y][1].second)*U[posA2.y][posA.x]+(Grid[1][posA2.x].first-A.x)*(Grid[posA2.y][1].second-A.y)*U[posA.y][posA.x]+(A.x-Grid[1][posA.x].first)*(Grid[posA2.y][1].second-A.y)*U[posA.y][posA2.x]+(A.x-Grid[1][posA.x].first)*(A.y-Grid[posA.y][1].second)*U[posA2.y][posA2.x])/((Grid[1][posA2.x].first-Grid[1][posA.x].first)*(Grid[posA2.y][1].second-Grid[posA.y][1].second));
            A.v = ((Grid[1][posA2.x].first-A.x)*(A.y-Grid[posA.y][1].second)*V[posA2.y][posA.x]+(Grid[1][posA2.x].first-A.x)*(Grid[posA2.y][1].second-A.y)*V[posA.y][posA.x]+(A.x-Grid[1][posA.x].first)*(Grid[posA2.y][1].second-A.y)*V[posA.y][posA2.x]+(A.x-Grid[1][posA.x].first)*(A.y-Grid[posA.y][1].second)*V[posA2.y][posA2.x])/((Grid[1][posA2.x].first-Grid[1][posA.x].first)*(Grid[posA2.y][1].second-Grid[posA.y][1].second));

            // Interpolate Particle B
            Particle& B = mParticleB[particle];
            const Index posB = Velocity.GetIndex(B);
            Index posB2(posB.x+1, posB.y-1);

            if (posB2.x == Velocity.Width) posB2.x = 0;
            if (posB2.y == -1) posB2.y = Velocity.Height-1;

            B.u = ((Grid[1][posB2.x].first-B.x)*(B.y-Grid[posB.y][1].second)*U[posB2.y][posB.x]+(Grid[1][posB2.x].first-B.x)*(Grid[posB2.y][1].second-B.y)*U[posB.y][posB.x]+(B.x-Grid[1][posB.x].first)*(Grid[posB2.y][1].second-B.y)*U[posB.y][posB2.x]+(B.x-Grid[1][posB.x].first)*(B.y-Grid[posB.y][1].second)*U[posB2.y][posB2.x])/((Grid[1][posB2.x].first-Grid[1][posB.x].first)*(Grid[posB2.y][1].second-Grid[posB.y][1].second));
            B.v = ((Grid[1][posB2.x].first-B.x)*(B.y-Grid[posB.y][1].second)*V[posB2.y][posB.x]+(Grid[1][posB2.x].first-B.x)*(Grid[posB2.y][1].second-B.y)*V[posB.y][posB.x]+(B.x-Grid[1][posB.x].first)*(Grid[posB2.y][1].second-B.y)*V[posB.y][posB2.x]+(B.x-Grid[1][posB.x].first)*(B.y-Grid[posB.y][1].second)*V[posB2.y][posB2.x])/((Grid[1][posB2.x].first-Grid[1][posB.x].first)*(Grid[posB2.y][1].second-Grid[posB.y][1].second));
        }

        for( size_t i = 0; i < CountAS.size(); i++ ){
            for( size_t j = 0; j < CountAS[0].size(); j++ ){
                CountAS[i][j] = (CountA[i+1][j] * ParticleMass) / (Concentration.dx * Concentration.dy);
                CountBS[i][j] = (CountB[i+1][j] * ParticleMass) / (Concentration.dx * Concentration.dy);
            }
        }

        double U2Mean = 0.0;
        for( size_t i = 0; i < CountAS.size(); i++ ){
            double PartMean = 0.0;
            for( size_t j = 0; j < CountAS[0].size(); j++ ){
                PartMean += std::pow(CountAS[i][j] - CountBS[i][j], 2);
            }
            U2Mean += PartMean / CountAS[0].size();
        }
        MeanU2[step] = U2Mean / CountAS.size();

        double CAMean = 0.0;
        for( size_t i = 0; i < CountAS.size(); i++ ){
            double PartMean = 0.0;
            for( size_t j = 0; j < CountAS[0].size(); j++ ){
                PartMean += CountAS[i][j];
            }
            CAMean += PartMean / CountAS[0].size();
        }
        MeanCA[step] = CAMean / CountAS.size();

        // Update Particle Positions
        for( size_t particle = 0; particle < Particles; particle++ ){
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
            for( size_t b = 0; b < Particles; b++ ){
                const double distance = mParticleA[a].PeriodicDistance(mParticleB[b], FieldWidth, FieldHeight);
                const double probability = ReactionProbability * 1.0 / (4.0 * 3.14159 * (2.0 * Diffusion) * TimeStep) * std::exp(-std::pow(distance, 2.0) / (4.0 * (2.0 * Diffusion) * TimeStep));
                const double random = probability - mRandom(gen);

                ReactionChance[b] = random;
            }

            size_t index = 0; double maximum = std::numeric_limits<double>::min();
            for( size_t b = 0; b < Particles; b++ ){
                if( ReactionChance[b] > maximum ){
                    index = b;
                    maximum = ReactionChance[b];
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