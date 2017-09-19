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

template <typename T>
const T LinearAccess(T* array, const size_t x, const size_t y, const size_t Width) {
    return array[x + y * Width];
}

template <typename T>
void LinearSet(T* array, const size_t x, const size_t y, const size_t Width, const T value) {
    array[x + y * Width] = value;
}

void UpdateConcentration(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double dx, const double dy, const unsigned int ConcentrationHeight, unsigned int *CountA, unsigned int *CountB, size_t CountWidth){
    for( size_t particle = 0; particle < Particles; particle++ ){
        const unsigned int ax = std::floor(mParticleA[particle].x / dx);
        const unsigned int ay = (ConcentrationHeight-1) - std::floor(mParticleA[particle].y / dy);
        LinearSet(CountA, ay, ax, CountWidth, LinearAccess(CountA, ay, ax, CountWidth) + 1);

        const unsigned int bx = std::floor(mParticleB[particle].x / dx);
        const unsigned int by = (ConcentrationHeight-1) - std::floor(mParticleB[particle].y / dy);
        LinearSet(CountB, by, bx, CountWidth, LinearAccess(CountB, by, bx, CountWidth) + 1);
    }
}

void Interpolate(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const size_t VelocityWidth, const size_t VelocityHeight, const double VelocityDX, const double VelocityDY, double *U, double *V) {
    for( size_t particle = 0; particle < Particles; particle++ ){
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

void UpdateParticles(const size_t Particles, Particle *mParticleA, Particle *mParticleB, const double TimeStep, const double Diffusion, const unsigned int FieldWidth, const unsigned int FieldHeight, std::uniform_real_distribution<double> &mRandom, std::mt19937_64 &gen ){
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
}

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
    memset(CountA, 0, sizeof(unsigned int) * ConcentrationHeight * ConcentrationWidth);

    unsigned int *CountB = (unsigned int*) malloc(sizeof(unsigned int) * ConcentrationHeight * ConcentrationWidth);
    memset(CountB, 0, sizeof(unsigned int) * ConcentrationHeight * ConcentrationWidth);

    DoubleGrid CountAS = CreateDoubleGrid(ConcentrationWidth-1, ConcentrationHeight-1);
    DoubleGrid CountBS = CreateDoubleGrid(ConcentrationWidth-1, ConcentrationHeight-1);

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

        UpdateConcentration(Particles, mParticleA, mParticleB, ConcentrationDX, ConcentrationDY, ConcentrationHeight, CountA, CountB, ConcentrationHeight);
        for( size_t i = 0; i < ConcentrationHeight; i++ ){
            for( size_t j = 0; j < ConcentrationWidth; j++ ){
                std::cout << LinearAccess(CountA, i, j, ConcentrationHeight) << ", ";
            }
            std::cout << std::endl;
        }

        Interpolate( Particles, mParticleA, mParticleB, VelocityWidth, VelocityHeight, VelocityDX, VelocityDY, U, V );

        double U2Mean = 0.0, CAMean = 0.0;
        for( size_t i = 0; i < CountAS.size(); i++ ){
            double PartMeanU2 = 0.0, PartMeanCA = 0.0;
            for( size_t j = 0; j < CountAS[0].size(); j++ ){
                CountAS[i][j] = (LinearAccess(CountA, i+1, j, ConcentrationWidth) * ParticleMass) / (ConcentrationDX * ConcentrationDY);
                CountBS[i][j] = (LinearAccess(CountB, i+1, j, ConcentrationWidth) * ParticleMass) / (ConcentrationDX * ConcentrationDY);

                PartMeanU2 += std::pow(CountAS[i][j] - CountBS[i][j], 2);
                PartMeanCA += CountAS[i][j];
            }
            U2Mean += PartMeanU2 / CountAS[0].size();
            CAMean += PartMeanCA / CountAS[0].size();
        }
        MeanU2[step] = U2Mean / CountAS.size();
        MeanCA[step] = CAMean / CountAS.size();

        std::cout << "MeanU2 " << MeanU2[step] << std::endl;
        std::cout << "MeanCA " << MeanCA[step] << std::endl;

        // Update Particle Positions
        UpdateParticles(Particles, mParticleA, mParticleB, TimeStep, Diffusion, FieldWidth, FieldHeight, mRandom, gen);

        for( size_t particle = 0; particle < Particles; particle++ ){
            std::cout << mParticleA[particle].x << " " << mParticleA[particle].y << std::endl;
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