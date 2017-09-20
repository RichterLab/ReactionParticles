#include "gtest/gtest.h"

#include "reaction.h"

class ReactionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
	// Setup A Particles
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

	// Setup B Particles
	mParticleB[0].x = 6.683352434689608;
	mParticleB[1].x = 8.307203074836540;
	mParticleB[2].x = 1.146952769896913;
	mParticleB[3].x = 2.803967363212206;
	mParticleB[4].x = 7.463246470539689;
	mParticleB[5].x = 4.570490854194263;
	mParticleB[6].x = 6.478712945735252;
	mParticleB[7].x = 8.659819362606868;
	mParticleB[8].x = 1.663233636699968;
	mParticleB[9].x = 7.272637561834018;

	mParticleB[0].y = 0.193934956352039;
	mParticleB[1].y = 0.410571125591421;
	mParticleB[2].y = 0.040404841745219;
	mParticleB[3].y = 0.903354372755528;
	mParticleB[4].y = 0.729667599112814;
	mParticleB[5].y = 0.245990568047725;
	mParticleB[6].y = 0.835649541816172;
	mParticleB[7].y = 0.341178356375825;
	mParticleB[8].y = 0.744876002240049;
	mParticleB[9].y = 0.955182389042570;
  }

  Particle mParticleA[10], mParticleB[10];
};

TEST_F( ReactionTest, UpdateConcentration ) {
	const size_t Particles = 10;
	const size_t FieldWidth = 10;
	const size_t FieldHeight = 1;
	const size_t ConcentrationWidth = 2 * FieldWidth;
    const size_t ConcentrationHeight = 2 * FieldHeight;
    const double ConcentrationDX = FieldWidth / (double)(ConcentrationWidth-1);
    const double ConcentrationDY = FieldHeight / (double)(ConcentrationHeight-1);

	unsigned int *CountA = (unsigned int*) malloc(sizeof(unsigned int) * ConcentrationHeight * ConcentrationWidth);
	unsigned int *CountB = (unsigned int*) malloc(sizeof(unsigned int) * ConcentrationHeight * ConcentrationWidth);

	UpdateConcentration(Particles, mParticleA, mParticleB, ConcentrationDX, ConcentrationDY, ConcentrationWidth, ConcentrationHeight, CountA, CountB);

	unsigned int CountAExpected[40] = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0
	};

	for( size_t i = 0; i < ConcentrationHeight * ConcentrationWidth; i++ ){
		ASSERT_EQ(CountA[i], CountAExpected[i]);
	}

	unsigned int CountBExpected[40] = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 2, 1, 1, 1, 1, 0, 0, 0
	};

	for( size_t i = 0; i < ConcentrationHeight * ConcentrationWidth; i++ ){
		ASSERT_EQ(CountB[i], CountBExpected[i]);
	}

}