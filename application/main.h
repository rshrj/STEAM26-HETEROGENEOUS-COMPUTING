#ifndef MAIN_H
#define MAIN_H

#ifdef USE_GPU
static constexpr const bool kGPUAvail = true;
#else
static constexpr const bool kGPUAvail = false;
#endif
#ifdef GPU_CODE
static constexpr const bool kGPUCode = true;
#else
static constexpr const bool kGPUCode = false;
#endif

#define NLAYERS 10
#define NSTEPS 3
enum STEPS {GENERATION = 1, DETECTOR = 2, RECONSTRUCTION = 3};
static_assert(NLAYERS >= 3, "NLAYERS must be above 3!");

#include "particles.h"
#include "gpuDefines.h"
#include "gpu.h"

#include <vector>
#include <random>
#include <cmath>
#include <memory>

#include "cmodules/qconfig.h"
#include "cmodules/vecpod.h"
#define MY_PI ((float) M_PI)

//Getters and setters for random numbers, GPU debug histograms in different versions for CPU, Vc, and GPU
//Overload "disUniform(mGen)" or alike for both float and float_v, so we can just call getRand on either of the types
#ifdef GPU_CODE
	class gpuDisPoisson {};
	class gpuDisUniform {};
	class gpuDisExponential {};
	class gpuDisNormal {};
#endif

//Type of std::RNG to use
#define RND_TYPE minstd_rand
//#define RND_TYPE mt19937_64

class randGen
{
public:
	randGen(int seed) : mGen(seed) {}
	std::RND_TYPE mGen;
};

#endif
