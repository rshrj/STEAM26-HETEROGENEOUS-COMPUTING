#ifndef GENERATOR_H
#define GENERATOR_H

#include "main.h"

struct generatorData
{
	unsigned int nDecayParticles = 0;

	//Flat structures to copy data from / to GPUs
	particleMC* mParticlesMC;
	eventInfo* mEventInfosMC;
	unsigned int mParticleBufferSize;
	unsigned int mNParticles;
	unsigned int mNEventCount;
};

class generator : public generatorData
{
public:
	GPUh() void estimateBufferSize(int nEvents);
	GPU() int generateEvent(int evNum);
	GPU() float genRandPower(float onset, int power);
	GPU() void generatePrimParticle(particleMC& part, float disPowerMin, float disPower, unsigned short pid, float maxEta);
	template <class T> GPU() int decay3Body(particleMC_T<T> &part, vec3_T<T>* decayProducts);

	gpuProcessor* gpu = nullptr;

#ifndef GPU_CODE
	generator(std::RND_TYPE& gen);
	generator() = delete;
	generator(const generator&) = delete;
	void ReserveMemory(unsigned int nEvents, unsigned int nParticles);
	int generateEvents(int nEvents);

private:
	//Random variable states and types on CPU
	std::RND_TYPE& mGen;
	std::vector<std::poisson_distribution<int>> disPoisson;
	std::exponential_distribution<float> disExp;
	std::uniform_real_distribution<float> disUniform;

	//Manage memory via vecpod (std::vector w.o. calling constructors for initialization) on CPU
	vecpod<particleMC> mParticlesMCMemory;
	vecpod<eventInfo> mEventInfosMemory;
#endif
private:
	GPUi() void insertSecondaries(particleMC& part, vec3* decayProducts, particleMC* nextParticle, unsigned int& nParticle);
};

#endif
