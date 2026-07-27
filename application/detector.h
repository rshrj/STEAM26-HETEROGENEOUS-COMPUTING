#ifndef DETECTOR_H
#define DETECTOR_H

#include "main.h"

class generator;

struct detectorData
{
	unsigned int nParticles = 0;
	unsigned int nHitsCalo = 0;

	//Flat structures to copy data from / to GPUs
	particleResponse* mParticles;
	hitTracker* mHitsTracker;
	hitCalo* mHitsCalo;
	eventInfo* mEventInfos;
	eventInfo* mEventInfosCalo;
	unsigned int mNEntriesInBuffer; //particles for fastSim, trackerHits for full sim
	unsigned int mNHitsCaloInBuffer;
	unsigned int mNEventCount;
	const particleMC* mParticlesMC;
	const eventInfo* mEventInfosMC;
};

class detector : public detectorData
{
public:
	GPU() int generateResponse(const particleMC* __restrict partMC, const eventInfo& evMC, eventInfo& ev, eventInfo& evCalo);

	gpuProcessor* gpu = nullptr;

#ifdef GPU_CODE
private:
	gpuDisNormal disSigma;
#else
	int generateResponse(const generator& __restrict gen);
	detector(std::RND_TYPE& gen) : mGen(gen), disSigma(0., 1.f) {}
	detector() = delete;
	detector(const detector&) = delete;
	void ReserveMemory(unsigned int nEvents, unsigned int nParticles);

private:
	//Random variable states and types on CPU
	std::RND_TYPE& mGen;
	std::normal_distribution<float> disSigma;

	//Manage memory via std::vector on CPU
	vecpod<particleResponse> mParticlesMemory;
	vecpod<hitTracker> mHitsTrackerMemory;
	vecpod<hitCalo> mHitsCaloMemory;
	vecpod<eventInfo> mEventInfosMemory;
	vecpod<eventInfo> mEventInfosCaloMemory;
#endif
};

#endif
