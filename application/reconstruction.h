#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include "main.h"

class detector;
template <class T> struct hitList_T;
typedef hitList_T<float> hitList;

struct reconstructionData
{
	unsigned int nParticles = 0;
	
	//Flat structures to copy data from / to GPUs
	particleResponse* mParticles;
	eventInfo* mEventInfos;
	unsigned int mNEntriesInBuffer;
	unsigned int mParticleBufferSize;
	unsigned int mTmpBufferSize;
	unsigned int mNEventCount;
	const particleMC* mParticlesMC;
	const eventInfo* mEventInfosMC;
	char* mTmpBuffer;
};

class reconstruction : public reconstructionData
{
public:
	GPU() int reconstruct(const hitTracker* __restrict hitsTracker, int nHitsTracker, const hitCalo* __restrict hitsCalo, int nHitsCalo, const eventInfo& evTracker, const eventInfo& evCalo, eventInfo& ev);
	GPU() void fit(const eventInfo& ev, const hitTracker* __restrict hitsTracker, const hitCalo* __restrict hitsCalo, const particleMC* __restrict particlesMC);
	template <class T> GPUi() void fitParticle(particleResponse_T<T>& part, hitList_T<T>& partHits, const hitTracker* __restrict hitsTracker, const hitCalo* __restrict hitsCalo, T& lastR);
	GPUi() void fitParticleCalo(particleResponse& part, hitList& partHits, const hitCalo* __restrict hitsCalo, float& lastR);
	GPUh() void estimateParticleBufferSize(unsigned int nHitsTracker);
	GPUh() void calculateScratchBufferSize(unsigned int nParticles, unsigned int nHitsTracker, unsigned int nHitsCalo);
	
	int mError;
	gpuProcessor* gpu = nullptr;
	
#ifdef GPU_CODE
#else
	reconstruction() {}
	reconstruction(const reconstruction&) = delete;
	void ReserveMemory(unsigned int nEvents, unsigned int nParticles, unsigned int nHitsTracker, unsigned int nHitsCalo);
	int reconstruct(const detector& __restrict det);

private:
	//Manage memory via std::vector on CPU
	vecpod<particleResponse> mParticlesMemory;
	vecpod<eventInfo> mEventInfosMemory;
	vecpod<char> mTmpMem;
#endif
};

#endif
