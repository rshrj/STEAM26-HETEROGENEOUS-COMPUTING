#include "detector.h"
#include "generator.h"

template <class T, class S, class U> GPUi() static int detectorResponsePos(T& resp, int mcId, S& disSigma, U& mGen, const particleMC& part, float* pos, float R)
{
	if (part.intersectWithRadius(R, pos)) return(1); //Outside of detector R
	if (std::abs(pos[2]) > config.detectorZ) return(1); //Outside of detector Z
	resp.mMCId = mcId;
	float sigma = 0.f;
	sigma = disSigma(mGen);
	resp.mX[2] = pos[2] + sigma * config.resZ; //1 mm z resolution
	float phi = std::atan2(pos[1], pos[0]);
	sigma = disSigma(mGen);
	phi += sigma * config.resRPhi / R; //1 mm rPhi resolution
	float sinPhi = std::sin(phi);
	float cosPhi = std::cos(phi); //no sqrt(1 - ...)!! --> numerical issues
	resp.mX[0] = cosPhi * R;
	resp.mX[1] = sinPhi * R;
	return(0);
}

GPUi() int detector::generateResponse(const particleMC* __restrict partMC, const eventInfo& evMC, eventInfo& ev, eventInfo& evCalo)
{
	unsigned int nResponse = 0;
	unsigned int nResponseCalo = 0;
	ev.mFirstParticle = nParticles;
	evCalo.mFirstParticle =  nHitsCalo;

	for (unsigned int i = 0;i < evMC.mNParticles;i++)
	{
		const particleMC& part = partMC[i];
		if (config.particles[part.mPID].mDecayLength > 0.f) continue;
		float pos[3];
		if (config.particles[part.mPID].hitTracker())
		{
			for (int j = 0;j < NLAYERS;j++)
			{
				hitTracker& resp = mHitsTracker[ev.mFirstParticle + nResponse];
				if (detectorResponsePos(resp, i, disSigma, mGen, part, pos, config.detectorR[j])) continue;
				nResponse++;
				resp.mLayer = j;
			}
		}
		if (config.particles[part.mPID].hitCalo())
		{
			hitCalo& resp = mHitsCalo[evCalo.mFirstParticle + nResponseCalo];
			if (detectorResponsePos(resp, i, disSigma, mGen, part, pos, config.detectorRCalo)) continue;
			nResponseCalo++;
			float sigma = 0.f;
			sigma = disSigma(mGen);
			float p = std::sqrt(part.p2() + config.particles[part.mPID].mMass * config.particles[part.mPID].mMass);
			p = p + sigma * config.resPt * p * p;
			if (p < 0) p = 0;
			if (config.particles[part.mPID].hitCaloHadron()) p = -p;//Use negative sign to flag hadronic particle
			resp.mEnergy = p;
		}
	}
	ev.mNParticles = nResponse;
	evCalo.mNParticles = nResponseCalo;
	nParticles += nResponse;
	nHitsCalo += nResponseCalo;

	return(0);
}

void detector::ReserveMemory(unsigned int nEvents, unsigned int nParticles)
{
	if (nParticles > mParticlesMemory.size())
	{
		mHitsTrackerMemory.resize(nParticles * NLAYERS);
		mHitsCaloMemory.resize(nParticles);
	}
	if (nEvents > mEventInfosMemory.size())
	{
		mEventInfosMemory.resize(nEvents);
		mEventInfosCaloMemory.resize(nEvents);
	}
	mParticles = mParticlesMemory.data();
	mHitsTracker = mHitsTrackerMemory.data();
	mHitsCalo = mHitsCaloMemory.data();
	mEventInfos = mEventInfosMemory.data();
	mEventInfosCalo = mEventInfosCaloMemory.data();
}

int detector::generateResponse(const generator& __restrict gen)
{
	ReserveMemory(gen.mNEventCount, gen.mNParticles);
	mNEventCount = gen.mNEventCount;
	mParticlesMC = gen.mParticlesMC;
	mEventInfosMC = gen.mEventInfosMC;
	nParticles = 0;
	nHitsCalo = 0;
	for (unsigned int i = 0;i < gen.mNEventCount;i++)
	{
		if (generateResponse(&mParticlesMC[mEventInfosMC[i].mFirstParticle], mEventInfosMC[i], mEventInfos[i], mEventInfosCalo[i])) return(1);
	}
	return(0);
}
