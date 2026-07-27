#include "detector.h"
#include "reconstruction.h"

template <class T> struct indexType;
template <> struct indexType<float> {typedef int type;};

template <class T> struct hitList_T
{
	typename indexType<T>::type hitId[NLAYERS + 1];
};
typedef hitList_T<float> hitList;

GPUi() int reconstruction::reconstruct(const hitTracker* __restrict hitsTracker, int nHitsTracker, const hitCalo* __restrict hitsCalo, int nHitsCalo, const eventInfo& evTracker, const eventInfo& evCalo, eventInfo& ev)
{
	ev.mFirstParticle = kGPUCode ? evCalo.mFirstParticle : mNEntriesInBuffer; //For GPUs, we reserve space for 1 calo hit, so we take that index here as it matches the space for 1 track
	unsigned int partFound = 0;
	//MC based stupid reconstruction algorithm
	int nCurMCParticle = 0;
	unsigned int nTotalHitsCovered = 0;
	while (nTotalHitsCovered < evTracker.mNParticles)
	{
		int nFound = 0;
		int recHits[NLAYERS];
		for (int i = 0;i < nHitsTracker;i++)
		{
			if (hitsTracker[i].mMCId == nCurMCParticle)
			{
				recHits[nFound++] = i;
				if (nFound >= NLAYERS) break;
			}
		}
		nTotalHitsCovered += nFound;
		if (nFound == NLAYERS)
		{
			int caloHit = -1;
			for (int i = 0;i < nHitsCalo;i++)
			{
				if (hitsCalo[i].mMCId == nCurMCParticle)
				{
					caloHit = i;
					break;
				}
			}
			if (caloHit != -1 || nFound == NLAYERS)
			{
				if (nFound < NLAYERS) nFound = 0; //Only assume tracks in tracker with all hits
				if (ev.mFirstParticle + partFound >= mParticleBufferSize)
				{
					ev.mFirstParticle = 0;
					ev.mNParticles = 0;
					mError = 2;
					return(1);
				}
				particleResponse& part = mParticles[ev.mFirstParticle + partFound];
				part.mMCId = nCurMCParticle;

				if (nFound == 0 && hitsCalo[caloHit].mEnergy < 0) part.mPID = config.recPID[recoPID_Pion0];
				else if (nFound == 0)                             part.mPID = config.recPID[recoPID_Photon];
				else if (caloHit == -1)                           part.mPID = config.recPID[recoPID_Muon];
				else if (hitsCalo[caloHit].mEnergy < 0)           part.mPID = config.recPID[recoPID_Pion];
				else                                              part.mPID = config.recPID[recoPID_Electron];

				hitList& partHits = *(((hitList*) mTmpBuffer) + ev.mFirstParticle + partFound);
				for (int i = 0;i < NLAYERS;i++) partHits.hitId[i] = nFound ? recHits[i] : -1;
				partHits.hitId[NLAYERS] = caloHit;
				partFound++;

			}
		}
		nCurMCParticle++;
	}
	ev.mNParticles = partFound;
	ATOMIC_ADD(&nParticles, partFound);
	if (!kGPUCode) mNEntriesInBuffer += partFound; //Count real number of created response particles on CPU (serial processing)
	return(0);
}

template <class T> GPU() static inline void rotate(T& x, T& y, T& tgp, T& c0, T& c1, T dalpha)
{
	const T sinAlpha = std::sin(dalpha);
	const T cosAlpha = std::cos(dalpha);
	const T x0 = x, y0 = y;
	x = cosAlpha * x0 -sinAlpha * y0;
	y = sinAlpha * y0 + cosAlpha * y0;
	tgp += dalpha;
	c0 *= cosAlpha * cosAlpha;
	c1 *= cosAlpha;
}

template <class T> GPU() static inline void propagate(T& x, T toX, T& y, T& tgp, T& z, T& tgl, T* c)
{
	T deltaX = toX - x;
	x = toX;
	y += deltaX * tgp;
	z += deltaX * tgl;
	const T c1 = c[1];
	c[1] += deltaX * c[2];
	c[0] += deltaX * (c1 + c[1]);
	const T c4 = c[4];
	c[4] += deltaX * c[5];
	c[3] += deltaX * (c4 + c[4]);
}

template <class T> GPU() static inline void filter(T& y, T& tgp, T& z, T& tgl, T* c, T hitZ, T sigma2Y, T sigma2Z) //System is assumed to be rotated such that hitY = 0
{
	const T dY = -y;
	const T dZ = hitZ - z;
	const T Sy = c[0] + sigma2Y;
	const T Sz = c[3] + sigma2Z;
	const T Sy1 = 1.f / Sy;
	const T Sz1 = 1.f / Sz;
	y += dY * c[0] * Sy1;
	tgp += dY * c[1] * Sy1;
	z += dZ * c[3] * Sz1;
	tgl += dZ * c[4] * Sz1;
	c[2] = (c[2] - c[1] * (c[0] + sigma2Y)) * Sy1;
	c[1] = c[1] * sigma2Y * Sy1;
	c[0] = c[0] * sigma2Y * Sy1;
	c[5] = (c[5] - c[4] * (c[3] + sigma2Y)) * Sy1;
	c[4] = c[4] * sigma2Y * Sy1;
	c[3] = c[3] * sigma2Y * Sy1;
}

template <class T> GPU() static inline void greaterZero(T& val);
template <> GPU() inline void greaterZero<float>(float& val)
{
	if (val < 0) val = 0;
}
template <class S> GPU() inline float fetchHit(const S* hits, typename indexType<float>::type index, const float* memberPtr)
{
	index = index * (sizeof(S) / sizeof(float)) + ((const char*) memberPtr - (const char*) hits) / sizeof(float);
	return ((const float*) hits)[index];
}

#define FETCHHIT(mem, index, member) fetchHit(mem, index, &mem->member)

GPUi() void reconstruction::fitParticleCalo(particleResponse& part, hitList& partHits, const hitCalo* __restrict hitsCalo, float& lastR)
{
	float scale = FETCHHIT(hitsCalo, partHits.hitId[NLAYERS], mEnergy); //hadrons are flagged by negative sign
	scale = scale * scale - config.particles[part.mPID].mMass * config.particles[part.mPID].mMass;
	greaterZero(scale);
	float x[3];
	for (int i = 0;i < 3;i++) x[i] = FETCHHIT(hitsCalo, partHits.hitId[NLAYERS], mX[i]);
	scale = std::sqrt(scale / (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));
	for (int i = 0;i < 3;i++)
	{
		part.mX[i] = x[i];
		part.mP[i] = x[i] * scale;
	}
	part.mPt = std::sqrt(part.mP[0] * part.mP[0] + part.mP[1] * part.mP[1]);
	lastR = config.detectorRCalo;
}

GPUi() static float getMass(particleResponse part)
{
	return config.particles[part.mPID].mMass;
}

template <class T> GPUi() void reconstruction::fitParticle(particleResponse_T<T>& part, hitList_T<T>& partHits, const hitTracker* __restrict hitsTracker, const hitCalo* __restrict hitsCalo, T& lastR)
{
	T hitAlpha[NLAYERS + 1];
	for (int j = 0;j < NLAYERS;j++) hitAlpha[j] = std::atan2(FETCHHIT(hitsTracker, partHits.hitId[j], mX[1]), FETCHHIT(hitsTracker, partHits.hitId[j], mX[0]));
	hitAlpha[NLAYERS] = std::atan2(FETCHHIT(hitsCalo, partHits.hitId[NLAYERS], mX[1]), FETCHHIT(hitsCalo, partHits.hitId[NLAYERS], mX[0]));

	T x, y, tgp, z, tgl; //Track model
	T c[6]; //covariance

	//Kalman filter 1st step - Initialization (and rotation for 2nd step) - explicit formula
	//Assume start in detector system of hit 1 --> hitX = config.detector, hitY = 0
	T dalpha = hitAlpha[0] - hitAlpha[1];
	T sinAlpha = std::sin(dalpha);
	T cosAlpha = std::cos(dalpha);
	x = cosAlpha * config.detectorR[0]; //Before rotation, hit was at y = 0!
	y = sinAlpha * config.detectorR[0];
	float sigma2Y = config.resRPhi * config.resRPhi;
	float sigma2Z = config.resZ * config.resZ;
	c[0] = cosAlpha * cosAlpha * sigma2Y;
	//Rotated to detector system of hit 2 --> hitX = config.detectorR[1], hitY = 0

	//Kalman filter 2nd step, propagate and update with explicit formula
	T delta1 = 1.f / (config.detectorR[1] - x);
	tgp = -y * delta1;
	tgl = (FETCHHIT(hitsTracker, partHits.hitId[1], mX[2]) - FETCHHIT(hitsTracker, partHits.hitId[0], mX[2])) * delta1;
	z = FETCHHIT(hitsTracker, partHits.hitId[1], mX[2]);
	//sigma2Y = config.resRPhi * config.resRPhi; //sigma stays the same
	c[2] = (c[0] + sigma2Y) * delta1 * delta1;
	c[1] = sigma2Y * delta1;
	c[0] = sigma2Y;
	c[5] = 2.f * sigma2Z * delta1 * delta1;
	c[4] = sigma2Z * delta1;
	c[3] = sigma2Z;
	//Rotate to detector system of hit 3
	dalpha = hitAlpha[1] - hitAlpha[2];
	sinAlpha = std::sin(dalpha);
	cosAlpha = std::cos(dalpha);
	c[0] *= cosAlpha * cosAlpha;
	c[1] *= cosAlpha;
	x = cosAlpha * config.detectorR[1];
	y = sinAlpha * config.detectorR[1];
	tgp += dalpha;

	//Now all values are initialized and the infinities canceled, can start to use normal rotation / propagation / filter
	lastR = config.detectorR[2];
	for (int nLayer = 2;nLayer < NLAYERS;nLayer++)
	{
		if (nLayer != 2) rotate<T>(x, y, tgp, c[0], c[1], hitAlpha[nLayer - 1] - hitAlpha[nLayer]);
		propagate<T>(x, T(config.detectorR[nLayer]), y, tgp, z, tgl, c);
		filter<T>(y, tgp, z, tgl, c, FETCHHIT(hitsTracker, partHits.hitId[nLayer], mX[2]), sigma2Y, sigma2Z);
	}
	rotate<T>(x, y, tgp, c[0], c[1], hitAlpha[NLAYERS - 1] - hitAlpha[NLAYERS]);
	propagate<T>(x, T(config.detectorRCalo), y, tgp, z, tgl, c);
	filter<T>(y, tgp, z, tgl, c, FETCHHIT(hitsCalo, partHits.hitId[NLAYERS], mX[2]), config.resRPhiCalo * config.resRPhiCalo, config.resZCalo * config.resZCalo);
	lastR = config.detectorRCalo;

	//Rotate back in global system
	sinAlpha = std::sin(hitAlpha[NLAYERS]);
	cosAlpha = std::cos(hitAlpha[NLAYERS]);
	part.mX[0] = cosAlpha * x - sinAlpha * y;
	part.mX[1] = sinAlpha * x + cosAlpha * y;
	part.mX[2] = z;
	T scale = FETCHHIT(hitsCalo, partHits.hitId[NLAYERS], mEnergy); //hadrons are flagged by negative sign
	const T mass = getMass(part);
	scale = scale * scale - mass * mass;
	greaterZero(scale);
	scale /= 1 + tgp * tgp + tgl * tgl;
	scale = std::sqrt(scale);
	part.mP[0] = (cosAlpha - sinAlpha * tgp) * scale;
	part.mP[1] = (sinAlpha + cosAlpha * tgp) * scale;
	part.mP[2] = tgl * scale;
	part.mPt = std::sqrt(1 + tgp * tgp) * scale;
}

GPUi() void reconstruction::fit(const eventInfo& ev, const hitTracker* __restrict hitsTracker, const hitCalo* __restrict hitsCalo, const particleMC* __restrict particlesMC)
{
	for (unsigned int i = 0;i < ev.mNParticles;i++)
	{
		particleResponse& part = mParticles[ev.mFirstParticle + i];
		hitList& partHits = *(((hitList*) mTmpBuffer) + ev.mFirstParticle + i);
			float lastR;
		if (partHits.hitId[0] == -1) fitParticleCalo(part, partHits, hitsCalo, lastR);
		else fitParticle(part, partHits, hitsTracker, hitsCalo, lastR);
	}
}

#ifndef GPU_CODE
GPUh() void reconstruction::estimateParticleBufferSize(unsigned int nHitsTracker)
{
	unsigned int nParticles = nHitsTracker * 3; //Just some arbitrary factor 3 to make the buffer big
	if (nParticles % 16) nParticles += 16 - nParticles % 16;
	mParticleBufferSize = nParticles;
}

GPUh() void reconstruction::calculateScratchBufferSize(unsigned int nParticles, unsigned int nHitsTracker, unsigned int nHitsCalo)
{
	mTmpBufferSize = nParticles * (NLAYERS + 1) * sizeof(int);
}

void reconstruction::ReserveMemory(unsigned int nEvents, unsigned int nParticles, unsigned int nHitsTracker, unsigned int nHitsCalo)
{
	if (mParticleBufferSize > nParticles) nParticles = mParticleBufferSize;
	calculateScratchBufferSize(nParticles, nHitsTracker, nHitsCalo);
	if (nParticles > mParticlesMemory.size()) mParticlesMemory.resize(nParticles);
	if (nEvents > mEventInfosMemory.size()) mEventInfosMemory.resize(nEvents);
	mTmpMem.resize(mTmpBufferSize);
	mParticles = mParticlesMemory.data();
	mEventInfos = mEventInfosMemory.data();
	mTmpBuffer = mTmpMem.data();
	mParticleBufferSize = nParticles;
}

int reconstruction::reconstruct(const detector& __restrict det)
{
	estimateParticleBufferSize(det.mNEntriesInBuffer);
	ReserveMemory(det.mNEventCount, mParticleBufferSize, det.mNEntriesInBuffer, det.mNHitsCaloInBuffer);
	mNEventCount = det.mNEventCount;
	mParticlesMC = det.mParticlesMC;
	mEventInfosMC = det.mEventInfosMC;
	mNEntriesInBuffer = 0;
	if (gpu)
	{
		return(gpu->reconstructionGPU(*this, det));
	}
	else
	{
		for (unsigned int i = 0;i < det.mNEventCount;i++)
		{
			if (reconstruct(&det.mHitsTracker[det.mEventInfos[i].mFirstParticle], det.mEventInfos[i].mNParticles, &det.mHitsCalo[det.mEventInfosCalo[i].mFirstParticle], det.mEventInfosCalo[i].mNParticles, det.mEventInfos[i], det.mEventInfosCalo[i], mEventInfos[i])) return(1);
			fit(mEventInfos[i], &det.mHitsTracker[det.mEventInfos[i].mFirstParticle], &det.mHitsCalo[det.mEventInfosCalo[i].mFirstParticle], nullptr);
		}
	}
	return(0);
}
#endif
