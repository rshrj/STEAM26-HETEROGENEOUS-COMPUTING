#include "detector.h"
#include "reconstruction.h"

#ifndef GPU_CODE
#include <random>
#include <thread>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include "cmodules/timer.h"
#endif

static constexpr size_t correction_map_size = 500;

template <class T> struct indexType;
template <> struct indexType<float> {typedef int type;};

template <class T> struct hitList_T
{
	typename indexType<T>::type hitId[NLAYERS + 1];
};
typedef hitList_T<float> hitList;

GPUi() int reconstruction::reconstruct(const hitTracker* __restrict hitsTracker, int nHitsTracker, const hitCalo* __restrict hitsCalo, int nHitsCalo, const eventInfo& evTracker, const eventInfo& evCalo, eventInfo& ev)
{
	ev.mFirstParticle = evCalo.mFirstParticle; //Use calo hit id for output buffer index
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

#ifdef GPU_CODE
	atomicAdd(&nParticles, partFound);
#else
	nParticles.fetch_add(partFound);
#endif
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

GPUi() float reconstruction::getCorrection(float x_, float y_, float z, float alpha) // This is some bogus random correction, to blow up the memory size and add some more compute operations
{
	float sa = std::sin(alpha);
	float ca = std::cos(alpha);
	float x = ca * x_ - sa * y_;
	float y = sa * x_ + ca * y_;
	unsigned int xx = std::max<int>(0, std::min<int>(correction_map_size - 1, (x + config.detectorRCalo) * ((float)correction_map_size / (2.f * config.detectorRCalo))));
	unsigned int yy = std::max<int>(0, std::min<int>(correction_map_size - 1, (y + config.detectorRCalo) * ((float)correction_map_size / (2.f * config.detectorRCalo))));
	unsigned int zz = std::max<int>(0, std::min<int>(correction_map_size - 1, (z + config.detectorZ) * ((float)correction_map_size / (2.f * config.detectorZ))));
	float corr = mCorrectionMap[(xx * correction_map_size + yy) * correction_map_size + zz];
	return std::asinh(corr);
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
	c[0] += getCorrection(x, y, FETCHHIT(hitsTracker, partHits.hitId[0], mX[2]), hitAlpha[0]);
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
	c[0] += getCorrection(x, y, z, hitAlpha[1]);
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
		c[0] += getCorrection(x, y, z, hitAlpha[nLayer]);
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
	mCorrectionMemory.resize(CorrectionBufferSize());
	mTmpMem.resize(mTmpBufferSize);
	mParticles = mParticlesMemory.data();
	mEventInfos = mEventInfosMemory.data();
	mTmpBuffer = mTmpMem.data();
	mCorrectionMap = mCorrectionMemory.data();
	mParticleBufferSize = nParticles;
}

GPUh() size_t reconstruction::CorrectionBufferSize() const
{
	return correction_map_size * correction_map_size * correction_map_size;
}

int reconstruction::reconstruct(const detector& __restrict det)
{
	estimateParticleBufferSize(det.nParticles);
	ReserveMemory(det.mNEventCount, mParticleBufferSize, det.nParticles, det.nHitsCalo);

	std::mt19937 rng(config.seed);
	std::normal_distribution<float> dist(0.0f, 0.01);
	for (size_t i = 0; i < mCorrectionMemory.size(); i++) {
    	mCorrectionMemory[i] = dist(rng);
	}

	// Lambda function running the actual reconstruction, optionally multi-threaded oder with multi-processing
	auto reconstruct_lambda = [this, &det](int myId) {
		for (int repeat = 0; repeat < config.nBenchmarkRepeat; repeat++) {
			if (myId == 0) {
			printf("Running reconstruction (run %d)\n", repeat);
			}
			mNEventCount = det.mNEventCount;
			mParticlesMC = det.mParticlesMC;
			mEventInfosMC = det.mEventInfosMC;
			if (gpu)
			{
				if (gpu->reconstructionGPU(*this, det)) {
					exit(1);
				}
			}
			else
			{
				unsigned int eventsPerWorker = (det.mNEventCount + config.nThreads - 1) / config.nThreads;
				unsigned int myStart = myId * eventsPerWorker;
				unsigned int myEnd = std::min<unsigned int>((myId + 1) * eventsPerWorker, det.mNEventCount);
				for (unsigned int i = myStart;i < myEnd;i++)
				{
					if (reconstruct(&det.mHitsTracker[det.mEventInfos[i].mFirstParticle], det.mEventInfos[i].mNParticles, &det.mHitsCalo[det.mEventInfosCalo[i].mFirstParticle], det.mEventInfosCalo[i].mNParticles, det.mEventInfos[i], det.mEventInfosCalo[i], mEventInfos[i])) exit(1);
					fit(mEventInfos[i], &det.mHitsTracker[det.mEventInfos[300].mFirstParticle], &det.mHitsCalo[det.mEventInfosCalo[300].mFirstParticle], nullptr);
				}
			}
		}
	};

	HighResTimer timer;
	timer.Start();

	int processId = 0;
	std::vector<pid_t> processPIDs;
	if (config.nThreads && config.multiProcess) { // Start nThreads - 1 additional processes, if config.multiProcess set
		for (int i = 0; i < config.nThreads - 1; i++) {
			pid_t pid = fork();
			if (pid < 0) {
				printf("Error creating process!\n");
				exit(1);
			} else if (pid == 0) {
				processId = i + 1;
				break;
			} else {
				processPIDs.emplace_back(pid);
			}
		}
	}

	std::vector<std::thread> threadWorkers;
	if (config.nThreads > 1 && !config.multiProcess) { // Start nThreads - 1 additional threads, if config.multiProcess is not set
		threadWorkers.resize(config.nThreads - 1);
		for (int i = 0; i < config.nThreads - 1; i++) {
			threadWorkers[i] = std::thread(reconstruct_lambda, i + 1); // All newly started threads execute the reconstruct_lambda
		}
	}
	reconstruct_lambda(processId); // The main thread / process, and also all additional processes in case of config.multiProcess, execute the reconstruct_lambda
	if (config.nThreads > 1 && !config.multiProcess) { // Wait for threads to finish
		for (int i = 0; i < config.nThreads - 1; i++) {
			threadWorkers[i].join();
		}
		threadWorkers.clear();
	}
	if (config.nThreads > 1 && config.multiProcess) { // Wait for processes to finish
		for (int i = 0; i < config.nThreads - 1; i++) {
		    int status;
        	waitpid(processPIDs[i], &status, 0);
		}
		processPIDs.clear();
	}

	if (processId == 0) {
	  printf("Pure reconstruction time: %f seconds\n", timer.GetCurrentElapsedTime());
	} else {
	  exit(0);
	}

	return(0);
}
#endif
