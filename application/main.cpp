#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cmodules/timer.h"

#include "main.h"
#include "generator.h"
#include "detector.h"
#include "reconstruction.h"
#include "reader.h"
#include <fenv.h>

template <class T> static inline void setupCheckAvail(T& chk, bool avail, const char* text, T def = 0)
{
	if (chk != def && !avail)
	{
		printf("%s\n", text);
		chk = def;
	}
}

int setup(int argc, char** argv)
{
	int retVal = qConfigParse(argc, (const char**) argv);
	if (retVal)
	{
		if (retVal != qConfig::qcrHelp) printf("Error parsing command line arguments\n");
		return(1);
	}
	config.detectorZ = std::sinh(config.maxEta) * config.detectorRCalo;
	if (config.seed == -1)
	{
		std::random_device rd;
		config.seed = rd();
	}

	setupCheckAvail(config.runOnGPU, kGPUAvail, "Cannot run on GPU because CUDA support is disabled");
	setupCheckAvail(config.nGPUThreads, config.nGPUThreads <= config.nThreads, "Cannot have more GPU threads than threads", config.nThreads);

	if (config.nEventsOnGPU > config.nEvents) config.nEventsOnGPU = config.nEvents;

	if (readSimulationFile()) return 1;

	if (config.quiet <= 0)
	{
		printf("Particle toy simulation v0.1, coded by David Rohr in 2018\n");
		printf("Parameters:\n");
		qConfigPrint();
	}

	return(0);
}

void main_processing(int& nEventsTotal, int& nParticlesMCTotal, int& nDecayParticlesTotal, int& nParticlesTotal, int& nHitsCaloTotal, int& nRecParticlesTotal, double& timeGen, double& timeResponse, double& timeReconstruction, HighResTimer& timerTotal)
{
	HighResTimer timer;
	int iThread = 0;
	unsigned long long int nEvents = config.nEvents / config.nThreads;
	bool useGPU = config.runOnGPU && iThread < config.nGPUThreads;
	if (config.runOnGPU && config.nEventsOnGPU && config.nGPUThreads != config.nThreads)
	{
		nEvents = useGPU ? (config.nEventsOnGPU / config.nGPUThreads) : ((config.nEvents - config.nEventsOnGPU) / (config.nThreads - config.nGPUThreads));
	}

	std::random_device rd;
	randGen rand(config.seed + iThread);
	generator gen(rand.mGen);
	detector det(rand.mGen);
	reconstruction rec;
	gpuProcessor gpu;
	if (useGPU && gpu.gpuInit(gen, det, rec))
	{
		printf("Error initializing GPU\n");
		exit(1);
	}

	timer.Start();
	if (gen.generateEvents(nEvents)) {printf("Error during event generation\n"); exit(1);}
	timeGen += timer.GetCurrentElapsedTime();

	timer.ResetStart();
	if (det.generateResponse(gen)) {printf("Error during detector simulation\n"); exit(1);}
	timeResponse += timer.GetCurrentElapsedTime();

	timer.ResetStart();
	if (rec.reconstruct(det)) {printf("Error during reconstruction\n"); exit(1);}
	timeReconstruction += timer.GetCurrentElapsedTime();

	nParticlesMCTotal += gen.mNParticles;
	nDecayParticlesTotal += gen.nDecayParticles;
	nParticlesTotal += det.nParticles;
	nHitsCaloTotal += det.nHitsCalo;
	nRecParticlesTotal += rec.nParticles;

	nEventsTotal += nEvents;
}

int main(int argc, char** argv)
{
	//feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
	if (setup(argc, argv)) return(1);

	if (config.quiet <= 0) printf("\nProcessing events...\n");
	int nEventsTotal = 0, nParticlesMCTotal = 0, nParticlesTotal = 0, nHitsCaloTotal = 0, nRecParticlesTotal = 0, nDecayParticlesTotal = 0;
	double timeGen = 0., timeResponse = 0., timeReconstruction = 0.;
	HighResTimer timerTotal;
	timerTotal.Start();
	if (config.nThreads == 1)
	{
		main_processing(nEventsTotal, nParticlesMCTotal, nDecayParticlesTotal, nParticlesTotal, nHitsCaloTotal, nRecParticlesTotal, timeGen, timeResponse, timeReconstruction, timerTotal);
	}
	else
	{
		exit(1);
	}
	double timeTotal = timerTotal.GetCurrentElapsedTime();
	if (config.quiet <= 0)
	{
		printf("\nGenerated %d events, %d MC particles, %d decaying particles, %d hits in tracker, %d hits in calo, %d reconstructed particles\n", nEventsTotal, nParticlesMCTotal, nDecayParticlesTotal, nParticlesTotal, nHitsCaloTotal, nRecParticlesTotal);
		printf("Time: Total: %f seconds, (Generate events: %f seconds, Simulate detectors: %f seconds, Reconstruction: %f seconds)\n", timeTotal, timeGen / config.nThreads, timeResponse / config.nThreads, timeReconstruction / config.nThreads);
	}
}
