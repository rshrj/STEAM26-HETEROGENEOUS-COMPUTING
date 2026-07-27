#ifndef QCONFIG_H
#define QCONFIG_H
#include "gpuDefines.h"

struct particleDef
{
	float mMass;
	float mDecayLength;
	float mdNdEtaMean;
	unsigned char mSpecPower;
	unsigned char mHit; //1 Tracker, 2 EM, 3 Hadron
	unsigned short mDecay[3];
	GPU() bool hitTracker() {return mHit & 1;}
	GPU() bool hitCalo() {return mHit & 6;}
	GPU() bool hitCaloEM() {return mHit & 2;}
	GPU() bool hitCaloHadron() {return mHit & 4;}
};

#include <vector>
#include <string>
extern std::vector<std::string> configParticleNames;
#endif

BeginConfig(structConfig, config)
AddHelpText("General settings")
AddOption(nEvents, unsigned long long int, 1000, "nEvents", 'n', "Number of events to generate", min(1ull))
AddOption(seed, int, -1, "seed", 0, "Random number seed (-1 for random)")
AddOption(quiet, int, 0, "quiet", 0, "Quiet mode, no console output (set to -1 to disable auto-quiet for MPI / ZMQ)")
AddHelp("help", 'h')

AddHelpText("Generator settings")
AddOption(simulationFile, const char*, "defaultSim.txt", "simulationFile", 's', "Text file with simulation settings")
AddOption(specPowerMin, float, 0.1f, "specPowerMin", 0, "Onset of pt spectra")

AddHelpText("Detector Geometry")
AddOptionArray(detectorR, float, 10, AddArrayDefaults(5.f,6.f,7.f,8.f,9.f,10.f,11.f,12.f,13.f,14.f), "detectorR", 'R', "Radius of detector (multiple layers", min(0.f))
AddOption(detectorRCalo, float, 15.f, "detectorRCalo", 0, "Radius of detector (calorimeter layer)", min(0.f))
AddVariable(detectorZ, float, 0.f)
AddOption(maxEta, float, 2.0f, "maxEta", 'e', "Maximum eta coverage of detector", min(0.f))
AddOption(maxEtaDecay, float, 8.0f, "maxEtaDecay", 'E', "Eta range for generated decaying particles", min(0.f))

AddHelpText("Detector Resolution")
AddOption(resPPhi, float, 1.f, "resPPhi", 0, "Phi Momentum resolution in degree")
AddOption(resPTheta, float, 1.f, "resPTheta", 0, "Theta Momentum resolution in degree")
AddOption(resZ, float, 0.1f, "resZ", 0, "Z position resolution in cm")
AddOption(resRPhi, float, 0.1f, "resRPhi", 0, "RPhi position resolution in cm")
AddOption(resZCalo, float, 1.f, "resZ", 0, "Z position resolution in cm (Calorimeter)")
AddOption(resRPhiCalo, float, 1.f, "resRPhi", 0, "RPhi position resolution in cm (Calorimeter)")
AddOption(resPt, float, 0.005, "resPt", 0, "Relative Pt resolution (sigmaPt / Pt*Pt), Calorimeter energy resolution in case of full simulation")

AddHelpText("Multithreading settings")
AddOption(nThreads, int, 1, "nThreads", 't', "Number of threads to run", min(1))

AddHelpText("GPU settings")
AddOption(runOnGPU, int, 0, "gpu", 'g', "Use the GPU for the generator (0: no, 1: yes, 2: for generation, 3: for generation and detectorm, 4 for generation, detector, and analysis", def(1))
AddOption(nGPUThreads, int, 1, "nGPUThreads", 0, "Number of OpenMP threads to use the GPU", min(0))
AddOption(nEventsOnGPU, unsigned int, 0, "nEventsOnGPU", 0, "Number of Events to compute on GPUs (subset of nEvents, in case nGPUThreads < nThreads)")
AddOption(gpuDevice, unsigned int, 0, "gpuDevice", 0, "GPU Device Number to use")
AddOption(gpuThreads, unsigned int, 256, "gpuThreads", 0, "Number of GPU threads per CUDA block")
AddOption(gpuBlocks, unsigned int, 40, "gpuBlocks", 0, "Number of GPU CUDA blocks")

AddVariable(nParticles, int, 0)
AddVariableAny(particles[128], particleDef)
AddVariableAny(recPID[6], unsigned short)

EndConfig()
