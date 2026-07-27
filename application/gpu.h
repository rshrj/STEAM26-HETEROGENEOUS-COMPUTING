#ifndef GPU_H
#define GPU_H

#include "particles.h"

class generator;
class detector;
class reconstruction;
class simpleHist;
#ifndef GPU_CODE
struct curandState_t;
#endif

class gpuProcessor
{
public:
	gpuProcessor();
	~gpuProcessor();
	int gpuInit(generator& gen, detector& det, reconstruction& rec);
	int generatorGPU(generator& gen);
	int detectorGPU(detector& det, const generator& gen);
	int reconstructionGPU(reconstruction& rec, const detector& det);
#ifdef USE_GPU
private:
	void setGpuPointers(char* ptr, generator& gen, reconstruction& rec);
	int gpuExit();

	detector* gpuDetCPU = nullptr;
	detector* gpuDetGPU = nullptr;
	eventInfo* gpuEventInfosDet = nullptr;
	eventInfo* gpuEventInfosDetCalo = nullptr;
	particleResponse* gpuParticles = nullptr;
	hitTracker* gpuHitsTracker = nullptr;
	hitCalo* gpuHitsCalo = nullptr;

	reconstruction* gpuRecCPU = nullptr;
	reconstruction* gpuRecGPU = nullptr;
	eventInfo* gpuEventInfosRec = nullptr;
	particleResponse* gpuParticlesRec = nullptr;
	char* gpuTmpRec = nullptr;

	char* gpuMemoryPointerBase = nullptr;
	size_t gpuMemorySize = 0;
	bool initialized = false;
#endif
};

//GPU Processing function definition
#ifndef USE_GPU
	inline gpuProcessor::gpuProcessor() {}
	inline gpuProcessor::~gpuProcessor() {}
	inline int gpuProcessor::gpuInit(generator& gen, detector& det, reconstruction& rec) {return 1;}
	inline int gpuProcessor::generatorGPU(generator& gen) {return 1;}
	inline int gpuProcessor::detectorGPU(detector& det, const generator& gen) {return 1;}
	inline int gpuProcessor::reconstructionGPU(reconstruction& rec, const detector& det) {return 1;}
#endif

#endif
