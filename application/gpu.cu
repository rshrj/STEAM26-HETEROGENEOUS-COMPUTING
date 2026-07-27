#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <array>

//Let everything that we include know, that we are compiling for GPU this time
#define GPU_CODE

//Define globally accessible constant config struct on GPU
#define QCONFIG_GPU
#include "cmodules/qconfig.h"
static __constant__ structConfig config;

#include "main.h"

//Make host config struct accessible as qConfigHost::config
namespace qConfigHost { extern "C" structConfig config; }
auto &configH = qConfigHost::config;

//Include CPP file with the actual implementation. These files will be compiled a second time in here to generate the GPU code
#include "generator.h"
#include "reconstruction.cpp"

struct gpuProcessor::internals
{
	std::array<cudaStream_t, 4> streams;
};

//Some debug output for CUDA errors
#define qon_mstr(a) #a
#define qon_mxstr(a) qon_mstr(a)
#define CUDACHK(cmd) if (cudaChk(cmd, qon_mxstr(cmd) "(" __FILE__ ": " qon_mxstr(__LINE__) ")")) return(1);
static inline int cudaChk(cudaError_t retVal, const char* msg)
{
	if (retVal != cudaSuccess)
	{
		printf("CUDA Error during command %s (Error code %d: %s)\n", msg, retVal, cudaGetErrorString(retVal));
		return(1);
	}
	return(0);
}

template <class T> static inline void gpuSetPointerWithSize(T* &ptr, char* &basePtr, size_t size, size_t alignment = 4096)
{
	size_t addr = (size_t) basePtr;
	if (addr % alignment) addr += alignment - addr % alignment;
	basePtr = (char*) addr;
	ptr = (T*) basePtr;
	basePtr += size;
}

template <class T> static inline void gpuSetPointer(T* &ptr, char* &basePtr, size_t count, size_t alignment = 4096)
{
	return gpuSetPointerWithSize(ptr, basePtr, count * sizeof(T), alignment);
}


void gpuProcessor::setGpuPointers(char* ptr, generator& gen, reconstruction& rec)
{
	char* ptrBase = ptr;
	//Detector
	gpuSetPointer(gpuDetGPU, ptr, 1);
	gpuSetPointer(gpuEventInfosDet, ptr, configH.nEvents);
	gpuSetPointer(gpuEventInfosDetCalo, ptr, configH.nEvents);
	gpuSetPointer(gpuHitsTracker, ptr, gen.mParticleBufferSize * NLAYERS);
	gpuSetPointer(gpuHitsCalo, ptr, gen.mParticleBufferSize);
	//Reconstruction
	gpuSetPointer(gpuRecGPU, ptr, 1);
	gpuSetPointer(gpuEventInfosRec, ptr, configH.nEvents);
	gpuSetPointer(gpuParticlesRec, ptr, rec.mParticleBufferSize);
	gpuSetPointer(gpuTmpRec, ptr, rec.mTmpBufferSize);
	gpuSetPointer(gpuCorrectionBuffer, ptr, rec.CorrectionBufferSize());
	gpuMemorySize = ptr - ptrBase;
}

gpuProcessor::gpuProcessor()
{
	mInternals = std::make_unique<internals>();
}

int gpuProcessor::gpuInit(generator& gen, detector& det, reconstruction& rec)
{
	int count;
	CUDACHK(cudaGetDeviceCount(&count));
	int device = configH.gpuDevice;
	if (device >= count)
	{
		printf("Invalid device requested %d, %d devices available\n", device, count);
		return 1;
	}
	CUDACHK(cudaSetDevice(device));

	gen.estimateBufferSize(configH.nEvents);
	rec.estimateParticleBufferSize(gen.mParticleBufferSize * NLAYERS);
	rec.calculateScratchBufferSize(gen.mParticleBufferSize, gen.mParticleBufferSize * NLAYERS, gen.mParticleBufferSize);
	if ((gpuDetCPU = (detector*) malloc(sizeof(detector))) == nullptr) return(1);
	if ((gpuRecCPU = (reconstruction*) malloc(sizeof(reconstruction))) == nullptr) return(1);
	setGpuPointers(nullptr, gen, rec);
	CUDACHK(cudaMalloc(&gpuMemoryPointerBase, gpuMemorySize));
	setGpuPointers(gpuMemoryPointerBase, gen, rec);
	printf("Allocated %'zu bytes of GPU memory\n", gpuMemorySize);

	CUDACHK(cudaMemcpyToSymbol(config, &configH, sizeof(structConfig), 0, cudaMemcpyHostToDevice));
	CUDACHK(cudaGetLastError());
	CUDACHK(cudaDeviceSynchronize());

	gpuDetCPU->mParticles = gpuParticles;
	gpuDetCPU->mHitsTracker = gpuHitsTracker;
	gpuDetCPU->mHitsCalo = gpuHitsCalo;
	gpuDetCPU->mEventInfos = gpuEventInfosDet;
	gpuDetCPU->mEventInfosCalo = gpuEventInfosDetCalo;

	gpuRecCPU->nParticles = 0;
	gpuRecCPU->mParticles = gpuParticlesRec;
	gpuRecCPU->mEventInfos = gpuEventInfosRec;
	gpuRecCPU->mTmpBuffer = gpuTmpRec;
	gpuRecCPU->mCorrectionMap = gpuCorrectionBuffer;
	gpuRecCPU->mParticleBufferSize = rec.mParticleBufferSize;

	CUDACHK(cudaMemcpy(gpuDetGPU, gpuDetCPU, sizeof(*gpuDetCPU), cudaMemcpyHostToDevice));
	CUDACHK(cudaMemcpy(gpuRecGPU, gpuRecCPU, sizeof(*gpuDetCPU), cudaMemcpyHostToDevice));

	for (unsigned int i = 0; i < mInternals->streams.size(); i++) {
		CUDACHK(cudaStreamCreate(&mInternals->streams[i]));
	}

	gen.gpu = det.gpu = rec.gpu = this;
	initialized = true;
	return(0);
}

int gpuProcessor::gpuExit()
{
	if (!initialized) {
		return 1;
	}
	for (unsigned int i = 0; i < mInternals->streams.size(); i++) {
		CUDACHK(cudaStreamDestroy(mInternals->streams[i]));
	}

	CUDACHK(cudaFree(gpuMemoryPointerBase));
	free(gpuDetCPU);
	free(gpuRecCPU);
	initialized = false;
	return(0);
}

gpuProcessor::~gpuProcessor()
{
	if (initialized && gpuExit())
	{
		printf("Error uninitializing GPU\n");
		exit(1);
	}
}

__global__ void reconstructionKernelGPU(reconstruction* __restrict rec, const detector* __restrict det, int start, int end)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	int size = blockDim.x * gridDim.x;
	for (int i = start + id;i < end;i += size)
	{
		if (rec->reconstruct(&det->mHitsTracker[det->mEventInfos[i].mFirstParticle], det->mEventInfos[i].mNParticles, &det->mHitsCalo[det->mEventInfosCalo[i].mFirstParticle], det->mEventInfosCalo[i].mNParticles, det->mEventInfos[i], det->mEventInfosCalo[i], rec->mEventInfos[i])) break;
		rec->fit(rec->mEventInfos[i], &det->mHitsTracker[det->mEventInfos[i].mFirstParticle], &det->mHitsCalo[det->mEventInfosCalo[i].mFirstParticle], nullptr);
	}
}

//Main host driven GPU processing functions - for running step by step
int gpuProcessor::reconstructionGPU(reconstruction& rec, const detector& det)
{
	CUDACHK(cudaMemcpy(gpuCorrectionBuffer, rec.mCorrectionMap, sizeof(*rec.mCorrectionMap) * rec.CorrectionBufferSize(), cudaMemcpyHostToDevice));

	gpuDetCPU->nParticles = det.nParticles;
	gpuDetCPU->nHitsCalo = det.nHitsCalo;
	gpuDetCPU->mNEventCount = det.mNEventCount;
	const auto& stream = mInternals->streams[0];
	CUDACHK(cudaMemcpyAsync(gpuDetGPU, gpuDetCPU, sizeof(*gpuDetCPU), cudaMemcpyHostToDevice, stream));
	// Run on GPU in batches of size config.gpuEventsAtOnce on CUDA stream stream
	for (int i = 0; i < det.mNEventCount; i += configH.gpuEventsAtOnce) {
		int end = std::min<int>(i + configH.gpuEventsAtOnce, det.mNEventCount);
		// Compute sizes and offsets for cudaMemcpy
		int iterationEventOffset = i;
		int iterationNEvents = end - i;
		int iterationParticleOffset = det.mEventInfos[i].mFirstParticle;
		int iterationNParticles = (end == det.mNEventCount ? det.nParticles : det.mEventInfos[i + configH.gpuEventsAtOnce].mFirstParticle) - iterationParticleOffset;
		int iterationCaloOffset = det.mEventInfosCalo[i].mFirstParticle;
		int iterationNCalo = (end == det.mNEventCount ? det.nHitsCalo : det.mEventInfosCalo[i + configH.gpuEventsAtOnce].mFirstParticle) - iterationCaloOffset;
		// Copy data to GPU
		CUDACHK(cudaMemcpyAsync(gpuEventInfosDet + iterationEventOffset, det.mEventInfos + iterationEventOffset, iterationNEvents * sizeof(eventInfo), cudaMemcpyHostToDevice, stream));
		CUDACHK(cudaMemcpyAsync(gpuEventInfosDetCalo + iterationEventOffset, det.mEventInfosCalo + iterationEventOffset, iterationNEvents * sizeof(eventInfo), cudaMemcpyHostToDevice, stream));
		CUDACHK(cudaMemcpyAsync(gpuHitsTracker + iterationParticleOffset, det.mHitsTracker + iterationParticleOffset, iterationNParticles * sizeof(hitTracker), cudaMemcpyHostToDevice, stream));
		CUDACHK(cudaMemcpyAsync(gpuHitsCalo + iterationCaloOffset, det.mHitsCalo + iterationCaloOffset, iterationNCalo * sizeof(hitCalo), cudaMemcpyHostToDevice, stream));
		// Execute GPU kernel
		reconstructionKernelGPU<<<configH.gpuBlocks, configH.gpuThreads, 0, stream>>>(gpuRecGPU, gpuDetGPU, i, end);
		// Copy data back to host
		CUDACHK(cudaGetLastError());
		CUDACHK(cudaMemcpyAsync(rec.mEventInfos + iterationEventOffset, gpuEventInfosRec + iterationEventOffset, iterationNEvents * sizeof(eventInfo), cudaMemcpyDeviceToHost, stream));
		CUDACHK(cudaMemcpyAsync(rec.mParticles + iterationCaloOffset, gpuParticlesRec + iterationCaloOffset, iterationNCalo * sizeof(particleResponse), cudaMemcpyDeviceToHost, stream));
	}
	CUDACHK(cudaMemcpyAsync(gpuRecCPU, gpuRecGPU, sizeof(*gpuRecCPU), cudaMemcpyDeviceToHost, stream));
	CUDACHK(cudaStreamSynchronize(stream));
	rec.nParticles = gpuRecCPU->nParticles;

	return(0);
}
