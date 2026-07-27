#ifndef GPUDEFINES_H
#define GPUDEFINES_H

#ifdef GPU_CODE
#define GPU() __device__
#define GPUi() __device__ inline
#define GPUh() __host__
#define GPUhi() __host__ inline
#define GPUhd() __host__ __device__
#else
#define GPU()
#define GPUi()
#define GPUh()
#define GPUhi()
#define GPUhd()
#endif

#endif
