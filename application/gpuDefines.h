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

//Atomic function wrappers for the host
#ifdef GPU_CODE
	template <class T> GPUhd() T ATOMIC_ADD(T* var, T val) {return(atomicAdd(var, val));}
#else
	template <class T> T ATOMIC_ADD(T* var, T val) {T retVal = *var; *var += val;return(retVal);}
#endif

#endif
