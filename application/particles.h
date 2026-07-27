#ifndef EVENTS_H
#define EVENTS_H

#include "gpuDefines.h"

#include <vector>
#include <cmath>

#if defined(USE_VC) && !defined(GPU_CODE)
#include "Vc/Vc"
using float_v  = Vc::float_v;
using int_v  = Vc::int_v;
#endif

//These particle classes have no constructor, such that creation of a large std::vector of them does not need to reinitilize the memory to 0 via constructors.
template <class T> struct vec3_T
{
	T mX[3];
	GPU() inline T x2() const {return(mX[0] * mX[0] + mX[1] * mX[1] + mX[2] * mX[2]);}
	GPU() inline T r() const {return(std::sqrt(mX[0] * mX[0] + mX[1] * mX[1]));}
	GPU() inline T x() const {return(std::sqrt(x2()));}
	GPU() void set(T x, T y, T z) {mX[0] = x; mX[1] = y; mX[2] = z;}
};

struct hitBase : vec3_T<float>
{
	int mMCId;
};

struct hitTracker : public hitBase
{
	int mLayer;
};

struct hitCalo : public hitBase
{
	float mEnergy;
};

template <class T> struct pidType;
template <> struct pidType<float> {typedef unsigned short type;};
#if !defined(GPU_CODE) && defined(USE_VC)
template <> struct pidType<float_v> {typedef int_v::IndexType type;};
#endif

template <class T> struct particleBase_T : public vec3_T<T>
{
	using vec3_T<T>::mX;
	T mP[3];
	T mPt;
	typename pidType<T>::type mPID;

	GPU() inline T p2() const {return(mPt * mPt + mP[2] * mP[2]);}
	GPU() inline T px() const {return(mX[0] * mP[0] + mX[1] * mP[1] + mX[2] * mP[2]);}
	GPU() inline T p() const {return(std::sqrt(p2()));}
	GPU() inline T eta() const
	{
		T cosTheta = mP[2] / p();
		return(-0.5 * std::log((1.0f - cosTheta) / (1.0f + cosTheta)));
	}

	GPU() T dca2() const
	{
		const T x2 = mX[0] * mX[0] + mX[1] * mX[1] + mX[2] * mX[2];
		const T xp = mX[0] * mP[0] + mX[1] * mP[1] + mX[2] * mP[2];
		const T p2 = mPt * mPt + mP[2] * mP[2];
		return(x2 - xp * xp / p2);
	}
};

template <class T> struct particleMC_T : public particleBase_T<T>
{
	using particleBase_T<T>::mX;
	using particleBase_T<T>::mP;
	using particleBase_T<T>::mPt;
	bool mPrim;

	GPU() int intersectWithRadius(T r, T* pos) const
	{
		if (mPrim)
		{
			T factor = r / mPt;
			for (int i = 0;i < 3;i++) pos[i] = mP[i] * factor;
			return(0);
		}
		T p2 = 2 * (mX[0] * mP[0] + mX[1] * mP[1]);
		T p3 = mX[0] * mX[0] + mX[1] * mX[1] - r * r;
		p2 /= 2 * mPt * mPt;
		p3 /= mPt * mPt;
		const T rad = p2 * p2 - p3;
		if (rad < 0) return(1);
		const T a = std::sqrt(rad) - p2;
		for (int k = 0;k < 3;k++) pos[k] = mX[k] + a * mP[k];
		return(0);
	}
};

template <class T> struct particleResponse_T : public particleBase_T<T>
{
	int mMCId;
};

typedef vec3_T<float> vec3;
typedef particleBase_T<float> particleBase;
typedef particleMC_T<float> particleMC;
typedef particleResponse_T<float> particleResponse;

#if defined(USE_VC) && !defined(GPU_CODE)
typedef vec3_T<float_v> vec3_v;
typedef particleBase_T<float_v> particleBase_v;
typedef particleMC_T<float_v> particleMC_v;
typedef particleResponse_T<float_v> particleResponse_v;
#endif

struct xParticleDebug : particleBase
{
	particleResponse mResponse[3];
	particleMC mMC[3];
	float mVtxChi2;
	float mDCA;
	float mInvMass;
	float mPt;
	float mCTau;
};

struct eventInfo
{
	unsigned int mFirstParticle;
	unsigned int mNParticles;
};

enum recoPID { recoPID_Pion = 0, recoPID_Pion0 = 1, recoPID_Electron = 2, recoPID_Muon = 3, recoPID_Photon = 4};

#endif
