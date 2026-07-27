#include "generator.h"

GPU() inline void generator::generatePrimParticle(particleMC& part, float disPowerMin, float disPower, unsigned short pid, float maxEta)
{
	float phi = 0.f;
	phi = disUniform(mGen);
	phi *= 2.f * MY_PI;
	float eta = 0.f;
	eta = disUniform(mGen);
	eta = eta * maxEta * 2.f - maxEta;
	part.mPt = genRandPower(disPowerMin, disPower);
	part.mPID = pid;
	part.mX[0] = part.mX[1] = part.mX[2] = 0.;
	part.mP[0] = std::cos(phi) * part.mPt;
	part.mP[1] = std::sin(phi) * part.mPt;
	part.mP[2] = std::sinh(eta) * part.mPt;
	part.mPrim = true;
}

GPUi() float generator::genRandPower(float onset, int power)
{
	float y = 0.f;
	do
	{
		y = disUniform(mGen);
	} while (y < 1e-20f); //Avoid FPE, cut at infinite momentum
	float ret = 1. / std::pow(y, 1.f / (power-1.f)) * onset;
	return(ret);
}

GPUi() static void getMass(float* mass, particleMC part)
{
	unsigned char pid = part.mPID;
	const unsigned short* decays = config.particles[pid].mDecay;
	for (int i = 0;i < 3;i++) mass[i] = config.particles[decays[i]].mMass;
	mass[3] = config.particles[pid].mMass;
}

template <class T> GPUi() int generator::decay3Body(particleMC_T<T> &part, vec3_T<T>* decayProducts)
{
	T randMomentum = 0.f, randY[2], randZ[2];
	randMomentum = disUniform(mGen);
	for (int k = 0;k < 2;k++)
	{
		randY[k] = disUniform(mGen);
		randZ[k] = disUniform(mGen);
	}


	T mass[4];
	getMass(mass, part);
	T totalEnergy = mass[3] - mass[0] - mass[1] - mass[2];
	T e[3] = {mass[0], randMomentum * totalEnergy + mass[0] + mass[1], mass[3]};

	T decayMomentum[3];
	for (int i = 0;i < 2;i++)
	{
		const T a = e[i+1], b = e[i], c = mass[2];
		const T x = (a - b - c) * (a + b + c) * (a - b + c) * (a + b - c);
		decayMomentum[i] = std::sqrt(x)/(2.f*a);
	}

	decayProducts[0].set(0.f, decayMomentum[0], 0.f);

	for (int i = 1;;i++)
	{
		decayProducts[i].set(0.f, -decayMomentum[i - 1], 0.f);

		T cosZ = 2.f * randZ[i - 1] - 1.f;
		T sinZ = std::sqrt(1 - cosZ * cosZ);
		T angY = 2.f * MY_PI * randY[i - 1];
		T cosY = std::cos(angY);
		T sinY = std::sin(angY);
		for (int j = 0;j <= i;j++)
		{
			T x = decayProducts[j].mX[0];
			T y = decayProducts[j].mX[1];
			T z = decayProducts[j].mX[2];;
			decayProducts[j].mX[1] = sinZ * x + cosZ * y;
			x = cosZ * x - sinZ * y;
			decayProducts[j].mX[0] = cosY * x - sinY * z;
			decayProducts[j].mX[2] = sinY * x + cosY * z;
		}

		if (i == 2) break;

		T beta = decayMomentum[i] / sqrt(decayMomentum[i] * decayMomentum[i] + e[i] * e[i]);
		beta = std::min(beta, 0.99999f);
		T b2 = beta * beta;
		T gamma = 1.f / std::sqrt(1.f - b2);
		for (int j = 0;j <= i;j++)
		{
			T bp = beta * decayProducts[j].mX[1];
			T gamma2 = (gamma - 1.0f) / b2;
			decayProducts[j].mX[1] += gamma2 * bp * beta + gamma * beta * std::sqrt(decayProducts[j].x2() + mass[j] * mass[j]);
		}
	}

	T p2 = part.p2();
	T gamma = std::sqrt(1.f + p2 / (mass[3] * mass[3]));
	T b2 = 1.f - 1.f / (gamma * gamma);
	T beta = std::sqrt(b2);

	T w = beta / std::sqrt(p2);
	T bx = w * part.mP[0], by = w * part.mP[1], bz = w * part.mP[2];
	for (int i = 0;i < 3;i++)
	{
		T e = std::sqrt(decayProducts[i].x2() + mass[i] * mass[i]);
		T bp = bx * decayProducts[i].mX[0] + by * decayProducts[i].mX[1] + bz * decayProducts[i].mX[2];
		T gamma2 = (gamma - 1.0f) / b2;

		decayProducts[i].mX[0] += gamma2 * bp * bx + gamma * bx * e;
		decayProducts[i].mX[1] += gamma2 * bp * by + gamma * by * e;
		decayProducts[i].mX[2] += gamma2 * bp * bz + gamma * bz * e;
	}

	return 0;
}

GPUi() void generator::insertSecondaries(particleMC& part, vec3* decayProducts, particleMC* nextParticle, unsigned int& nParticle)
{
	const float p2 = part.p2();
	const float gamma = std::sqrt(1.f + p2 / (config.particles[part.mPID].mMass * config.particles[part.mPID].mMass));
	const float beta = std::sqrt(1.f - 1.f / (gamma * gamma));

	float ctau = 0.f;
	ctau = disExp(mGen);
	ctau *= config.particles[part.mPID].mDecayLength;
	const float len1p = ctau * beta * gamma / std::sqrt(p2);

	for (int j = 0;j < 3;j++)
	{
		particleMC& partSec = nextParticle[nParticle++];
		partSec.mPID = config.particles[part.mPID].mDecay[j];
		for (int k = 0;k < 3;k++)
		{
			partSec.mX[k] = len1p * part.mP[k];
			partSec.mP[k] = decayProducts[j].mX[k];
		}
		partSec.mPt = std::sqrt(partSec.mP[0] * partSec.mP[0] + partSec.mP[1] * partSec.mP[1]);
		partSec.mPrim = false;
	}
}

GPUi() int generator::generateEvent(int evNum)
{
	unsigned int nCreate[config.nParticles];
	unsigned int myNParticles = 0;
	for (int i = 1;i < config.nParticles;i++)
	{
		nCreate[i] = disPoisson[i](mGen);
		myNParticles += nCreate[i] * (config.particles[i].mDecayLength > 0.f ? 4 : 1);
	}
	unsigned int myFirstParticle = mNParticles += myNParticles;
	if (myFirstParticle + myNParticles > mParticleBufferSize)
	{
		//We can increase the memory if we are on the host
		mParticleBufferSize *= 2;
		mParticlesMCMemory.resize(mParticleBufferSize);
		mParticlesMC = mParticlesMCMemory.data();
	}
	particleMC* nextParticle = &mParticlesMC[myFirstParticle];
	mEventInfosMC[evNum].mFirstParticle = myFirstParticle;
	mEventInfosMC[evNum].mNParticles = myNParticles;

	unsigned int countDecays = 0;
	for (int iType = 1;iType < config.nParticles;iType++)
	{
		unsigned int nParticle = 0;
		for (unsigned int i = 0;i < nCreate[iType];i++)
		{
			particleMC& part = nextParticle[nParticle++];
			generatePrimParticle(part, config.specPowerMin, config.particles[iType].mSpecPower, iType, config.particles[iType].mDecayLength > 0.f ? config.maxEtaDecay : config.maxEta);
		}

		if (config.particles[iType].mDecayLength > 0.f)
		{
			for (unsigned int i = 0;i < nCreate[iType];i++)
			{
				particleMC& part = nextParticle[i];
				vec3 decayProducts[3];
				decay3Body(part, decayProducts);
				insertSecondaries(part, decayProducts, nextParticle, nParticle);
			}
			countDecays += nCreate[iType];
		}
		nextParticle += nParticle;
	}
	nDecayParticles += countDecays;

	return(0);
}

generator::generator(std::RND_TYPE& gen) : mGen(gen), disPoisson(0), disExp(1.f), disUniform(0.f, 1.f)
{
	disPoisson.resize(config.nParticles);
	for (int i = 0;i < config.nParticles;i++)
	{
		disPoisson[i] = std::poisson_distribution<int>(config.particles[i].mdNdEtaMean * 2.f * (config.particles[i].mDecayLength <= 0.f ? config.maxEta : config.maxEtaDecay));
	}
}

GPUh() void generator::estimateBufferSize(int nEvents)
{
	mParticleBufferSize = 0;
	for (int i = 1;i < config.nParticles;i++)
	{
		mParticleBufferSize += config.particles[i].mdNdEtaMean * 2.f * (config.particles[i].mDecayLength <= 0.f ? config.maxEta : (4 * config.maxEtaDecay));
	}
	mParticleBufferSize *= 1.2 * (nEvents + 100);
}

void generator::ReserveMemory(unsigned int nEvents, unsigned int nParticles)
{
	if (nParticles < mParticleBufferSize) nParticles = mParticleBufferSize;
	if (nParticles > mParticlesMCMemory.size())
	{
		mParticlesMCMemory.resize(nParticles);
		mParticleBufferSize = nParticles;
	}
	if (nEvents > mEventInfosMemory.size()) mEventInfosMemory.resize(nEvents);
	mParticlesMC = mParticlesMCMemory.data();
	mEventInfosMC = mEventInfosMemory.data();

}

int generator::generateEvents(int nEvents)
{
	//Compute buffer size to allocate, reserve memory in std::vector, and set raw poiters which are used for processing
	estimateBufferSize(nEvents);
	ReserveMemory(nEvents, mParticleBufferSize);
	mNParticles = 0;
	mNEventCount = nEvents;
	for (int i = 0;i < nEvents;i++)
	{
		if (generateEvent(i)) return(1);
	}
	return(0);
}
