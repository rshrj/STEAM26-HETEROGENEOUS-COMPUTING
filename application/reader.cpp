#include <stdio.h>
#include <memory>
#include <string.h>
#include <algorithm>
#include "reader.h"
#include "cmodules/qconfig.h"
#include "particles.h"

std::vector<std::string> configParticleNames;

int getParticlePID(const char* name)
{
	int retVal = 0;
	for (int j = 1;j < config.nParticles;j++)
	{
		if (configParticleNames[j].compare(name) == 0)
		{
			retVal = j;
			break;
		}
	}
	return retVal;
}

int readSimulationFile()
{
	printf("Reading simulation config file: %s\n", config.simulationFile);
	FILE* fp = fopen(config.simulationFile, "rb");
	if (fp == NULL)
	{
		printf("Error opening simulation config file\n");
		return 1;
	}

	fseek(fp, 0, SEEK_END);
	size_t fsize = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	std::unique_ptr<char[]> buffer_p(new char[fsize + 1]);
	char* buffer = buffer_p.get();
	if (fread(buffer, 1, fsize, fp) != fsize)
	{
		printf("Error reading file\n");
		fclose(fp);
		return 1;
	}
	fclose(fp);
	buffer[fsize] = 0;

	configParticleNames.clear();
	config.nParticles = 1;
	configParticleNames.emplace_back("Unknown");
	memset(&config.particles[0], 0, sizeof(config.particles[0]));
	std::vector<bool> particleUsed(1);

	int nLines = 0;
	char *ptr = buffer, *ptr2;
	constexpr int maxwords = 11;
	while (*ptr)
	{
		ptr2 = ptr;
		while (*ptr && *ptr != 10 && *ptr != 13) ptr++;
		while (*ptr == 10 || *ptr == 13) *(ptr++) = 0;
		nLines++;
		if (nLines == 1 || *ptr == '#') continue;

		char* words[maxwords];
		int nWords = 0;
		while (*ptr2)
		{
			while (*ptr2 == ' ' || *ptr2 == 9) ptr2++;
			if (*ptr2 == 0) break;
			words[nWords++] = ptr2;
			while (*ptr2 && *ptr2 != ' ' && *ptr2 != 9) ptr2++;
			if (*ptr2) *(ptr2++) = 0;
		}
		if (nWords == 0) continue;
		if (nWords != maxwords)
		{
			printf("Wrong number of fields (%d v.s. %d) in line %d\n", nWords, maxwords, nLines);
			return 1;
		}
		if (config.nParticles >= (config.runOnGPU ? 8 : 128))
		{
			printf("Maximum number of particles exceeded (%d)\n", config.nParticles);
			return 1;
		}
		configParticleNames.emplace_back(words[0]);
		auto& part = config.particles[config.nParticles];
		part.mMass = atof(words[1]);
		part.mDecayLength = atof(words[2]);
		if (part.mDecayLength <= 1e-6) part.mDecayLength = -1;
		part.mdNdEtaMean = atof(words[3]);
		particleUsed.push_back(part.mdNdEtaMean > 0.f);
		part.mSpecPower = atoi(words[4]);
		part.mHit = (atoi(words[5]) ? 1 : 0) + (atoi(words[6]) ? 1 : 0) * 2 + (atoi(words[7]) ? 1 : 0) * 4;
		if (part.mDecayLength <= 0.f)
		{
			part.mDecay[0] = part.mDecay[1] = part.mDecay[2] = -1;
		}
		else
		{
			for (int i = 0;i < 3;i++)
			{
				part.mDecay[i] = getParticlePID(words[8 + i]);
				if (part.mDecay[i] == 0)
				{
					printf("Invalid decay particle: %s\n", words[8 + i]);
					return 1;
				}
				particleUsed[part.mDecay[i]] = true;
			}
		}
		config.nParticles++;
	}

	if (config.nParticles <= 1)
	{
		printf("No particles defined\n");
		return 1;
	}

	int newCount = 1;
	for (int i = 1;i < config.nParticles;i++)
	{
		if (particleUsed[i])
		{
			if (i != newCount)
			{
				config.particles[newCount] = config.particles[i];
				configParticleNames[newCount] = configParticleNames[i];
				for (int j = 0;j < config.nParticles;j++)
				{
					for (int k = 0;k < 3;k++)
					{
						if (config.particles[j].mDecay[k] == i) config.particles[j].mDecay[k] = newCount;
					}
				}
			}
			newCount++;
		}
	}
	config.nParticles = newCount;

	for (int i = 1;i < config.nParticles;i++)
	{
		if (config.particles[i].mDecayLength > 0.f) std::sort(config.particles[i].mDecay, config.particles[i].mDecay + 3);
	}

	config.recPID[recoPID_Pion] = getParticlePID("Pion");
	config.recPID[recoPID_Pion0] = getParticlePID("Pion0");
	config.recPID[recoPID_Electron] = getParticlePID("Electron");
	config.recPID[recoPID_Muon] = getParticlePID("Muon");
	config.recPID[recoPID_Photon] = getParticlePID("Photon");

	printf("Read the following particles:\n");
	for (int i = 1;i < config.nParticles;i++)
	{
		printf("\t (%2d) %s: m = %f GeV, ct = %f cm, dNdEta = %f, specPower = %d, hitTracker = %d, hitCaloEM = %d, hitCaloHadron = %d", i, configParticleNames[i].c_str(), config.particles[i].mMass, config.particles[i].mDecayLength, config.particles[i].mdNdEtaMean, config.particles[i].mSpecPower,
			(int) ((config.particles[i].mHit & 1) > 0), (int) ((config.particles[i].mHit & 2) > 0), (int) ((config.particles[i].mHit & 4) > 0));
		if (config.particles[i].mDecayLength > 0.f)
		{
			printf(", decays into");
			for (int j = 0;j < 3;j++) printf(" %s", configParticleNames[config.particles[i].mDecay[j]].c_str());
		}
		printf("\n");
	}

	return 0;
}
