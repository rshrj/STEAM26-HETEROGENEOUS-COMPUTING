# Heterogeneous Computing with CPUs and GPUs

CERN STEAM Academy 2026 — course repository.

**Conveners / speakers:** David Rohr (CERN)

David Rohr was born in Mannheim Germany, and is a physicist and developer of fast algorithms for event reconstruction at CERN. He works for the ALICE experiment and currently has the positions as lead developer for the TPC detector, leader of the working group for software for hardware accelerators, and computing coordinator of the ALICE 3 upgrade of the experiment.
David studied physics in Heidelberg and finished his PhD at the University of Frankfurt in 2015. He has always been working interdisciplinary in scientific computing. Since 2020, David holds a staff position at CERN.
Before coming to CERN, David was working for designing, benchmarking, and commissioning several computing centers in Germany and also the High Level Trigger farm of the ALICE experiment at CERN.
One of his focus is the low-level optimization of algorithms and applications, for both CPUs and GPUs, ranging from assembler code, explicit vectorization, and GPU computing to modern C++ code.

## Overview

The course will give an overview of heterogeneous computing for event reconstruction in high energy physics.
The topics will range from multi-threading v.s. multi-processing, over pipelines to hide latencies both on host and on GPU computing to dynamic workload distributions in heterogeneous environments.
In particular, it will show that different problems will require different computing solutions and that there is no magic bullet.
And it will point out the differences between classical high performance computing and computing in high energy physics.
The course will be complemented with some simple exercises.
The purpose of the exercise will not be lots of C++ developments, but rather to modify a provided reference application and perform some benchmarks.

## Environment setup

This course runs on the academy laptops (AlmaLinux 9). For the shared environment —
accounts, WiFi, editors, lxplus, AFS, CERNBox — see the
[STEAM Academy documentation](https://stac.docs.cern.ch/).

The course is supposed to be runnable in the default `lxplus-like:9` container at the NGT cluster, with the login as described here: [stac.docs.cern.ch/ngt/gpu-resources/](https://stac.docs.cern.ch/ngt/gpu-resources/)

Required compilers: GCC >= 11, CUDA >= 13, CMake

```sh
# make a compiler available (e.g. a CERN LCG build via environment modules)
module load <compiler>

# get the code
git clone https://github.com/CERN-STEAM-Academy/26-HETEROGENEOUS-COMPUTING.git
cd 26-HETEROGENEOUS-COMPUTING

# configure / build / run
```

## Materials

- Please find the slides at the indico timetable ([link](https://indico.cern.ch/event/1697464/timetable) below).
- The description of the exercises is at the end of the slides.
- The reference application of the exercises is in this repository.

## Schedule

See the [Indico timetable](https://indico.cern.ch/event/1697464/timetable).
