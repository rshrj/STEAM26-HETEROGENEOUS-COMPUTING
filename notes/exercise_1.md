Run 1 `./project --nEvents 100000 --nBenchmarkRepeat 3`
```
Reading simulation config file: defaultSim.txt
Read the following particles:
         ( 1) Electron: m = 0.000511 GeV, ct = -1.000000 cm, dNdEta = 2.500000, specPower = 8, hitTracker = 1, hitCaloEM = 1, hitCaloHadron = 0
         ( 2) X: m = 50.000000 GeV, ct = 0.100000 cm, dNdEta = 0.250000, specPower = 5, hitTracker = 1, hitCaloEM = 1, hitCaloHadron = 1, decays into Electron Electron Electron
Particle toy simulation v0.1, coded by David Rohr in 2018
Parameters:
    General settings:
        nEvents: 100000
        nBenchmarkRepeat: 3
        seed: -685402990
        quiet: 0
    Generator settings:
        simulationFile: defaultSim.txt
        specPowerMin: 0.1
    Detector Geometry:
        detectorR[10]: 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
        detectorRCalo: 15
        maxEta: 2
        maxEtaDecay: 8
    Detector Resolution:
        resPPhi: 1
        resPTheta: 1
        resZ: 0.1
        resRPhi: 0.1
        resZCalo: 1
        resRPhiCalo: 1
        resPt: 0.005
    Multithreading settings:
        nThreads: 1
        multiProcess: 0
    GPU settings:
        runOnGPU: 0
        gpuEventsAtOnce: 100000
        gpuDevice: 0
        gpuThreads: 256
        gpuBlocks: 40

Processing events...
Running reconstruction (run 0)
Running reconstruction (run 1)
Running reconstruction (run 2)
Pure reconstruction time: 5.556137 seconds

Generated 100000 events, 2603790 MC particles, 400794 decaying particles, 21443822 hits in tracker, 2109571 hits in calo, 6349644 reconstructed particles
Time: Generate events: 0.130735 seconds, Simulate detectors: 1.061920 seconds, Reconstruction: 6.885794 seconds
```


Run 2 `./project --nEvents 100000 --nBenchmarkRepeat 3 --nThreads 8`

```
Reading simulation config file: defaultSim.txt
Read the following particles:
         ( 1) Electron: m = 0.000511 GeV, ct = -1.000000 cm, dNdEta = 2.500000, specPower = 8, hitTracker = 1, hitCaloEM = 1, hitCaloHadron = 0
         ( 2) X: m = 50.000000 GeV, ct = 0.100000 cm, dNdEta = 0.250000, specPower = 5, hitTracker = 1, hitCaloEM = 1, hitCaloHadron = 1, decays into Electron Electron Electron
Particle toy simulation v0.1, coded by David Rohr in 2018
Parameters:
    General settings:
        nEvents: 100000
        nBenchmarkRepeat: 3
        seed: 1633332152
        quiet: 0
    Generator settings:
        simulationFile: defaultSim.txt
        specPowerMin: 0.1
    Detector Geometry:
        detectorR[10]: 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
        detectorRCalo: 15
        maxEta: 2
        maxEtaDecay: 8
    Detector Resolution:
        resPPhi: 1
        resPTheta: 1
        resZ: 0.1
        resRPhi: 0.1
        resZCalo: 1
        resRPhiCalo: 1
        resPt: 0.005
    Multithreading settings:
        nThreads: 8
        multiProcess: 0
    GPU settings:
        runOnGPU: 0
        gpuEventsAtOnce: 100000
        gpuDevice: 0
        gpuThreads: 256
        gpuBlocks: 40

Processing events...
Running reconstruction (run 0)
Running reconstruction (run 1)
Running reconstruction (run 2)
Pure reconstruction time: 0.963835 seconds

Generated 100000 events, 2597553 MC particles, 399366 decaying particles, 21407893 hits in tracker, 2106628 hits in calo, 6340098 reconstructed particles
Time: Generate events: 0.133712 seconds, Simulate detectors: 1.066962 seconds, Reconstruction: 2.442170 seconds
```


Run 3 `./project --nEvents 100000 --nBenchmarkRepeat 3 --nThreads 8 --multiProcess`

```
Reading simulation config file: defaultSim.txt
Read the following particles:
         ( 1) Electron: m = 0.000511 GeV, ct = -1.000000 cm, dNdEta = 2.500000, specPower = 8, hitTracker = 1, hitCaloEM = 1, hitCaloHadron = 0
         ( 2) X: m = 50.000000 GeV, ct = 0.100000 cm, dNdEta = 0.250000, specPower = 5, hitTracker = 1, hitCaloEM = 1, hitCaloHadron = 1, decays into Electron Electron Electron
Particle toy simulation v0.1, coded by David Rohr in 2018
Parameters:
    General settings:
        nEvents: 100000
        nBenchmarkRepeat: 3
        seed: 1134653993
        quiet: 0
    Generator settings:
        simulationFile: defaultSim.txt
        specPowerMin: 0.1
    Detector Geometry:
        detectorR[10]: 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
        detectorRCalo: 15
        maxEta: 2
        maxEtaDecay: 8
    Detector Resolution:
        resPPhi: 1
        resPTheta: 1
        resZ: 0.1
        resRPhi: 0.1
        resZCalo: 1
        resRPhiCalo: 1
        resPt: 0.005
    Multithreading settings:
        nThreads: 8
        multiProcess: 1
    GPU settings:
        runOnGPU: 0
        gpuEventsAtOnce: 100000
        gpuDevice: 0
        gpuThreads: 256
        gpuBlocks: 40

Processing events...
Running reconstruction (run 0)
Running reconstruction (run 1)
Running reconstruction (run 2)
Pure reconstruction time: 10.621118 seconds

Generated 100000 events, 2599033 MC particles, 399700 decaying particles, 21412638 hits in tracker, 2106789 hits in calo, 793020 reconstructed particles
Time: Generate events: 0.133074 seconds, Simulate detectors: 1.072255 seconds, Reconstruction: 13.320217 seconds
```


Run 4 `./project --nEvents 100000 --nBenchmarkRepeat 3 –g`

```
Reading simulation config file: defaultSim.txt
Read the following particles:
         ( 1) Electron: m = 0.000511 GeV, ct = -1.000000 cm, dNdEta = 2.500000, specPower = 8, hitTracker = 1, hitCaloEM = 1, hitCaloHadron = 0
         ( 2) X: m = 50.000000 GeV, ct = 0.100000 cm, dNdEta = 0.250000, specPower = 5, hitTracker = 1, hitCaloEM = 1, hitCaloHadron = 1, decays into Electron Electron Electron
Particle toy simulation v0.1, coded by David Rohr in 2018
Parameters:
    General settings:
        nEvents: 100000
        nBenchmarkRepeat: 3
        seed: 240975101
        quiet: 0
    Generator settings:
        simulationFile: defaultSim.txt
        specPowerMin: 0.1
    Detector Geometry:
        detectorR[10]: 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
        detectorRCalo: 15
        maxEta: 2
        maxEtaDecay: 8
    Detector Resolution:
        resPPhi: 1
        resPTheta: 1
        resZ: 0.1
        resRPhi: 0.1
        resZCalo: 1
        resRPhiCalo: 1
        resPt: 0.005
    Multithreading settings:
        nThreads: 1
        multiProcess: 0
    GPU settings:
        runOnGPU: 1
        gpuEventsAtOnce: 100000
        gpuDevice: 0
        gpuThreads: 256
        gpuBlocks: 40

Processing events...
Allocated 4699899136 bytes of GPU memory
Running reconstruction (run 0)
Running reconstruction (run 1)
Running reconstruction (run 2)
Pure reconstruction time: 0.456458 seconds

Generated 100000 events, 2599593 MC particles, 400038 decaying particles, 21417002 hits in tracker, 2106956 hits in calo, 6341502 reconstructed particles
Time: Generate events: 0.133460 seconds, Simulate detectors: 1.082083 seconds, Reconstruction: 1.915132 seconds
```