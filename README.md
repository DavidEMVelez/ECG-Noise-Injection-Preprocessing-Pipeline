# ECG Noise Injection Preprocessing Pipeline


This repository provides Python tools for generating noisy biosignal datasets with focus on ECG signals from the MIT-BIH Arrhythmia Database (MIT-BIH ADB).  This code aims to support biomedical data augmentation, biometrics algorithm robustness evaluation, and reproducible experiments in ECG biometrics. This work is part of the experimental framework described in the article: **"Understanding the Impact of Noise on ECG Biometrics: A Comparative Theoretical and Experimental Analysis"**  
Authors: David Velez

## Features
- Inject noise into clean biosignals signals at user-defined signal-to-noise ratio (SNR) levels.
- Load ECG records from MIT-BIH Arrhythmia Database
- Load noise records from MIT-BIH Noise Stress Test Database:
  - Baseline wander (BW)
  - Muscle artifacts (MA)
  - Electromyographic noise (EMG)
- Generate synthetic noise sources:
  - Power-line interference (PLI, 50Hz)
  - White Gaussian noise (AWGN)
- Extract ECG templates, resample to uniform length, and save to .CSV dataset.


