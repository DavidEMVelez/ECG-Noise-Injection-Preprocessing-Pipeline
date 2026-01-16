import numpy as np
import wfdb 
from scipy.signal import detrend 

"""
ecg_noise_gen.py

This module provides tools aimed at biomedical data augmentation and robustness
evaluation, with focus on ECG signals.

The following functionalities are provided:
    - Loading ECG and noise records from the MIT-BIH Arrhythmia Database and
      the MIT-BIH Noise Stress Test Database.
    - Generation of synthetic noise sources (e.g., power-line interference
      and Gaussian noise).
    - Noise injection at a controlled signal-to-noise ratio (SNR) for
      reproducible experiments in biometric and biomedical signal processing.

Author: David Velez
License: 
"""

noise_types = "AWGN", "BW", "MA", "EMG", "PLI" # Available Noise Types

def get_noise_sample(noise_type:str, t0:int, tf:int, fs:int):
    '''
    Generates or loads a noise signal segment of a specified type.

    This function provides a unified interface for retrieving different
    noise sources.

        Loaded from MIT-BIH Noise Stress Test Database:
        - 'BW'  : Baseline Wander
        - 'MA'  : Motion Artifact
        - 'EMG' : Electromyographic noise
        
        Synthetically generated:
        - 'PLI' : Power-Line Interference (50 Hz sinusoid)
        - other : Additive White Gaussian Noise (AWGN)

    Returns noise segment of duration defined by [t0, tf] and
    sampling frequency fs.
    
    
    Parameters
    ----------
    noise_type : str
        Noise identifier.
    t0 : int
        Start time in seconds.
    tf : int
        End time in seconds.
    fs : int
        Sampling frequency in Hz.

    Returns
    -------
    noise : np.ndarray
        1-D array containing the generated or loaded noise segment.
    '''
    signal_length = int((tf-t0)*fs) # get signal lenght
    match noise_type:
        case 'BW':
            noise = get_record('bw',t0, tf, fs)
        case 'MA':
            noise = get_record('ma', t0, tf, fs)
        case 'EMG':
            noise = get_record('em', t0, tf, fs)
        case 'PLI':
            noise  = generate_50hz_noise(signal_length, fs) 
        case _:   # AWGN and Default
            noise = np.random.normal(loc=0, scale=1, size=signal_length) # Generate white Gaussian noise 
    return noise

def get_record(source:str, t0:int, tf:int, fs:int):
    '''    
    Load a record segment from a PhysioNet record.

    This function selects the appropriate database based on the record identifier:
        - Numeric identifiers → MIT-BIH Arrhythmia Database
        - Otherwise           → MIT-BIH Noise Stress Test Database
    
    The signal is cropped between t0 and tf (in seconds).

    Parameters
    ----------
    source : str
        Record identifier or filename.
    t0 : int
        Start time in seconds.
    tf : int
        End time in seconds.
    fs : int
        Sampling frequency in Hz.

    Returns
    -------
    ecg_signal : np.ndarray
        1-D array containing the signal segment from the chosen record.
    '''
    if str(source).isdigit():
        dir_path = 'mit-bih-arrhythmia-database-1.0.0/'
    else:
        dir_path = 'mit-bih-noise-stress-test-database-1.0.0/'
    record = wfdb.rdrecord(dir_path + source)
    signal = record.p_signal[ int(t0*fs) : int(tf*fs),1 ].squeeze()
    return signal

def generate_noisy_signal(signal:np.ndarray, noise:np.ndarray, snr_db:float):
    '''
    Add noise to a signal at a specified Signal-to-Noise Ratio (SNR).
    Input noise and signal lenghts are expect to match. Otherwise noise is 
    automatically resized to match the signal length.
    
    
    Parameters
    ----------
    signal : np.ndarray
        Input Signal.
    noise : np.ndarray
        Noise to be added.
    snr_db : float
        Desired signal-to-noise ratio in decibels (dB).

    Returns
    -------
    noisy_signal : np.ndarray
        Signal corrupted with scaled noise at the specified SNR.
    '''
    if len(noise) < len(signal):                         # Ensure noise is same length as signal
        repeats = int(np.ceil(len(signal) / len(noise))) # Repeat noise to match signal length
        noise = np.tile(noise, repeats)[:len(signal)]
    elif len(noise) > len(signal):
        noise = noise[:len(signal)]                      # Trim noise to same length
        
    noise = detrend(noise)          
    
    if np.mean(noise) == 0:
        return signal

    signal_power = np.mean(signal**2)            # Calculate the power of the signal // linear 
    snr_linear = 10**(snr_db / 10)               # Convert SNR in dB to linear ratio
    noise_power = signal_power / snr_linear      # Calculate the desired noise power 
    noise_std = np.sqrt(noise_power)             # Scale noise to achieve the desired noise power  
    noise = noise - np.mean(noise)               # Center the noise
    noise = noise * (noise_std / np.std(noise))  #
    noisy_signal = signal + noise                # Add the scaled Perlin noise to the signal
    return noisy_signal

def generate_50hz_noise(size:int, fs:int):
    '''
    This function generates  a sinusoidal signal at 50 Hz, commonly used to
    simulate electrical mains interference.
    

    Parameters
    ----------
    size : int
        Number of samples to generate..
    fs : int
        Sampling frequency in Hz..

    Returns
    -------
    pli_noise : np.ndarray
        1-D array containing synthetic PLI noise (50 Hz).
    '''
    t = np.arange(size) / fs                        # Time array
    pli_noise = 100 * np.sin(2 * np.pi * 50 * t)    # Generate 50Hz noise
    return pli_noise
