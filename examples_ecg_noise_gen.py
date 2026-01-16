import ecg_noise_gen as ecgng 
import matplotlib.pyplot as plt
import numpy as np



"""
examples_ecg_noise_gen.py

Example usage of 'ecg_noise_gen':
    - Loads a raw ECG record from the MIT-BIH Arrhythmia Database.
    - Visualizes the clean ECG signal.
    - Generates different noise types.
    - Injects each noise source into the chosen ECG signal at a fixed SNR.
    - Displays the resulting noisy signals.

This file is illustrative of 'ecg_noise_gen' capabilities.

Author: David Velez
License: MIT
"""



plt.close('all')

# ECG Signal 
recording_id = '103'    # ECG record identifier
length = 30.0           # Signal duration [s]
fs = 360                # Sampling frequency [Hz]

# Choose a random start time for the segment
max_time = 120.0                                   # Example: max recording length (adjust to dataset)
t0 = np.random.randint(0, max_time - length + 1)  

# Load ECG Signal 
ecg_record = ecgng.get_record(recording_id, t0, t0 + length, fs)

plt.figure(figsize=(10,6))
plt.title(f"Raw ECG signal (Record {recording_id})")
plt.plot(ecg_record)
plt.show()

# Plot Each Noise Type
snr_db = 0
for noise_type in ecgng.noise_types:
    # Load Noise
    noise = ecgng.get_noise_sample(noise_type, t0, t0 + length, fs)
    # Add noise to ECG Signal at desired SNR 
    noisy_signal = ecgng.generate_noisy_signal(ecg_record, noise, snr_db)
    plt.figure(figsize=(10,6))
    plt.plot(noisy_signal)
    plt.title(
    f"Noisy ECG (Record {recording_id}, "
    f"Noise={noise_type}, SNR={snr_db} dB)")
    plt.show()

