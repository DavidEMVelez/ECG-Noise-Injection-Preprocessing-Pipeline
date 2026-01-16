import wfdb
import os
import glob
import neurokit2 as nk
from scipy.signal import detrend
from scipy.signal import resample
import pandas as pd
import biosppy
import numpy as np
import ecg_noise_gen as ecgng
import random
import csv

'''
gen_dataset.py

This script is part of the experimental framework presented in the article:
"Understanding the Impact of Noise on ECG Biometrics: A Comparative Theoretical
and Experimental Analysis". The scritp generates noisy ECG datasets from clean
MIT-BIH Arrhythmia Database signals via the following operations:

1. Loads ECG recordings from clean MIT-BIH Arrhythmia Database
2. Generates various types of noise (as defined in `ecg_noise_gen`):
   - Power-line interference (PLI, 50Hz)
   - Baseline wander (BW)
   - Muscle artifacts (MA)
   - Electromyographic noise (EMG)
   - White Gaussian noise (AWGN)
3. Injects noise into each segment at multiple signal-to-noise ratio (SNR)
   levels.
4. Resamples ECG templates to a uniform length.
5. Saves processed noisy templates to CSV files.
6. Maintains a metadata CSV tracking:
   - File name
   - Database source
   - Noise type
   - SNR level

Parameters and settings such as `length`, `fs`, SNR range, and offset after R-peaks
can be configured at the beginning of the script.

Author: David Velez
License: 
'''

def detect_first_rpeaks(signal:np.ndarray, sampling_rate:int, n:int=20):
    """
    Detect the first N R-peaks in a ECG signal.
    
    Parameters
    ----------
    signal : np.ndarray
        ECG signal.
    sampling_rate : int
        Sampling frequency in Hz.
    n : int, optional
        Number of R-peaks to return (default=20)
    
    Returns
    -------
    r_peaks : np.ndarray
        Indices of the first N detected R-peaks.
    """
    ecg_cleaned = nk.ecg_clean(signal, sampling_rate=sampling_rate)
    _, info = nk.ecg_peaks(ecg_cleaned, sampling_rate=sampling_rate)
    return info["ECG_R_Peaks"][:n]

def get_signal_MIT(record_name:str, segment:int):
    """
    Load an ECG record from MIT-BIH database, split into thirds (for training,
    testing and validation), and extract a segment.
    
    Parameters
    ----------
    record_name : str
        Path to record
    segment : int
        Index of segment to select (0, 1, 2)
    
    Returns
    -------
    segment : np.ndarray
        Detrended ECG segment up to offset after last R-peak.
    """
    record = wfdb.rdrecord(record_name)     # Load record
    fs = record.fs
    signal = record.p_signal[:, 0]          # Lead 1 only
    n_samples = len(signal)                 # Split signal into segment
    slice3 = n_samples // 3
    halves = [signal[:slice3], signal[slice3:2*slice3],signal[2*slice3:]]
    signal = detrend(halves[segment])
    rpeaks = detect_first_rpeaks(signal, sampling_rate=fs, n=12)
    cutoff_idx = rpeaks[-1] + int((offset_ms / 1000) * fs)
    segment = signal[:cutoff_idx]
    return segment


def save_segments_to_csv(all_segments: np.ndarray, path:str):
    """
    Save a list of ECG segments to CSV with columns per segment.
    
    Parameters
    ----------
    all_segments : list of np.ndarray
        List containing ECG signal arrays.
    path : str
        Path to output CSV file.
    """
    # Convert list into dict with names seg_01, seg_02, ...
    seg_dict = {f"subj_{i+1:03d}": pd.Series(sig) for i, sig in enumerate(all_segments)}
    # Create DataFrame (pandas auto-pads shorter signals with NaN)
    df = pd.DataFrame(seg_dict)
    # Save to CSV
    df.to_csv(path, index=False)
    print(f"Saved {len(all_segments)} segments to {path}")
    
    
def clean_output_folder(path:str):
    """
    Delete all .csv files in the given directory.
    
    Parameters
    ----------
    path : str
        Folder path to clean.
    """
    files = glob.glob(os.path.join(out_dir, "*.csv")) # List all .csv files in the folder
    for f in files:
        if os.path.isfile(f):
            os.remove(f)  # delete .csv file
    
    
def process_and_save_segment(all_segments:np.ndarray, fs:int, target_len:int, 
                             append_String:str, out_dir:str, 
                             max_templates:int =10):
    """
    Process each segment through BioSPPy ECG pipeline, resample templates, and save to .csv.

    Parameters
    ----------
    all_segments : list of np.ndarray
        ECG segments to process
    fs : int
        Sampling rate in Hz
    target_len : int
        Number of samples per template after resampling
    append_String : str
        String to append to output CSV filename
    out_dir : str
        Output directory
    max_templates : int, optional
        Maximum number of templates per segment (default=10)
    """

    all_dfs = []
    for i, signal in enumerate(all_segments, start=1):  
        seg_name = f"{i:03d}"
        # Run BioSPPy ECG pipeline
        out = biosppy.signals.ecg.ecg(signal=signal, sampling_rate=fs, show=False)
        templates = out['templates']
        #if len(templates)<10:
        #    print("azar")
        # === Limit number of templates per segment ===
        templates = templates[:max_templates]
        # Resample each template to target_len samples
        resampled_templates = np.array([resample(beat, target_len) for beat in templates])
        # Build dataframe for this segment
        df = pd.DataFrame(resampled_templates, columns=[f"sample_{j+1}" for j in range(target_len)])
        df["seg_name"] = seg_name
        #df["template_id"] = range(1, len(df) + 1)
        all_dfs.append(df)
    # Concatenate everything
    big_df = pd.concat(all_dfs, ignore_index=True)
    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)
    # Save one single CSV
    filename = os.path.join(out_dir, f"{append_String}.csv")
    big_df.to_csv(filename, index=False)
    print(f"Saved {big_df.shape[0]} templates from {len(all_segments)} segments to {filename}")
    
# ----------------------------------------------------------------------------- 
# Script Main Execution
# ----------------------------------------------------------------------------- 

#----------------------------------------------     Pre-processing
# Set output directory 
out_dir="ecg_segments_csv"      
clean_output_folder(out_dir)    

# Define noise injection range and step
snr_db_min = -50
snr_db_max = 30
step = 5

# Keep X ms after last R-peak 
offset_ms = 500  

# Initialize segment lists for 3 segments of signal
all_segments_mitbih_adb_rec1 = []
all_segments_mitbih_adb_rec2 = []
all_segments_mitbih_adb_rec3 = []

# MIT sampling rate
fs = 360
target_len = 256

# Load MIT-BIH Arrhythmia Database header files
data_path_mitbih_adb= "mit-bih-arrhythmia-database-1.0.0"
hea_files_mitbih_adb = glob.glob(os.path.join(data_path_mitbih_adb, "**", "*.hea"), recursive=True)

# Extract first third of each record: Training
for hea_file in hea_files_mitbih_adb:
     record_name = os.path.splitext(hea_file)[0]
     segment = get_signal_MIT(record_name,0)
     all_segments_mitbih_adb_rec1.append(segment)

# Extract first third of each record: Testing
for hea_file in hea_files_mitbih_adb:
     record_name = os.path.splitext(hea_file)[0]
     segment = get_signal_MIT(record_name,1)
     all_segments_mitbih_adb_rec2.append(segment)

# Extract first third of each record: Validation
for hea_file in hea_files_mitbih_adb:
     record_name = os.path.splitext(hea_file)[0]
     segment = get_signal_MIT(record_name,2)
     all_segments_mitbih_adb_rec3.append(segment)


# Export clean ECG templates to .csv
process_and_save_segment(all_segments_mitbih_adb_rec1, fs, target_len, "mit_set_clean_training",     out_dir)
process_and_save_segment(all_segments_mitbih_adb_rec2, fs, target_len, "mit_set_clean_testing",      out_dir)
process_and_save_segment(all_segments_mitbih_adb_rec3, fs, target_len, "mit_set_clean_validation",   out_dir)


#----------------------------------------------     Noise Injection
length = 20;

csv_filename = os.path.join(out_dir, "file_list.csv")


# Create CSV file to list all noisy segments
with open(csv_filename, mode='w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["File Name", "Data Base", "Noise Type", "SNR (dB)"])

# Define SNR levels
snr_db_levels = list(range(snr_db_min, snr_db_max + 1, step))


# Define the three segment sets with their corresponding output prefix
segment_sets = [
    (all_segments_mitbih_adb_rec1, "training"),
    (all_segments_mitbih_adb_rec2, "testing"),
    (all_segments_mitbih_adb_rec3, "validation")
]

# Loop through all sets, noise types, and SNR levels
for segments, set_name in segment_sets:
    for noise_type in ecgng.noise_types:
        # Random integer start time t0
        t0 = random.randint(0, 1200)
        # Generate noise sample for this noise type
        noise = ecgng.get_noise_sample(noise_type, t0, t0 + length, fs)
        for snr_db in snr_db_levels:
            temp_segments = []
            # Inject noise into each segment
            for segment in segments:
                noisy_segment = ecgng.generate_noisy_signal(segment, noise, snr_db)
                temp_segments.append(noisy_segment)
            # Construct filename
            filename = f"mit_{set_name}_set_{noise_type}_SNR_{snr_db}_dB"
            # Process and save noisy segments
            process_and_save_segment(temp_segments, fs, 256, filename, out_dir)
            # Append metadata to CSV
            with open(csv_filename, mode='a', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([filename, "MIT-BIH-ADB", noise_type, snr_db])