I did this because of one code sample that I used in Python for the same, from scipy.signal.
import numpy as np
from scipy.signal import butter, filtfilt, iirnotch

# -----------------------------
# Parameters
# -----------------------------
Fs = 500.0  # Sampling frequency (Hz)
ecg = np.loadtxt("ecg_dataset.txt")

# -----------------------------
# Filter design
# -----------------------------
def bandpass(low, high, fs, order=4):
    nyq = fs / 2
    b, a = butter(order, [low/nyq, high/nyq], btype='band')
    return b, a

def lowpass(cutoff, fs, order=2):
    nyq = fs / 2
    b, a = butter(order, cutoff/nyq, btype='low')
    return b, a

def highpass(cutoff, fs, order=2):
    nyq = fs / 2
    b, a = butter(order, cutoff/nyq, btype='high')
    return b, a

# -----------------------------
# PQRST separation
# -----------------------------
b_pqrst, a_pqrst = bandpass(0.5, 40, Fs)
pqrst = filtfilt(b_pqrst, a_pqrst, ecg)

# -----------------------------
# Artefact separation
# -----------------------------
b_base, a_base = lowpass(0.5, Fs)
baseline = filtfilt(b_base, a_base, ecg)

b_emg, a_emg = highpass(40, Fs)
emg = filtfilt(b_emg, a_emg, ecg)

b_notch, a_notch = iirnotch(50, 30, Fs)
powerline = filtfilt(b_notch, a_notch, ecg)

The above code was intended to process an ECG, but the writer of the code was unaware that data could be in a .dat or a .csv format. 
Hence, in order to miniaturize this code and process it on an ESP or some other microcontroller, it would be highly convenient to have actual proper code instead of this high-level Python BS. 

So I wanted to build this library up from scratch, as I felt it would not only help my Capstone cause, but also make me experienced enough in building proper libraries.