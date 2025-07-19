# linear-filter-base-jl

# Signal Processing Tools in Julia

This repository contains a collection of Julia functions for preprocessing and analyzing electrophysiological signals (e.g., MEG/EEG). These tools cover a wide range of operations including filtering, analytic signal transformation, signal orthogonalization, source leakage correction, and more.

The toolbox is intended to support tasks in neuroscience and biomedical signal processing pipelines, particularly in time-frequency and source-space analysis of electrophysiological recordings.

---
## Overview

This package provides modular, flexible functions for:
- Computing and applying custom bandpass filters
- Extracting the analytic signal using the Hilbert transform
- Removing leakage and orthogonalizing signals
- Spectral decomposition and power extraction
- Concatenating or cutting epoched data
- Preprocessing signals from MFF and FIFF file formats

---

## Functions

### `msfun_prepare_cosine_filter.jl`
**Purpose:** Prepares a cosine taper filter in the frequency domain.  
**Inputs:** 
- `cfg`: A dictionary with filter type, window function, frequencies, and width.
- `N`: Signal length
- `sfreq`: Sampling frequency  
**Outputs:** 
- `win`: Time-domain window
- `F`: Frequency-domain filter coefficients

---

### `msfun_sig_filter.jl`
**Purpose:** Applies a cosine filter to continuous or epoched data.  
**Inputs:** 
- `X`: Signal (2D or 3D array)
- `sfreq`: Sampling frequency
- `cfg`: Filter configuration (same as `prepare_cosine_filter`)  
**Outputs:** 
- Filtered signal of same shape as input

---

### `msfun_sig_spectrum.jl`
**Purpose:** Computes spectral power for input signal.  
**Inputs:** 
- `X`: Signal (2D or 3D array)
- `sfreq`: Sampling frequency
- `cfg`: Struct specifying filtering parameters and output options  
**Outputs:** 
- `P`: Power matrix
- `cfg`: Updated configuration struct

---

### `msfun_sig_slow_modulation.jl`
**Purpose:** Removes fast oscillations to extract slow amplitude/phase modulations from narrowband analytic signal.  
**Inputs:** 
- `Z`: Analytic or real-valued signal (2D or 3D)
- `fcenter`: Frequency index to divide out  
**Outputs:** 
- `Zslow`: Slow-modulated signal

---

### `msfun_sig_downsample.jl`
**Purpose:** Downsamples a signal along its time dimension.  
**Inputs:** 
- `X`: Signal (2D or 3D array)
- `r`: Downsampling factor  
**Outputs:** 
- Downsampled signal

---

### `msfun_sig_concat_epoch.jl`
**Purpose:** Concatenates or epochs signal depending on input format and parameters.  
**Inputs:** 
- `sig`: Input data (2D or 3D)
- `L`: Number or length of epochs
- `type`: Either `"epochnum"` or `"epochlength"`  
**Outputs:** 
- Reshaped or concatenated signal

---

### `msfun_sig_inst_orth.jl`
**Purpose:** Removes instantaneous linear component of `X` predicted by `Y` (complex-valued signals).  
**Inputs:** 
- `X`: Target signal (complex)
- `Y`: Predictor signal (complex)  
**Outputs:** 
- Orthogonalized residual signal `Z`

---

### `msfun_sig_analytic.jl`
**Purpose:** Computes the analytic signal using the Hilbert transform along a given dimension.  
**Inputs:** 
- `X`: Real-valued signal
- `dim`: Dimension to transform (optional)  
**Outputs:** 
- Complex analytic signal

---

### `msfun_sig_leakcorr.jl`
**Purpose:** Removes signal leakage using various orthogonalization or regression methods.  
**Inputs:** 
- `X`: Signal to be corrected
- `Y`: Source of potential leakage
- `cfg`: Dict specifying method: `'gcs'`, `'orthinst'`, `'orthstat'`, or `'custom'`  
**Outputs:** 
- Corrected signal `Z`

---

### `sig_preprocess_fiff.jl`
**Purpose:** Reads and preprocesses signals from a FIFF file (MNE format).  
**Inputs:** 
- `raw`: FIFF raw data structure
- `times`: Array of time windows
- `cfg`: Configuration dict (includes channels, filtering, BLC, etc.)  
**Outputs:** 
- Preprocessed signal `sig`
- Updated configuration `cfg`

---

### `msfun_sig_preprocess_mff.jl`
**Purpose:** Reads and preprocesses signals from an MFF file (EGI format) using FieldTrip.  
**Inputs:** 
- `times`: Array of time windows
- `sfreq`: Sampling frequency
- `cfg`: Configuration dict (includes MFF path, filtering, BLC, etc.)  
**Outputs:** 
- Preprocessed signal `sig`
- Updated configuration `cfg`

---

## Usage

Each script is written as a standalone function in Julia and can be imported or included in your own data analysis workflows. For example:

```julia
include("msfun_sig_filter.jl")
filtered = msfun_sig_filter(data, sfreq, cfg)
```

## Dependencies
Julia 1.9+
FFTW.jl
For msfun_sig_preprocess_mff.jl: FieldTrip must be accessible via MAT.jl or equivalent wrapper
For msfun_sig_preprocess_fiff.jl: MNE-compatible .fif file reader if used
