# Automated Dual-Comb Spectroscopy Pipeline

This repository contains a MATLAB-based Data Acquisition (DAQ) and Digital Signal Processing (DSP) pipeline for dual-comb spectroscopy. The code is designed to interface with a Siglent SSA3021X Spectrum Analyzer via TCP/IP to extract raw comb data, autonomously identify absorption features, and generate high-resolution datasets for RF-to-optical frequency mapping.

### Comb Profile Extraction
<img width="1897" height="1212" alt="comb_profile" src="https://github.com/user-attachments/assets/62dccc13-4a9e-428c-947d-36c5d39e3b5b" />

### Dynamic Dip Detection
<img width="1897" height="1096" alt="dip_detection" src="https://github.com/user-attachments/assets/134f46b9-88a9-429c-b473-8defc2b3013a" />

* **1. Comb Profile Extraction (`SSA_3021x_find_comb_profile.m`)**
  Executes a broad-span sweep (e.g., 20 MHz span, 100 Hz RBW) to calculate dynamic SNR thresholds, isolate the RF comb envelope, and determine the frequency boundaries required for the segmented survey.

* **2. Dynamic Dip Detection (`SSA_3021x_find_dips.m`)**
  Performs automated segmented surveying (e.g., 100 kHz span, 30 Hz RBW) across the isolated comb envelope. Utilizes a real-time moving-median DSP algorithm to track the local comb envelope and flag the center frequencies of absorption features.

* **3. Targeted High-Resolution Scanning (`SSA_3021x_dips_high_res.m`)**
  Executes narrow (e.g., 4 kHz span, 3 Hz RBW) sweeps strictly localized to the coordinates discovered in the previous step. It generates cleanly isolated datasets of the specific absorption features for curve-fitting and RF-to-optical frequency mapping.
  
* **4. File I/O Utility (`createUniqueFile.m`)**
  A custom data-handling function that injects real-time timestamps into output filenames to prevent accidental overwriting.

## Acknowledgements
Core instrument control framework and the `createUniqueFile` utility originated within the QSENSE group
 at the University of California, Los Angeles (UCLA).
