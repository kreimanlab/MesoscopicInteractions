# MesoscopicInteractions
Code and data release for "Mesoscopic functional interactions in the human brain reveal small-world properties", manuscript

<h2>Contents</h2>

1. /txt2h5eeg
    - Preprocessing, starting from NeuroWorks .txt export  
2. /synth
    - Preprocessing, data consolidation  
3. /v2
    - Artifact annotation  
4. /opencl
    - Coherence calculations  
5. /stats
    - Analysis, figure generation  

<h3>txt2h5eeg</h3>

Text to H5eeg conversion


Jiarui Wang :: jwang04@g.harvard.edu
Last edited: December 23, 2016


1   Introduction

    This package converts Neuroworks exported text to H5eeg .hdf5 format.


2   Requirements

    Matlab (never tested any verions earlier than 2016b)
    Python 3


4   Usage

    1   Configure where your exported data directories are located by editing
        the file "config"

    2   Any files within EXPORTED_DIR and its subdirectories starting with
        "Export", in the format "Export-suffix". The hyphen is required.

    3   Execute: python3 main.py <name of folder to convert>


5   Credits

    The H5eeg format specification was designed by the Nathan Crone Lab
