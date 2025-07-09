# Expectile-Periodogram

This paper introduces a novel periodogram-like function, called the expectile periodogram, for modeling spectral features of time series and detecting hidden periodicities. The expectile periodogram is constructed from trigonometric expectile regression, in which a specially designed check function is used to substitute the squared $l_2$ norm that leads to the ordinary periodogram. The expectile periodogram retains the key properties of the ordinary periodogram as a frequency-domain representation of serial dependence in time series, while offering a more comprehensive understanding by examining the data across the entire range of expectile levels. We  establish the asymptotic theory and investigate the relationship between the expectile periodogram and the so-called expectile spectrum. Simulations demonstrate the efficiency of the expectile periodogram in the presence of hidden periodicities. Finally, by leveraging the inherent two-dimensional nature of the expectile periodogram, we train a deep learning (DL) model to classify earthquake waveform data. Remarkably, our approach outperforms alternative periodogram-based methods in terms of classification accuracy.

One can download the code and data to reproduce the results in the paper uploaded.

Code:
1.  "ep.r" contains functions to compute the expectile periodogram. 
2.  "simu.r" reproduces the results in the simulations. 
3.  "ep.py" reproduces the results in the earthquake data classification.

Data:
1.  "Patient-3_seizure-2_Q5-Q6_17-21.mat", EEG data.
2.  "sp.csv", S&P500 data.
3.  "ep500.csv", earthquake data.   (download the data from https://drive.google.com/file/d/1u1gHarQinYWFAvXO0f5radzSEOGuhTQ9/view?usp=drive_link before running)

For further inquiries, please contact Tianbo Chen (chentianbo@ahu.edu.cn).
