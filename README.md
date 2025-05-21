# Expectile-Periodogram

This paper introduces a novel periodogram-like function, called the expectile periodogram, for modeling spectral features of time series and detecting hidden periodicities. The expectile periodogram is constructed from trigonometric expectile regression, in which a specially designed check function is used to substitute the squared $l_2$ norm that leads to the ordinary periodogram. The expectile periodogram retains the key properties of the ordinary periodogram as a frequency-domain representation of serial dependence in time series, while offering a more comprehensive understanding by examining the data across the entire range of expectile levels. We  establish the asymptotic theory and investigate the relationship between the expectile periodogram and the so-called expectile spectrum. Simulations demonstrate the efficiency of the expectile periodogram in the presence of hidden periodicities. Finally, by leveraging the inherent two-dimensional nature of the expectile periodogram, we train a deep learning (DL) model to classify earthquake waveform data. Remarkably, our approach outperforms alternative periodogram-based methods in terms of classification accuracy.

One can download the code and data to reproduce the results in the following paper:
[1] Chen, T. (2024). Expectile Periodograms. arXiv preprint arXiv:2403.02060.

For further inqueries, please contact Tianbo Chen (chentianbo@ahu.edu.cn)
