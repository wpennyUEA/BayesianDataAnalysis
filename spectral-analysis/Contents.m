% 
% SPECTRAL ESTIMATION
%
% Autoregressive (AR) and Multivariate Autoregressive (MAR) models 
% have parameters that can be transformed into the spectral domain.
% Spectra can be created at user-specified frequencies of interest
% (unlike FFTs which use integer multiples of fundamental frequencies).
% They support inferences about Granger Causality and 
% frequency-resolved measures of correlation and dependence.
%
% -------------------------------------------------------------------------
% spm_ar_freq.m         AR-based spectral estimation - see [1]
% demo_ar_spec.m        Demo
%
% spm_mar_spectra.m     MAR-based Spectral estimation
% demo_mar_spectra.m    Demo - see also [2]
%
% spm_granger.m         MAR-based Granger causality 
% demo_granger.m        Demo
%
% spm_wavspec.m         Wavelet based spectrogram
% demo_wavspec          Demo

% demo_gew.m            Demo of Spectrally Resolved Granger - see [3]
% demo_pdc.m            Demo of partial directed coherence (PDC)
%                       and Directed Transfer Function (DTF). See [4]
%
% spm_mmtspec.m         Moving multitape spectrogram
% spm_dpss.m            Use in spm_mmtspec.m
%
% -------------------------------------------------------------------------
% REFERENCES
%
% [1] J Pardey, S Roberts and L Tarassenko (1995). A review of parametric modelling
% techniques for EEG analysis. Med Eng. Phys, 18(1), 2-11
% 
% [2] Cassidy and Brown. Spectral Phase Estimates in the setting of 
% multidirectional coupling. J Neurosci Methods. 2003 Aug 1-15;37(3):299.
%
% [3] Brovelli et al. (2004) PNAS 101(26), 9849-9854 and Geweke (1982) 
% JASA 77 (378), 304-313.
%
% [4] L. Baccala and K. Sameshima (2001) Biol Cyb 84, 463-474.
%



