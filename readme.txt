Synchrosqueezing Toolbox v1.21
Authors: Eugene Brevdo, Gaurav Thakur

If you use this code in your research, please include the citations to the 
papers below.

Acknowledgements
------------------

Hau-tieng Wu and Ingrid Daubechies were instrumental in the creation
of the underlying ideas behind the Synchrosqueezing transform code in
this toolbox.

The file curve_ext.c was authored by Jianfeng Lu (now at Duke).


Introduction
--------------

This toolbox implements several time-frequency and time-scale analysis
methods, as described in [1,2,3,4].  These include:

1. Forward and inverse discretized Continuous Wavelet Transform (CWT).
2. Forward and inverse CWT-based Synchrosqueezing (synsq_cwt)
   and STFT-based Synchrosqueezing (synsq_stft).
3. Various plotting and curve extraction utilities to go with
   the above functions.
4. A GUI (the Synchrosqueezing GUI) for analysis, filtering,
   denoising, and signal extraction.  This allows for simple, fast,
   initial analysis via Synchrosqueezing.

Installation
-------------

0. You will need the MATLAB Image Processing toolbox for some of the
   GUI filtering tools used by gsynsq (see basic command reference
   below).

1. Put the contents of the synchrosqueezing toolbox somewhere (say,
   $HOME/matlab/).

2. Add the new directories to your path permanently; e.g., add the
   following to your startup.m:
     addpath ~/matlab/synchrosqueezing;
     addpath ~/matlab/util;

Basic command reference
------------------------

gsynsq: The synchrosqueezing GUI.  Run it without any parameters and
use the menu items and keyboard shortcuts to analyze signals.  To
import time/signal vectors from your workspace, use file->import
(control+I).  To export analysis output when done, use file->export
(control+E).  This GUI calls, and provides examples for using, many of the
functions described below.

cwt_fw, cwt_iw: continuous wavelet forward/inverse transform

synsq_cwt_fw, synsq_cwt_iw: CWT Synchrosqueezing
  forward/inverse transform.

synsq_stft_fw, synsq_stft_iw: STFT Synchrosqueezing
  forward/inverse transform.

est_riskshrink_thresh: RiskShrink threshold estimator
  for denoising signals.

tplot, tplot_power: plotting of output of cwt_fw, stft_fw, synsq_cwt, synsq_stft

curve_ext_*: different types of curve extraction from output of synsq_cwt or synsq_stft
plot_ext_curves*: associated plotting functions

synsq_filter_pass: frequency-region filtering in synchrosqueezing
domain.

examples: several example scripts illustrating how the above functions work

References
-------------
1. Mallat, S., Wavelet Tour of Signal Processing 3rd ed.

2. I. Daubechies, J. Lu, H.T. Wu, "Synchrosqueezed Wavelet Transforms: an
   empricial mode decomposition-like tool", Applied and Computational Harmonic Analysis
 30(2):243-261, 2011.

3. G. Thakur, E. Brevdo, N.-S. Fučkar, and H.-T. Wu,
 "The Synchrosqueezing algorithm for time-varying spectral analysis: robustness
  properties and new paleoclimate applications", Signal Processing 93:1079-1094, 2013.

4. G. Thakur and H.-T. Wu,
 "Synchrosqueezing-based Recovery of Instantaneous Frequency from Nonuniform Samples",
 SIAM Journal on Mathematical Analysis, 43(5):2078-2095, 2011.