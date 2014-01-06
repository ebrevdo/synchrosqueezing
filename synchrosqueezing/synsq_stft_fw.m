% function [Tx, fs, Sx, Sfs, w] = synsq_stft_fw(t, x, opt)
%
% Calculates the STFT synchrosqueezing transform of vector x, with
% samples taken at times given in vector t. opt
% is a struct of options for the way synchrosqueezing is done, the
% type of wavelet used, and the wavelet parameters.  This
% implements the algorithm described in Sec. III of [1].
%
% 1. G. Thakur, E. Brevdo, N.-S. Fučkar, and H.-T. Wu,
% "The Synchrosqueezing algorithm for time-varying spectral analysis: robustness
%  properties and new paleoclimate applications," Signal Processing, 93:1079-1094, 2013.
%
% 2. G. Thakur and H.-T. Wu,
% "Synchrosqueezing-based Recovery of Instantaneous Frequency from Nonuniform Samples",
% SIAM Journal on Mathematical Analysis, 43(5):2078-2095, 2011.
%
% Input:
%  t: vector of times samples are taken (e.g. linspace(0,1,n))
%  x: vector of signal samples (e.g. x = cos(20*pi*t))
%  opt: struct of options
%    opt.type: type of wavelet (see help wfiltfn)
%      opt.s, opt.mu, etc (wavelet parameters: see opt from help wfiltfn)
%    opt.disp: display debug information?
%    opt.gamma: wavelet hard thresholding value (see help cwt_freq_direct)
%    opt.dtype: direct or numerical differentiation (default: 1,
%               uses direct).  if dtype=0, uses MEX differentiation
%               instead (see help cwt_freq), which is faster and
%               uses less memory, but may be less accurate.
%    
% Output:
%  Tx: synchrosqueezed output of x (columns associated with time t)
%  fs: frequencies associated with rows of Tx
%  Sx: STFT of x (see stft_fw)
%  Sfs: frequencies associated with rows of Sx
%  w: phase transform of Sx
%
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Eugene Brevdo, Gaurav Thakur
%---------------------------------------------------------------------------------
function [Tx, fs, Sx, Sfs, w, dSx] = synsq_stft_fw(t, x, opt)

if nargin<3, opt = struct(); end
if nargin<2, error('Too few input arguments'); end

% Choose some default values
if ~isfield(opt, 'type'), opt.type = 'gauss'; end
if ~isfield(opt, 'rpadded'), opt.rpadded = false; end

% Calculate sampling period, assuming regular spacing
dt = t(2)-t(1);

% Check for uniformity of spacing
if any(diff(t,2)/(t(end)-t(1))>1e-5)
    error('time vector t is not uniformly sampled');
end

% Calculate the modified STFT, using window of opt.winlen in frequency domain
x = x(:);
opt.stfttype = 'modified';
[Sx,Sfs,dSx] = stft_fw(x, dt, opt, t);

% Calculate phase transform
w = phase_stft(Sx, dSx, Sfs, opt, t);

% Calculate the synchrosqueezed frequency decomposition
% the parameter alpha from reference [2] is given by Sfs(2)-Sfs(1)
opt.transform = 'STFT';
[Tx,fs] = synsq_squeeze(Sx, w, t, [], opt);

end
