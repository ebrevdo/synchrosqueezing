% function w = phase_stft(Sx, dSx, fs, opt)
%
% Calculate the phase transform of modified STFT at each (freq,time) pair:
%   w(a,b) = Im( eta - d/dt [ Sx(t,eta) ] / Sx(t,eta) / (2*pi*i) )
% Uses direct differentiation by calculating dSx/dt in frequency
% domain (the secondary output of stft_fw, see help stft_fw)
%
% 1. G. Thakur and H.-T. Wu,
% "Synchrosqueezing-based Recovery of Instantaneous Frequency from Nonuniform Samples",
% SIAM Journal on Mathematical Analysis, 43(5):2078-2095, 2011.
%
% 2. G. Thakur, E. Brevdo, N.-S. Fučkar, and H.-T. Wu,
% "The Synchrosqueezing algorithm for time-varying spectral analysis: robustness
%  properties and new paleoclimate applications," Signal Processing, 93:1079-1094, 2013.
%
%
% Input:
%  Sx: wavelet transform of x (see help stft_fw)
%  dSx: samples of time derivative of STFT of x (see help stft_fw)
%  opt: options struct,
%    opt.gamma: wavelet threshold (default: sqrt(machine epsilon))
%
% Output:
%  w: phase transform, size(w) = size(Sx)
%
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Eugene Brevdo, Gaurav Thakur
%---------------------------------------------------------------------------------
function w = phase_stft(Sx, dSx, Sfs, opt, t)
    if nargin<3, opt = struct(); end
    % epsilon, gamma from [1] and [2]
    if ~isfield(opt, 'gamma'); opt.gamma = sqrt(eps); end

    % Calculate phase transform; modified STFT amounts to an extra frequency term
	w = repmat(Sfs.',[1,length(t)]) - imag(dSx./Sx/(2*pi));
	
%	threshold out small points
	w(abs(Sx)<opt.gamma) = Inf;
end
