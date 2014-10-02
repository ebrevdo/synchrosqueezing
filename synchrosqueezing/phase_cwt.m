% function w = phase_cwt(Wx, dWx, opt)
%
% Calculate the phase transform at each (scale,time) pair:
%   w(a,b) = Im( (1/2pi) * d/db [ Wx(a,b) ] / Wx(a,b) )
% Uses direct differentiation by calculating dWx/db in frequency
% domain (the secondary output of cwt_fw, see help cwt_fw)
%
% This is the analytic implementation of Eq. (7) of [1].
%
% 1. G. Thakur, E. Brevdo, N.-S. Fučkar, and H.-T. Wu,
% "The Synchrosqueezing algorithm for time-varying spectral analysis: robustness
%  properties and new paleoclimate applications," Signal Processing, 93:1079-1094, 2013.
%
% % 2. I. Daubechies, J. Lu, H.T. Wu, "Synchrosqueezed Wavelet Transforms: an
% empricial mode decomposition-like tool", Applied and Computational Harmonic Analysis
% 30(2):243-261, 2011.
%
% Input:
%  Wx: wavelet transform of x (see help cwt_fw)
%  dWx: samples of time derivative of wavelet transform of x (see help cwt_fw)
%  opt: options struct,
%    opt.gamma: wavelet threshold (default: sqrt(machine epsilon))
%
% Output:
%  w: phase transform, size(w) = size(Wx)
%
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Eugene Brevdo, Gaurav Thakur
%---------------------------------------------------------------------------------
function w = phase_cwt(Wx, dWx, opt)
    if nargin<3, opt = struct(); end
    % epsilon, gamma from [1]
    if ~isfield(opt, 'gamma'); opt.gamma = sqrt(eps); end

    % Calculate phase transform for each ai, normalize by (2*pi)
	if strcmpi(opt.dtype,'phase')
		u = unwrap(angle(Wx)).';
		w = [diff(u);u(end,:)-u(1,:)].'/(2*pi);
	else
		w = abs(imag(dWx./Wx/(2*pi)));
	end
    w(abs(Wx)<opt.gamma) = Inf;
end
