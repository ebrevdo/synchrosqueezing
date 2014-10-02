% function [Tx,fs] = synsq_squeeze(Wx, w, t, nv, opt)
% function [Tx,fs] = synsq_squeeze(Sx, w, t, nv, opt)
%
% Calculates synchrosqueezed CWT or STFT of x. Used internally by synsq_cwt_fw and synsq_stft_fw.
%
% Input:
%   Wx or Sx: CWT or STFT of x
%   w: phase transform at same locations in T-F plane
%   t: time vector
%   nv: number of voices (CWT only; leave empty for STFT)
%   opt: options struct
%		transform: underlying time-frequency transform, 'CWT' or 'STFT'
%		freqscale: 'log' or 'linear' frequency bins/divisions
%		findbins: 'min' or 'direct' method to find bins, 'direct' is faster
%		squeezing: 'full' or 'measure', the latter corresponding to the approach in the paper [3],
%					which is not invertible but has better robustness properties in some cases
%
% Output:
%   Tx: synchrosqueezed output
%   fs: associated frequencies
%
% Note the multiplicative correction term x in synsq_cwt_squeeze_mex (and in
% the matlab equivalent code), required due to the fact that
% the squeezing integral of Eq. (2.7), in, [1], is taken w.r.t. dlog(a).
% This correction term needs to be included as a factor of Eq. (2.3), which
% we implement here.
%
% A more detailed explanation is available in Sec. III of [2].
% Note the constant multiplier log(2)/nv has been moved to the
% inverse of the normalization constant, as calculated in synsq_adm.m
%
% 1. I. Daubechies, J. Lu, H.T. Wu, "Synchrosqueezed Wavelet Transforms: an
% empricial mode decomposition-like tool"
%  Applied and Computational Harmonic Analysis, 30(2):243-261, 2011.
%
% 2. G. Thakur, E. Brevdo, N.-S. Fučkar, and H.-T. Wu,
% "The Synchrosqueezing algorithm for time-varying spectral analysis: robustness
%  properties and new paleoclimate applications," Signal Processing, 93:1079-1094, 2013.
%
% 3. G. Thakur and H.-T. Wu,
% "Synchrosqueezing-based Recovery of Instantaneous Frequency from Nonuniform Samples",
% SIAM Journal on Mathematical Analysis, 43(5):2078-2095, 2011.
%  
%
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Eugene Brevdo, Gaurav Thakur
%---------------------------------------------------------------------------------
function [Tx,fs] = synsq_squeeze(Wx, w, t, nv, opt)
    dt = t(2)-t(1);
    dT = t(end)-t(1);

	%default options
	if ~isfield(opt, 'freqscale') && strcmpi(opt.transform,'CWT'), opt.freqscale = 'log'; end
	if ~isfield(opt, 'freqscale') && strcmpi(opt.transform,'STFT'), opt.freqscale = 'linear'; end
	if ~isfield(opt, 'findbins'), opt.findbins = 'direct'; end
	if ~isfield(opt, 'squeezing'), opt.squeezing = 'full'; end
	
    % Maximum measurable (nyquist) frequency of data
    fM = 1/(2*dt);
    % Minimum measurable (fundamental) frequency
    fm = 1/dT;

	%na is number of scales for CWT, number of freqs for STFT
	[na, N] = size(Wx);
	% frequency divisions w_l to search over in Synchrosqueezing
	if strcmpi(opt.freqscale,'log')
		lfm = log2(fm); lfM = log2(fM);
%		fs = 2.^linspace(lfm, lfM, na);
%		fs = logspace(log10(fm), log10(fM), na);
		fs = [fm * (fM/fm).^([0:na-2]/(floor(na)-1)), fM];
	elseif strcmpi(opt.freqscale,'linear')
		if strcmpi(opt.transform,'CWT')
			fs = linspace(fm, fM, na);
		elseif strcmpi(opt.transform,'STFT')
			fs = linspace(0,1,N)/dt;
			fs = fs(1:floor(N/2));
		end
		dfs = 1/(fs(2)-fs(1));
	end
	
	if strcmpi(opt.transform,'CWT')
		as = 2^(1/nv) .^ [1:1:na]';
		scaleterm = as.^(-1/2);
	elseif strcmpi(opt.transform,'STFT')
		as = linspace(fm,fM,na);
		scaleterm = ones(size(as));
	end
	
	%measure version from reference [3]
	if strcmpi(opt.squeezing,'measure')
		Wx=ones(size(Wx))/size(Wx,1);
	end

	%incorporate threshold by zeroing out Inf values, so they get ignored below
	Wx(isinf(w)) = 0;

	Tx = zeros(length(fs),size(Wx,2));

	%now do squeezing by finding which frequency bin each phase transform point w(ai,b) lands in
	%look only at points where w(ai,b) is positive and finite
	if (strcmpi(opt.findbins,'direct') & strcmpi(opt.freqscale,'linear'))
		for b=1:N
		   for ai=1:length(as)
				% Find w_l nearest to w(a_i,b)
				k = min(max(round(w(ai,b)*dfs),1),length(fs));
				Tx(k, b) = Tx(k, b) + Wx(ai, b) * scaleterm(ai);
			end
		end 
	elseif (strcmpi(opt.findbins,'direct') & strcmpi(opt.freqscale,'log'))
		for b=1:N
		   for ai=1:length(as)
				% Find w_l nearest to w(a_i,b)
				% uses approximation w(a_i,b) ~= 2.^(lfm + k*(lfM-lfm)/na)
				k = min(max(1 + round(na/(lfM-lfm)*(log2(w(ai,b))-lfm)),1),na);
				Tx(k, b) = Tx(k, b) + Wx(ai, b) * scaleterm(ai);
			end
		end
	elseif (strcmpi(opt.findbins,'min'))
		for b=1:N
		   for ai=1:length(as)
				% Find w_l nearest to w(a_i,b)
				[V,k] = min(abs(w(ai,b)-fs));
				Tx(k, b) = Tx(k, b) + Wx(ai, b) * scaleterm(ai);
			end
		end
	end

	if strcmpi(opt.transform,'CWT')
		Tx = 1/nv * Tx;
	elseif strcmpi(opt.transform,'STFT')
		Tx = (fs(2)-fs(1)) * Tx;
	end

	% MEX version, deprecated (above code has been reworked to attain similar speed with JIT compiler)
	%Tx = 1/nv * synsq_cwt_squeeze_mex(Wx, w, as, fs, ones(size(fs)), lfm, lfM);
end
