% function x = synsq_stft_fw(Tx, fs, opt, Cs, freqband)
%
% Inverse STFT Synchrosqueezing transform of Tx with associated
% frequencies in fs and curve bands in time-frequency plane
% specified by Cs and freqband.  This implements Eq. 5 of [1].
%
% 1. G. Thakur, E. Brevdo, N.-S. Fučkar, and H.-T. Wu,
% "The Synchrosqueezing algorithm for time-varying spectral analysis: robustness
%  properties and new paleoclimate applications," Signal Processing, 93:1079-1094, 2013.
%
% Input:
%   Tx, fs: See help synsq_cwt_fw
%   opt: options structure (see help synsq_cwt_fw)
%      opt.type: type of wavelet used in synsq_cwt_fw
%
%      other wavelet options (opt.mu, opt.s) should also match
%      those used in synsq_cwt_fw
%	Cs: (optional) curve centerpoints
%	freqs: (optional) curve bands
%
% Output:
%   x: components of reconstructed signal, and residual error
%
% Example:
%   [Tx,fs] = synsq_cwt_fw(t, x, 32); % Synchrosqueezing
%   Txf = synsq_filter_pass(Tx, fs, -Inf, 1); % Pass band filter
%   xf = synsq_cwt_iw(Txf, fs);  % Filtered signal reconstruction
%
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Eugene Brevdo, Gaurav Thakur
%---------------------------------------------------------------------------------
function x = synsq_stft_iw(Tx, fs, opt, Cs, freqband)
    if nargin<3, opt = struct(); end
	if nargin<4, Cs = ones(size(Tx,2),1); end
	if nargin<5, freqband = size(Tx,1); end
	
	%compute L2 norm of window to normalize inverse STFT with
	windowfunc = wfiltfn(opt.type,opt,false);
	C = quadgk(@(x) windowfunc(x).^2, -Inf, Inf);
	%quadgk is a bit inaccurate with the bump function, this scales it correctly
	if strcmpi(opt.type,'bump')
		C=C*0.8675;
	end
	
	% Invert Tx around curve masks in the time-frequency plane to recover
	% individual components; last one is the remaining signal
    % Integration over all frequencies recovers original signal
	% factor of 2 is because real parts contain half the energy
	x = zeros(size(Cs,1),size(Cs,2)+1);
	TxMask = zeros(size(Tx));
	TxRemainder = Tx;
	for n=[1:size(Cs,2)]
		TxMask = zeros(size(Tx));
		UpperCs=min(max(Cs(:,n)+freqband(:,n),1),length(fs));
		LowerCs=min(max(Cs(:,n)-freqband(:,n),1),length(fs));
		%Cs==0 corresponds to no curve at that time, so this removes such points from the inversion
		UpperCs(find(Cs(:,n)<1))=1;
		LowerCs(find(Cs(:,n)<1))=2;
		for m=[1:size(Tx,2)]
			TxMask(LowerCs(m):UpperCs(m),m) = Tx(LowerCs(m):UpperCs(m), m);
			TxRemainder(LowerCs(m):UpperCs(m),m) = 0;
		end
		x(:,n) = 1/(pi*C)*sum(real(TxMask),1).';
	end
	x(:,n+1) = 1/(pi*C)*sum(real(TxRemainder),1).';
	x = x.';
end
