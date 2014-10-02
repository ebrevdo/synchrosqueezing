% function [Sx,Sfs,dSx] = stft_fw(x, dt, opt)
%
% Compute the short-time Fourier transform and modified short-time Fourier transform from [1].
% The former is very closely based on Steven Schimmel's stft.m and istft.m from his
% SPHSC 503: Speech Signal Processing course at Univ. Washington.
%
% 1. G. Thakur and H.-T. Wu,
% "Synchrosqueezing-based Recovery of Instantaneous Frequency from Nonuniform Samples",
% SIAM Journal on Mathematical Analysis, 43(5):2078-2095, 2011.
%
% Inputs:
%     x: input signal vector, length n (need not be dyadic length)
%  type: wavelet type, string (see help wfiltfn)
%    winlen: length of window in samples; Nyquist frequency is winlen/2
%    dt: sampling period (default, dt = 1)
%   opt: options structure
%	 opt.stfttype: whether to do 'normal' or 'modified' STFT (default = 'normal')
%    opt.padtype: type of padding, options: 'symmetric',
%                 'replicate', 'circular' (default = 'symmetric')
%    opt.rpadded: return padded Sx and dSx?  (default = 1)
%    opt.type, opt.s, opt,mu, etc: window options (see help wfiltfn)
%	
%
% Outputs:
%    Sx: [na x n] size matrix (rows = scales, cols = times)
%        containing samples of the CWT of x.
%    Sfs: vector containing the associated frequencies
%   dSx (calculated if requested): [na x n] size matrix containing
%       samples of the time-derivatives of the STFT of x.
%
%
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Eugene Brevdo, Gaurav Thakur
%---------------------------------------------------------------------------------
function [Sx,Sfs,dSx] = stft_fw(x, dt, opt, t)

if nargin<3, opt = struct(); end


%opt.window is window length; opt.type overrides the default hamming window
if ~isfield(opt, 'stfttype'), opt.stfttype = 'normal'; end
if ~isfield(opt, 'winlen'), opt.winlen = round(length(x)/8); end

% options: symmeric, replicate, circular
if ~isfield(opt, 'padtype'), opt.padtype = 'symmetric'; end
if ~isfield(opt, 'rpadded'), opt.rpadded = false; end

if isfield(opt, 'type')
	windowfunc = wfiltfn(opt.type,opt,false);
	diffwindowfunc = wfiltfn(opt.type,opt,true);
end

% Pre-pad signal; this only works well for normal STFT
n = length(x);
if strcmpi(opt.stfttype,'normal')
	[x,Nold,n1,n2] = padsignal(x, opt.padtype, opt.winlen);
	n1 = floor(n1/2);
else
	n1 = 0;
end
N = length(x);

if strcmpi(opt.stfttype,'normal')
	%set up window
	if isfield(opt, 'type')
		window = windowfunc(linspace(-1,1,opt.winlen));
		diffwindow = diffwindowfunc(linspace(-1,1,opt.winlen));
	else
		window = hamming(opt.winlen);
		diffwindow = [diff(hamming(opt.winlen));0];
	end
	%window = window / norm(window);
	%diffwindow = diffwindow / norm(diffwindow)/dt;
	
	%frequency range
	Sfs = linspace(0,1,opt.winlen+1);
	Sfs = Sfs(1:floor(opt.winlen/2)+1) / dt;

	x = x(:).';			% Turn into row vector
	%compute STFT and keep only the positive frequencies
	xbuf = buffer(x, opt.winlen, opt.winlen-1, 'nodelay');
	xbuf = diag(sparse(window)) * xbuf;
	Sx = fft(xbuf,[],1);
	Sx = Sx(1:floor(opt.winlen/2)+1, :) / sqrt(N);

	%same steps for STFT derivative
	dxbuf = buffer(x, opt.winlen, opt.winlen-1, 'nodelay');
	dxbuf = diag(sparse(diffwindow)) * dxbuf;
	dSx = fft(dxbuf,[],1);
	dSx = dSx(1:floor(opt.winlen/2)+1, :) / sqrt(N);
	dSx = dSx/dt;
	
elseif strcmpi(opt.stfttype,'modified')
%	modified STFT is more accurately done in the frequency domain,
%	like a filter bank over different frequency bands
%	uses a lot of memory, so best used on small blocks (<5000 samples) at a time
	Sfs = linspace(0,1,N)/dt;
	Sx = zeros(N,N);
	dSx = zeros(N,N);
	
	halfN = round(N/2);
	halfwin = floor((opt.winlen-1)/2);
	window = windowfunc(linspace(-1,1,opt.winlen)).';
	diffwindow = diffwindowfunc(linspace(-1,1,opt.winlen)).' * 2/opt.winlen/dt;
	for k=[1:N]
		freqs = [-min([halfN-1,halfwin,k-1]):min([halfN-1,halfwin,N-k])];
		indices = mod(freqs,N)+1;
		Sx(indices,k) = x(k+freqs).*window(halfwin+freqs+1);
		dSx(indices,k) = x(k+freqs).*diffwindow(halfwin+freqs+1);
	end
	Sx = fft(Sx)/sqrt(N);
	dSx = fft(dSx)/sqrt(N);

	%only keep the positive frequencies
	Sx = Sx(1:halfN,:);
	dSx = dSx(1:halfN,:);
	Sfs = Sfs(1:halfN);
end

% 	Shorten Sx to proper size (remove padding)
if (~opt.rpadded)
	Sx = Sx(:, n1+1:n1+n);
	dSx = dSx(:, n1+1:n1+n);
end