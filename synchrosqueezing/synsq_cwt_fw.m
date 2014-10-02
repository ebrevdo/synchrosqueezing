% function [Tx, fs, Wx, as, w] = synsq_cwt_fw(t, x, nv, opt)
% Quick alternatives:
%    [Tx, fs] = synsq_cwt_fw(t, x, nv)
%    [Tx, fs, Wx, as] = synsq_cwt_fw(t, x, nv)
%    [Tx, fs, Wx, as, w] = synsq_cwt_fw(t, x, nv)
%
% Calculates the synchrosqueezing transform of vector x, with
% samples taken at times given in vector t.  Uses nv voices.  opt
% is a struct of options for the way synchrosqueezing is done, the
% type of wavelet used, and the wavelet parameters.  This
% implements the algorithm described in Sec. III of [1].
%
% 1. G. Thakur, E. Brevdo, N.-S. Fučkar, and H.-T. Wu,
% "The Synchrosqueezing algorithm for time-varying spectral analysis: robustness
%  properties and new paleoclimate applications," Signal Processing, 93:1079-1094, 2013.
%
% 2. I. Daubechies, J. Lu, H.T. Wu, "Synchrosqueezed Wavelet Transforms: an
% empricial mode decomposition-like tool", Applied and Computational Harmonic Analysis
% 30(2):243-261, 2011.
%
% Input:
%  t: vector of times samples are taken (e.g. linspace(0,1,n))
%  x: vector of signal samples (e.g. x = cos(20*pi*t))
%  nv: number of voices (default: 16, recommended: 32 or 64)
%  opt: struct of options
%    opt.type: type of wavelet (see help wfiltfn)
%      opt.s, opt.mu, etc (wavelet parameters: see opt from help wfiltfn)
%    opt.disp: display debug information?
%    opt.gamma: wavelet hard thresholding value (see help cwt_freq_direct)
%    opt.dtype: 'direct' (default), 'phase' or 'numerical' differentiation.
%				'numerical' uses MEX differentiation, which is faster and
%               uses less memory, but may be less accurate.
%    
% Output:
%  Tx: synchrosqueezed output of x (columns associated with time t)
%  fs: frequencies associated with rows of Tx
%  Wx: wavelet transform of x (see cwt_fw)
%  as: scales associated with rows of Wx
%  w:  phase transform for each element of Wx
%
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Eugene Brevdo, Gaurav Thakur
%---------------------------------------------------------------------------------
function [Tx, fs, Wx, as, w] = synsq_cwt_fw(t, x, nv, opt)

if nargin<4, opt = struct(); end
if nargin<3, nv = 16; end
if nargin<2, error('Too few input arguments'); end

% Choose some default values
% Wavelet type
if ~isfield(opt, 'type'), opt.type = 'morlet'; end
% Display bugging information?
if ~isfield(opt, 'disp'), opt.disp = 0; end
% Direct or numerical differentiation? (default: direct)
if ~isfield(opt, 'dtype'), opt.dtype = 'direct'; end

if ~isfield(opt, 'riskshrink'), opt.riskshrink = false; end

% Calculate sampling period, assuming regular spacing
dt = t(2)-t(1);

% Check for uniformity of spacing
if any(diff(t,2)/(t(end)-t(1))>1e-5)
    error('time vector t is not uniformly sampled');
end


% Calculate the wavelet transform - padded via symmetrization
x = x(:);
N = length(x);
[Nup,n1,n2] = p2up(N);

if strcmpi(opt.dtype,'direct')
    % Calculate derivative directly in the wavelet domain before taking
    % wavelet transform. (default)

    opt.rpadded = 0;
    [Wx,as,dWx] = cwt_fw(x, opt.type, nv, dt, opt);
    w = phase_cwt(Wx, dWx, opt);

    clear dWx;
elseif strcmpi(opt.dtype,'phase')
	% take derivative of unwrapped CWT phase directly in phase transform

    opt.rpadded = 0;
	[Wx,as] = cwt_fw(x, opt.type, nv, dt, opt);
    w = phase_cwt(Wx, [], opt);
else
    % Calculate derivative numerically after calculating wavelet
    % transform.  This requires less memory and is more accurate at
    % small a.

    opt.rpadded = 1;
    [Wx,as] = cwt_fw(x, opt.type, nv, dt, opt);
    Wx = Wx(:, n1-4+1:n1+N+4);
    w = phase_cwt_num(Wx, dt, opt);
end

% choose gamma automatically using riskshrink method, overriding any specified gamma
if (opt.riskshrink)
	opt.gamma = est_riskshrink_thresh(Wx, nv);
end

% Calculate the synchrosqueezed frequency decomposition
opt.transform='CWT';
[Tx,fs] = synsq_squeeze(Wx, w, t, nv, opt);

if strcmpi(opt.dtype,'numerical')
    Wx = Wx(:, 4+1:end-4);
    w = w(:, 4+1:end-4);
    Tx = Tx(:, 4+1:end-4);
end

end
