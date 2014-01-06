% ISTFT  Inverse short-time Fourier transform
%
% Very closely based on Steven Schimmel's stft.m and istft.m from
% his SPHSC 503: Speech Signal Processing course at Univ. Washington.
% Adapted for use with Synchrosqueezing Toolbox
%
% Inputs:
%  Sx: wavelet transform of a signal, see help stft_fw
%  type: wavelet used to take the STFT,
%        see help stft_fw and help wfiltfn
%  opt: options structure used for forward wavelet transform.
%
% Output:
%  x: the signal, as reconstructed from Sx
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Eugene Brevdo, Gaurav Thakur
%---------------------------------------------------------------------------------
function x = stft_iw(Sx, fs, opt)

if nargin<2, opt = struct(); end

%opt.window is window length, opt.type overrides default hamming window
if ~isfield(opt, 'winlen'), opt.winlen = round(size(Sx,2)/16); end
if ~isfield(opt, 'overlap'), opt.overlap = opt.winlen-1; end
if ~isfield(opt, 'rpadded'), opt.rpadded = false; end

if isfield(opt, 'type')
  A = wfiltfn(opt.type,opt);
  window = A(linspace(-1,1,opt.winlen));
else
  window = hamming(opt.winlen);
end
%window = window / norm(window, 2); % Unit norm
Nwin = length(window);

%find length of padding, similar to outputs of padsignal
n = size(Sx, 2);
if (~opt.rpadded)
  xLen=n;
else
  xLen=n-Nwin;
end
nup = xLen+2*Nwin;
n1 = Nwin-1;
n2 = Nwin;
newn1 = floor((n1-1)/2);

% add STFT padding if it doesn't exist
if (~opt.rpadded)
  Sxp = zeros(size(Sx));
  Sxp(:, newn1+1:newn1+n) = Sx;
  Sx = Sxp; clear Sxp;
else
  n=xLen;
end

% regenerate the full spectrum 0...2pi (minus zero Hz value)
Sx = [Sx; conj(Sx(floor((Nwin+1)/2):-1:2,:))];

% take the inverse fft over the columns
xbuf = real(ifft(Sx,[],1));

% apply the window to the columns
xbuf = xbuf .* repmat(window(:),1,size(xbuf,2));

% overlap-add the columns
x = unbuffer(xbuf,Nwin,opt.overlap);

% Keep the unpadded part only
x = x(n1+1: n1+n);

%compute L2 norm of window to normalize inverse STFT with
windowfunc = wfiltfn(opt.type,opt,false);
C = quadgk(@(x) windowfunc(x).^2, -Inf, Inf);
%quadgk is a bit inaccurate with the bump function, this scales it correctly
if strcmpi(opt.type,'bump')
  C=C*0.8675;
end
x = 2/(pi*C)*x;

%%% subfunction
function y = unbuffer(x,w,o)
% UNBUFFER  undo the effect of 'buffering' by overlap-add (see BUFFER)
%    A = UNBUFFER(B,WINDOWLEN,OVERLAP) returns the signal A that is
%    the unbuffered version of B.

y    = [];
skip = w - o;
N    = ceil(w/skip);
L    = (size(x,2) - 1) * skip + size(x,1);

% zero pad columns to make length nearest integer multiple of skip
if size(x,1)<skip*N, x(skip*N,end) = 0; end;

% selectively reshape columns of input into 1-d signals
for i = 1:N
    t = reshape(x(:,i:N:end),1,[]);
    l = length(t);
    y(i,l+(i-1)*skip) = 0;
    y(i,[1:l]+(i-1)*skip) = t;
end;

% overlap-add
y = sum(y,1);
y = y(1:L);
