% The inverse wavelet transform of signal Wx
%
% Implements Eq. (4.67) of [1].
%
% 1. Mallat, S., Wavelet Tour of Signal Processing 3rd ed.
%
% 2. G. Thakur, E. Brevdo, N.-S. Fuckar, and H.-T. Wu,
% "The Synchrosqueezing algorithm for time-varying spectral analysis: robustness
%  properties and new paleoclimate applications," Signal Processing, 93:1079-1094, 2013.
%
% Inputs:
%  Wx: wavelet transform of a signal, see help cwt_fw
%  type: wavelet used to take the wavelet transform,
%        see help cwt_fw and help wfiltfn
%  opt: options structure used for forward wavelet transform.
%  xLen: length of original x, before padding
%  xMean: mean of original x (not picked up in CWT since it's infinite scale component)
%
% Output:
%  x: the signal, as reconstructed from Wx
%
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Eugene Brevdo, Gaurav Thakur
%---------------------------------------------------------------------------------
function x = cwt_iw(Wx, type, opt, xLen, xMean)
	if nargin<4, xMean = 0; end
    if nargin<3, opt = struct(); end

    [na, n] = size(Wx);
	if nargin<4, xLen = n; end
	[N,n1,n2] = p2up(xLen);
	
    % options: symmeric, replicate, circular
    if ~isfield(opt, 'padtype'), opt.padtype = 'symmetric'; end
    if ~isfield(opt, 'rpadded'), opt.rpadded = 0; end

%	add CWT padding if it doesn't exist
	if ~(opt.rpadded)
		Wxp = zeros(na, N);
		Wxp(:, n1+1:n1+n) = Wx;
		Wx = Wxp; clear Wxp;
	else
		n=xLen;
	end

    % Following the same value in cwt_fw
    noct = log2(N)-1;
    nv = na/noct;
    as = 2^(1/nv) .^ (1:1:na);
    
    assert(mod(noct,1) == 0);
    assert(nv>0 && mod(nv,1)==0); % integer

    % Find the admissibility coefficient Cpsi
    Cpsi = cwt_adm(type, opt);

    x = zeros(1, N);
    for ai=1:na
        a = as(ai);
        Wxa = Wx(ai, :);

        psih = wfilth(type, N, a, opt);

        % Convolution theorem here
        Wxah = fft(Wxa);
        xah = Wxah .* psih;
        xa = ifftshift(ifft(xah));
        
        x = x + xa/a;
    end

     % Take real part and normalize by log_e(a)/Cpsi
     x = log(2^(1/nv))/Cpsi * real(x);

	 %add on mean
 	 x = x+xMean;
	 
     % Keep the unpadded part
     x = x(n1+1: n1+n);

end % cwt_iw
