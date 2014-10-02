% function [xpad, nup, n1, n2] = padsignal(x, padtype, padlength)
%
% Pads signal and returns indices of original signal
%
% Input:
%	x: original signal
%	padtype (optional): either 'symmetric' (default) or 'replicate'
%   padlength (optional): number of samples to pad on each side; default is nearest power of 2
% Output:
%	x: padded signal
%   nup: next power of 2
%   n1: length on left
%   n2: length on right
%
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Eugene Brevdo, Gaurav Thakur
%---------------------------------------------------------------------------------
function [xpad, nup, n1, n2] = padsignal(x, padtype, padlength)

[nup, n1, n2] = p2up(length(x));
if nargin<3
	[nup, n1, n2] = p2up(length(x));	%if padlength not given, pad up to nearest power of 2
else
	nup = length(x)+2*padlength;
	n1 = padlength-1;
	n2 = padlength;
end

%xl = padarray(x(:), n1, padtype, 'pre');
%xr = padarray(x(:), n2, padtype, 'post');
%xpad = [xl(1:n1); x(:); xr(end-n2+1:end)];
%padarray needs image processing toolbox; below is equivalent code to avoid that

x=x(:);
n=length(x);
if strcmpi(padtype,'symmetric')
	xl = repmat([x;flipud(x)],[ceil(n1/(2*n)),1]);
	xr = repmat([flipud(x);x],[ceil(n2/(2*n)),1]);
elseif strcmpi(padtype,'replicate')
	xl = x(1)*ones(n1,1);
	xr = x(end)*ones(n2,1);
end

xpad = [xl(end-n1+1:end); x; xr(1:n2)];