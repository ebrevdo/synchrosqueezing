% function gamma = est_riskshrink_thresh(Wx, nv)
%
% Estimate the RiskShrink hard thresholding level
%
% Implements Eq. of [1].
%
% 1. TODO, donoho & johnstone.
%
% Inputs:
%  Wx: wavelet transform of a signal, see help cwt_fw
%  opt: options structure used for forward wavelet transform.
%
% Output:
%  gamma: the RiskShrink hard threshold estimate
%
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/)
%---------------------------------------------------------------------------------
function gamma = est_riskshrink_thresh(Wx, nv)
    % TODO, error if opt has no nv, or no opt.
    if nargin<2, opt = struct(); end

    [na, n] = size(Wx);

    Wx_fine = abs(Wx(1:nv, :));
    
    gamma = sqrt(2*log(n)) * mad(Wx_fine(:)) * 1.4826;
    
end % est_riskshrink_thresh
