% function [xd,p] = rtrend(t, x, ord)
%
% Remove trends from x(t) of order ord.
%
% input
%    t: time vector
%    x: data vector
%  ord: polynomial fit order (integer >= 0)
%
% output
%   xd: the detrended version of x(t)
%    p: the polynomial coefficients, order ord (see help polyfit)
function [xd,p] = rtrend(t, x, ord)
    p = polyfit(t, x, ord);
    xp = polyval(p, t);
    xd = x-xp;
end
