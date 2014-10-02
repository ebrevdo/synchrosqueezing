function [t,x,tint,xintrt,p,Txint,fsint,Wxint,asint] = interp_synsq_file(fn, opt)

% Skip repeated values?
if ~isfield(opt, 'skiprep'), opt.skiprep = 0; end

% Data
data = textread(fn,'','delimiter',opt.delim);

% time difference in years from today
t = data(:,1);

% signal
x = data(:,2);

% reverse order
if (opt.trev)
    t = t(end:-1:1);
    x = x(end:-1:1);
end

if (opt.skiprep)
  rep = find(diff(t)==0);
  reprm = [rep(:); rep(:)+1];
  t(reprm) = [];
  x(reprm) = [];
end

% make year data regular
%
% if opt.tend is set, make that the base for the value intervals,
% and cut off the time at the true end.
t0 = t(1);
if isfield(opt, 'tend')
  tend0 = t(end);
  tend = opt.tend;
else
  tend0 = t(end);
  tend = t(end);
end
tint = tend:-opt.dt:t0;
tint(find(tint>tend0)) = [];
tint = tint(end:-1:1);

% interpolate
xint = interp1(t, x, tint, opt.intmethod);

% remove trend
[xintrt,p] = rtrend(tint, xint, opt.trord);

% take synchrosqueezing transform of temp
[Txint, fsint, Wxint, asint] = synsq_cwt_fw(tint, xintrt, opt.nv, opt);

end
