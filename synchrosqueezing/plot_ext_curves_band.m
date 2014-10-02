% Plot the curves extracted via curve_ext or curve_ext_multi
%
%  hc = plot_ext_curves(t, x, Tx, fs, Cs, Es, opt, clwin)
%
% Input:
%  t, x, opt: same as input to synsq_cwt_fw/iw
%  Tx, fs, clwin: same as input to curve_ext_multi
%  Cs, Es: same as output of curve_ext_multi
% Output:
%  hc: Object handle for the resulting figure
%
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Eugene Brevdo, Gaurav Thakur
%---------------------------------------------------------------------------------
function hc = plot_ext_curves_band(t, x, Tx, fs, Cs, freqband, opt)

if nargin<6, freqband = 4; end
if nargin<7, opt = struct(); end
if size(freqband,1)==1, fullfreqband=freqband*ones(size(Cs)); else fullfreqband=freqband; end;

if ~isfield(opt, 'markers')
  opt.markers = {'+','o','*','x','s','d','^','v','>','<','p','h' };
end

if ~isfield(opt,'band')
  opt.band = 1;
end

Nc = size(Cs, 2);
[na, N] = size(Tx);

assert(N == length(t));

hca = gcf;

UpperCs=fs(min(max(Cs+fullfreqband,1),length(fs)));
LowerCs=fs(min(max(Cs-fullfreqband,1),length(fs)));
for m=1:Nc
	UpperCs(find(Cs(:,m)==0),m)=NaN;
	LowerCs(find(Cs(:,m)==0),m)=NaN;
end

hold on
	hc=tplot(Tx, t, fs, opt);
	plot(t,UpperCs, 'LineWidth', 2);
	plot(t,LowerCs, 'LineWidth', 2);
hold off

