% function [Cs,freqband,bdTimes,LeftoverTx] = curve_ext_max(Tx, NumCurves, Startband, MaxStepSize, opt)
%
% Alternate, heuristic method for iterative curve extraction in the time-frequency plane,
% based on finding peak energies in the entire T-F plane, building curves by moving
% forward/backward from peaks and stopping when the components fade in or out
%
% 1. G. Thakur, E. Brevdo, N.-S. Fučkar, and H.-T. Wu,
% "The Synchrosqueezing algorithm for time-varying spectral analysis: robustness
%  properties and new paleoclimate applications," Signal Processing, 93:1079-1094, 2013.
%
% Input:
%   Tx, See help synsq_cwt_fw
%	NumCurves: Number of oscillatory components to extract
%	StartBand: initial width of band around center frequency, in number of frequency bins
%	MaxStepSize: maximum change in center frequency when it goes forward or backward in time
%   opt: options structure (see help synsq_cwt_fw)
%
% Output:
%   Cs: center frequencies of extracted components
%	freqband: width of frequency band of each component at each time
%	bdTimes: starting/stopping time of each component
%	LeftoverTx: residual Synchrosqueezing plot, after removing the extracted components
%
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: Gaurav Thakur
%---------------------------------------------------------------------------------

function [Cs,freqband,bdTimes,LeftoverTx] = curve_ext_max(Tx, NumCurves, Startband, MaxStepSize, opt)

[nFreqs, nTimes] = size(Tx);
if nargin<3, NumCurves=6; end
if nargin<4, Startband=floor(nFreqs/20); end
if nargin<5, MaxStepSize=2; end
if nargin<6, opt=struct(); end
if ~isfield(opt,'gamma') opt.gamma=10^-13; end

Cs=zeros(nTimes,NumCurves);
Es=zeros(nTimes,1);
bdTimes=repmat([1,nTimes].',[1,NumCurves]);	%initially set to whole domain
freqband=ones(nTimes,NumCurves)*Startband;
Tx(1,:)=0;				%lowest frequency may be inaccurate due to CWT, so ignore it

LeftoverTx=Tx;
for m=1:NumCurves
	%find max amplitude point to start from
	[V,B1]=max(max(abs(LeftoverTx(1:end,:)),[],2),[],1);
	[V,B2]=max(max(abs(LeftoverTx(1:end,:)),[],1),[],2);
	StartTime=B2; Cs(StartTime,m)=B1;
	LeftoverTx([max(Cs(StartTime,m)-freqband(StartTime,m),1):min(Cs(StartTime,m)+freqband(StartTime,m),nFreqs)],StartTime)=0;

	%go forward/right
	for n=[StartTime+1:1:nTimes]
		SelectedTx = zeros(nFreqs,1);
		SelectedRange = [max(Cs(n-1,m)-MaxStepSize,1):min(Cs(n-1,m)+MaxStepSize,nFreqs)];
		SelectedTx(SelectedRange) = LeftoverTx(SelectedRange,n);
		[V,B]=max(abs(SelectedTx));
		if (V<opt.gamma)			%if all amplitudes are below noise floor at this time, stop this curve
			Cs(n,m)=0;
			bdTimes(2,m)=n;
			freqband(n:end,m)=0;
			break;
		end
		Cs(n,m)=B;
		LeftoverTx([max(Cs(n,m)-freqband(n,m),1):min(Cs(n,m)+freqband(n,m),nFreqs)],n)=0;	%later curves should ignore all freqs of this one
	end

	%go backward/left
	for n=[StartTime-1:-1:1]
		SelectedTx = zeros(nFreqs,1);
		SelectedRange = [max(Cs(n+1,m)-MaxStepSize,1):min(Cs(n+1,m)+MaxStepSize,nFreqs)];
		SelectedTx(SelectedRange) = LeftoverTx(SelectedRange,n);
		[V,B]=max(abs(SelectedTx));
		if (V<opt.gamma)			%if all amplitudes are below noise floor at this time, stop this curve
			Cs(n,m)=0;
			bdTimes(1,m)=n;
			freqband(1:n,m)=0;
			break;
		end
		Cs(n,m)=B;
		LeftoverTx([max(Cs(n,m)-freqband(n,m),1):min(Cs(n,m)+freqband(n,m),nFreqs)],n)=0;	%later curves should ignore all freqs of this one
	end
%{
	Es(1)=abs(LeftoverTx(Cs(1,m),1));
	meanEs=Es(1);
	bandfactor=2;
	bandsize=round(bandfactor*meanEs/Es(n-1));		%intensity of last point compared to historical average
		Es(n)=abs(LeftoverTx(Cs(n,m),n));
		meanEs=(n*meanEs+Es(n))/(n+1);
	freqband(n)=round(freqband(n-1)*meanEs/Es(n-1));
%}

end