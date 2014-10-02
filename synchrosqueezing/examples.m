%Several quick examples illustrating how the toolbox functions are used
%Uncomment each set of lines to use them

%x is the signal
%xNew is the inverted (reconstructed) signal
t=linspace(0,10,2000);
x=cos(2*pi*(0.1*t.^2.6+3*sin(2*t)+10*t)) + exp(-0.2*t).*cos(2*pi*(40+t.^1.3).*t);
x=x(:);
dt=t(2)-t(1);


%various options and parameters
CWTopt=struct('gamma',eps,'type','morlet','mu',6,'s',2,'om',0,'nv',64,'freqscale','linear');
STFTopt=struct('gamma',eps,'type','gauss','mu',0,'s',0.05,'om',0,'winlen',256,'squeezing','full');

% Continuous wavelet transform (CWT)
%{
[Wx,as,dWx] = cwt_fw(x, CWTopt.type, CWTopt.nv, dt, CWTopt);
xNew = cwt_iw(Wx, CWTopt.type, CWTopt, length(x), mean(x)).';
tplot(Wx, t, as); colorbar; title('CWT','FontSize',14); xlabel('Time (seconds)','FontSize',14); ylabel('Frequency (hz)', 'FontSize',14);
figure(); plot(t,[x,xNew]);
%}

% Short-time Fourier transform (STFT)
%{
[Sx,fs,dSx] = stft_fw(x, dt, STFTopt);
xNew = stft_iw(Sx, fs, STFTopt).';
tplot(Sx, t, fs); colorbar; title('STFT','FontSize',14); xlabel('Time (seconds)','FontSize',14); ylabel('Frequency (hz)', 'FontSize',14);
figure(); plot(t,[x,xNew]);
%}

% CWT Synchrosqueezing transform
%{
[Tx, fs, Wx, as, Cw] = synsq_cwt_fw(t, x-mean(x), CWTopt.nv, CWTopt);
xNew = synsq_cwt_iw(Tx, fs, CWTopt).';
tplot(Tx, t, fs); colorbar; title('CWT Synchrosqueezing','FontSize',14); xlabel('Time (seconds)','FontSize',14); ylabel('Frequency (hz)', 'FontSize',14);
figure(); plot(t,[x,xNew(:,1)]);
%}

% STFT Synchrosqueezing transform
%{
[Tx, fs, Sx, Sfs, Sw, dSx] = synsq_stft_fw(t, x, STFTopt);	
xNew = synsq_stft_iw(Tx, fs, STFTopt).';
figure(); tplot(Tx, t, fs); colorbar; title('STFT Synchrosqueezing','FontSize',14); xlabel('Time (seconds)','FontSize',14); ylabel('Frequency (hz)', 'FontSize',14);
figure(); plot(t,[x,xNew(:,1)])
%}


% Blind source separation using STFT Synchrosqueezing transform
% computes Synchrosqueezing with high threshold, extracts curves in time-frequency plane
% then repeats transform with lower threshold and inverts the transform
%{
STFTopt.gamma = 10^-3;
[Tx, fs, Sx, Sfs, Sw, dSx] = synsq_stft_fw(t, x, STFTopt);	
[Cs,freqband,bdTimes,LeftoverTx] = curve_ext_max(Tx, 2, 15, 4, STFTopt);

STFTopt.gamma = eps;
[Tx, fs, Sx, Sfs, Sw, dSx] = synsq_stft_fw(t, x, STFTopt);
figure(); plot_ext_curves_band(t, x, Tx, fs, Cs, freqband, STFTopt); colorbar; title('STFT Synchrosqueezing','FontSize',14); xlabel('Time (seconds)','FontSize',14); ylabel('Frequency (hz)', 'FontSize',14);

% invert around masks formed from curves
xNew = synsq_stft_iw(Tx, fs, STFTopt, Cs, freqband).';
figure(); plot(t,[xNew(:,1:end-1)]); title('Separated components','FontSize',14); axis tight; xlabel('Time (seconds)','FontSize',14);
figure(); plot(t,[xNew(:,1:end-1),x]); title('Separated components with original signal','FontSize',14); axis tight; xlabel('Time (seconds)','FontSize',14);
P=get(gca, 'children'); set(gca,'children',[P(2:3);P(1)]);
%}