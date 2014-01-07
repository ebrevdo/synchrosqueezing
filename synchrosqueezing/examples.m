%Several quick examples illustrating how the toolbox functions are used
%Uncomment each set of lines to use them

%x is the signal
%xNew is the inverted (reconstructed) signal
t=linspace(0,20,2000);
x=cos(2*pi*(0.1*t.^2.6+3*sin(2*t)+10*t));
x=x(:);
dt=t(2)-t(1);


%various options and parameters
CWTopt=struct('gamma',10^-15,'type','bump','mu',6,'s',3,'om',0,'nv',64,'freqscale','linear');
STFTopt=struct('gamma',10^-15,'type','gauss','mu',0,'s',0.1,'om',0,'winlen',128,'squeezing','full');

% Continuous wavelet transform (CWT)
%{
[Wx,as,dWx] = cwt_fw(x, CWTopt.type, CWTopt.nv, dt, CWTopt);
xNew = cwt_iw(Wx, CWTopt.type, CWTopt, length(x), mean(x)).';
figure(); tplot(Wx, t, as); colorbar; title('CWT');
%figure(); plot(t,[x,xNew]);
%}

% Short-time Fourier transform (STFT)
%{
[Sx,fs,dSx] = stft_fw(x, dt, STFTopt);
xNew = stft_iw(Sx, fs, STFTopt).';
figure(); tplot(Sx, t, fs); colorbar; title('STFT');
%figure(); plot(t,[x,xNew]);
%}

% CWT Synchrosqueezing transform

[Tx, fs, Wx, as, Cw] = synsq_cwt_fw(t, x-mean(x), CWTopt.nv, CWTopt);
xNew = synsq_cwt_iw(Tx, fs, CWTopt).';
figure(); tplot(Tx, t, fs); colorbar; title('CWT Synchrosqueezing');
%figure(); plot(t,[x,xNew]);


% STFT Synchrosqueezing transform
%{
[Tx, fs, Sx, Sfs, Sw, dSx] = synsq_stft_fw(t, x, STFTopt);
xNew = synsq_stft_iw(Tx, fs, STFTopt).';
figure(); tplot(Tx, t, fs); colorbar; title('STFT Synchrosqueezing');
%figure(); plot(t,[x,xNew])
%}
