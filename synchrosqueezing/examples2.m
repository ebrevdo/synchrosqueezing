%Several quick examples illustrating how the toolbox functions are used
%Uncomment each set of lines to use them

%x is the signal
%xNew is the inverted (reconstructed) signal
figure('Position',[400,100,1280,700]);
ha=tight_subplot(2,2,[.08 .05], [.08 .05], [.06 .01]);

t=linspace(0,20,2000);
x=cos(2*pi*(0.1*t.^2.6+3*sin(2*t)+10*t));
x=x(:);
dt=t(2)-t(1);

%{
samples=40000:50000;
x=real(data(samples)).';
t=linspace(0,0.001,50001); t=t(2:end); dt=t(2)-t(1);
t=t(samples);
%}

%various options and parameters
CWTopt=struct('gamma',10^-15,'type','morlet','mu',6,'s',5,'om',0,'nv',64,'freqscale','linear');
STFTopt=struct('gamma',10^-15,'type','gauss','mu',0,'s',0.02,'om',0,'winlen',512,'squeezing','full');

% Continuous wavelet transform (CWT)
subplot(ha(1));
[Wx,as,dWx] = cwt_fw(x, CWTopt.type, CWTopt.nv, dt, CWTopt);
xNew = cwt_iw(Wx, CWTopt.type, CWTopt, length(x), mean(x)).';
tplot(Wx, t, as); colorbar; title('CWT','FontSize',12); ylabel('Frequency (hz)', 'FontSize',12);
%figure(); plot(t,[x,xNew]);



% Short-time Fourier transform (STFT)
subplot(ha(2));
[Sx,fs,dSx] = stft_fw(x, dt, STFTopt);
xNew = stft_iw(Sx, fs, STFTopt).';
tplot(Sx, t, fs); colorbar; title('STFT','FontSize',12);
%figure(); plot(t,[x,xNew]);


% CWT Synchrosqueezing transform
subplot(ha(3));
[Tx, fs, Wx, as, Cw] = synsq_cwt_fw(t, x-mean(x), CWTopt.nv, CWTopt);
xNew = synsq_cwt_iw(Tx, fs, CWTopt).';
tplot(Tx, t, fs); colorbar; title('CWT Synchrosqueezing','FontSize',12); xlabel('Time (seconds)','FontSize',12); ylabel('Frequency (hz)', 'FontSize',12);
%figure(); plot(t,[x,xNew]);


% STFT Synchrosqueezing transform
subplot(ha(4));
[Tx, fs, Sx, Sfs, Sw, dSx] = synsq_stft_fw(t, x, STFTopt);
xNew = synsq_stft_iw(Tx, fs, STFTopt).';
tplot(Tx, t, fs); colorbar; title('STFT Synchrosqueezing','FontSize',12); xlabel('Time (seconds)','FontSize',12);
%figure(); plot(t,[x,xNew])
