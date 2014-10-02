%x=log(Daily(1:end,2));
%t=DailyTimeUniform(1:end);
%t=[1:3000];
%x=cos(2*pi/10*(t));
t=linspace(0,10,1000);
%x=(1+0.6*cos(2*t)).*cos(4*pi*t+1.2*t.^2);
x=(2+cos(t)).*cos(2*pi*(7*t+0.2*t.^2.8));
IFx=7+0.2*2.8*t.^(1.8);

%x=CTModulate(randi([0,1],[1,100]),struct('Type','BFSK','ModIndex',1,'BT',1,'df',320000/2),10,3200000,0.5,0.5);
%x=exp(i*round(unwrap(angle(x))/1)*1);
%t=[1:length(x)]/3200000;

%t=linspace(1,5,3000)/10^3;
%x=exp(2*pi*i*t.^2*10^6);
%x=exp(2*pi*i*round(unwrap(angle(x))*8/(2*pi))/8);

x=x(:);
CWTopt=struct('gamma',10^-15,'type','bump','mu',pi,'s',pi,'om',0,'dtype','direct','padtype','symmetric','rpadded',true,'squeezing','full','freqscale','log');
STFTopt=struct('gamma',10^-4,'type','gauss','mu',0,'s',.1,'om',0,'padtype','symmetric','domain','time','winlen',512,'rpadded',false,'squeezing','full','freqdiv',8);
nv=64;
dt=t(2)-t(1);


%[Wx,as,dWx] = cwt_fw(x, CWTopt.type, nv, dt, CWTopt);
%xCWT = cwt_iw(Wx, CWTopt.type, CWTopt, length(x), mean(x)).';

[Tx, fs, Wx, as, Cw] = synsq_cwt_fw(t, x-mean(x), nv, CWTopt);
xSST = synsq_cwt_iw(Tx, fs, CWTopt,t).';

%[Sx,Sfs,dSx] = stft_fw(x, dt, STFTopt);
%xSTFT = stft_iw(Sx, STFTopt).';

%[Tx, fs, Sx, Sfs, Sw, dSx] = synsq_stft_fw(t, x, STFTopt);

plot([x,xSST])
%figure(); plot([x,xSTFT])
%figure(); plot([x-mean(x),xSST])


% surf(abs(Sx),'edgecolor','none'); view(0,90); axis tight;
%[V,B]=min(abs(fs-20));
%figure(); tplot(Tx, t, fs);
%colorbar('location','north');


%b = - (sum(x)*sum(x.*xSST)-sum(x.^2)*sum(xSST))/(sum(x)^2-length(x)*sum(x.^2));
%a = - sum(xSST)/sum(x) - length(x)/sum(x)*b;
