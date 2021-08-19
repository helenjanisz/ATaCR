function [ ] = plot_corrspectra( T1,T2,Zraw,corrseis,f,NFFT,dt,taxis,eventid,netsta,lp,hp )
% Plot spectra and corrected seismograms;

figure(1); clf;
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);
ilgd = 0;

subplot(3,1,1);  box on;
psd_Z = disp2accel_psd(Zraw,NFFT,dt);
ilgd = ilgd+1;
h1(ilgd) = semilogx(f,10*log10(psd_Z),'-k','linewidth',0.5); hold on;
lgd{ilgd} = 'Z';

clrs = lines(length(corrseis));
psd_band_avg = zeros(length(corrseis),1);
for icor = 1:length(corrseis)
    ilgd = ilgd+1;
    psd_Zcorr = disp2accel_psd(corrseis(icor).timeseries,NFFT,dt);

    h1(ilgd) = semilogx(f,10*log10(psd_Zcorr),'-','color',clrs(icor,:),'linewidth',0.5); hold on;
    lgd{ilgd} = corrseis(icor).label;
    
    % Calculate PSD within band of interest
    Ifband = (f>=1/T2 & f<=1/T1);
    psd_band_avg(icor) = mean(10*log10(psd_Zcorr(Ifband)));
end

% ---------- plot HLNM Peterson, 1993 ----------
[LL,HH,FF] = noise_models(100); 
semilogx(FF,HH,'color',[0.9 0.9 0.9]); hold on; 
semilogx(FF,LL,'color',[0.9 0.9 0.9]); hold on;
xx=[FF' fliplr(FF')];
yy=[HH' fliplr(LL')];
fill(xx,yy,[0.9 0.9 0.9],'FaceAlpha', 0.5,'EdgeColor','none');
%fill(xx,yy,[0.9 0.9 0.9],'EdgeColor','none');	  
grid on;set(gca,'layer','top');
xlabel('Frequency (Hz)');
ylabel('Power[DB], 10 \ast log[m^{2}/s^{4}/Hz]');
title(sprintf ('%s %s',netsta,eventid ));
xlim([1/250,max(f)]);
ylim([-200 -80]);

% Plot bounds of TF filter
semilogx([lp lp],[-200 -80],'-r')
semilogx([hp hp],[-200 -80],'-r')

legend(h1,lgd,'location','eastoutside');


% PLOT DATA
fn = 1/2/dt;
[b,a]=butter(2,[1/fn/T2,1/fn/T1]);
Z_filt  = filtfilt(b,a,Zraw);

% Find component with lowest PSD in frequency band of interest
[~,i_best] = min(psd_band_avg);
ZP21_filt  = filtfilt(b,a,corrseis(i_best).timeseries);

% Zraw data
subplot(3,1,2); box on;
plot(taxis,Z_filt,'-k')
xlabel('Time (s)');
title('Z raw');
xlim([min(taxis),max(taxis)]);

% ZP-21 data
subplot(3,1,3); box on;
plot(taxis,ZP21_filt,'-k');
xlabel('Time (s)');
title(['Lowest Noise: ',corrseis(i_best).label]);
xlim([min(taxis),max(taxis)]);

end

%%
% time_differentiate.m
% like the diff command to compute the derivative in the time domain
% and the freq_differentiate.m functions, but has the same effect in the time domain
% usage:
% out=time_differentiate(in,sps)
% where sps is the number of samples per second
% c version:
% //time derivative (as sac 2 point algorithm)
% void time_derivative(double *in, double *out, double delta, int npts) {
%	int i;
%	for (i=0; i<npts-1;i++) {
%		out[i] = (in[i+1] - in[i]) / delta;
%	}
%	//last value is not filled by sac. Here I fill it with the same as the previous value to have a minimal effect on future usage
%	out[npts-1] = out[npts-2];
%	
%	return;
% }
% //END

function out = time_differentiate(in,sps)
out = zeros(1,length(in));
for i=1:length(in) - 1
    out(i) = (in(i+1) - in(i)) * sps;
end
out(length(in)) = out(length(in) - 1);

end

%%
function [PSD_lo,PSD_hi,f] = noise_models(npts)

% MATLAB function "noise_models.m" generates Peterson's low
% and high noise models as a function of frequency.  
% Output values are the power spectral density in decibels
% re. (1 m^2/s^4/Hz), i.e., in units of acceleration.  The 
% following function is based on Fortran code provided
% by Bob Woodward.
%
% 	The low-noise model is "a composite of station spectra obtained 
% from many different instruments, vaults, geologic environments, 
% and geographic regions.  It is a hypothetical background 
% spectrum that is unlikely to be duplicated at any single location 
% on the Earth."
%
% 	The high-noise model is "a spectrum of *average* high background
% noise power in the (global) network."
%
% Reference:
%  Observations and Modeling of Seismic Background Noise, Jon Peterson,
%  Open-File Report 93-322, U.S. Department of Interior, U.S. Geological
%  Survey, Albuquerque, NM, USA, 1993,
%
% USAGE: [PSD_lo,PSD_hi,ff] = noise_models(npts);
%                                                             j.a.collins
%------------------------------------------------------------------------

PI = 4*atan(1);
T_lo = 0.1;                     % minimum period (seconds)
T_hi = 10000;                     % maximum period (seconds)
f_lo = 1/T_hi;                  % minimum frequency (Hz)
f_hi = 1/T_lo;                  % maximum frequency (Hz)

exp_lo = log10(f_lo);
exp_hi = log10(f_hi);
f = logspace(exp_lo, exp_hi, npts);
f = f(:);
w = 2*PI*f;
T = 1./f;

lnm = zeros(3,21);

%          period
lnm =  [     0.10     -162.36     5.64
             0.17     -166.70     0.00
             0.40     -170.00    -8.30
             0.80     -166.40    28.90
             1.24     -168.60    52.48
             2.40     -159.98    29.81
             4.30     -141.10     0.00
             5.00      -71.36   -99.77
             6.00      -97.26   -66.49
            10.00     -132.18   -31.57
            12.00     -205.27    36.16
            15.60      -37.65  -104.33
            21.90     -114.37   -47.10
            31.60     -160.58   -16.28
            45.00     -187.50     0.00
            70.00     -216.47    15.70
           101.00     -185.00     0.00
           154.00     -168.34    -7.61
           328.00     -217.43    11.90
           600.00     -258.28    26.60
         10000.00     -346.88    48.75    ]';


hnm = zeros(3,11);

%          period
hnm =  [     0.10      -108.73   -17.23
             0.22      -150.34   -80.50
             0.32      -122.31   -23.87
             0.80      -116.85    32.51
             3.80      -108.48    18.08
             4.60       -74.66   -32.95
             6.30         0.66  -127.18
             7.90       -93.37   -22.42
            15.40        73.54  -162.98
            20.00      -151.52    10.01
           354.80      -206.66    31.63    ]';


% find maximum value of period s.t. T >= period
for (n = 1:npts)
    ndx_lo = max( find(( T(n) >= lnm(1,:) )) );    
    ndx_hi = max( find(( T(n) >= hnm(1,:) )) );    
    PSD_lo(n) = lnm(2,ndx_lo) + lnm(3,ndx_lo)*log10(T(n));
    PSD_hi(n) = hnm(2,ndx_hi) + hnm(3,ndx_hi)*log10(T(n));
end

PSD_lo = PSD_lo(:);
PSD_hi = PSD_hi(:);


%PSD_lo_a = PSD_lo;                  % (dB re. 1 m**2/s**4/Hz)
%PSD_hi_a = PSD_hi;

% convert to units of velocity        (dB re. 1 m**2/s**2/Hz)
%PSD_lo_v = PSD_lo_a - 20*log10(w);
%PSD_hi_v = PSD_hi_a - 20*log10(w);
%PSD_lo = PSD_lo_v;
%PSD_hi = PSD_hi_v;

% convert to units of displacement    (dB re. 1 m**2/Hz)
%PSD_lo_d = PSD_lo_a - 40*log10(w);
%PSD_hi_d = PSD_hi_a - 40*log10(w);
%PSD_lo = PSD_lo_d;
%PSD_hi = PSD_hi_d;
% plot for DB 
%figure
%semilogx(f,PSD_hi); hold on 
%semilogx(f,PSD_lo);
%xlabel('F(hz)'); xlim([0.008 10])
%semilogx(1./f,PSD_lo)



% ===== non log. units   NOT DB ======================
%PSD_lo = 10.^(PSD_lo/10.0);         % (m**2/s**4/Hz)
%PSD_hi = 10.^(PSD_hi/10.0);


%PSD_lo_v = 10.^(PSD_lo_v/10.0);     % (m**2/s**2/Hz)
%PSD_hi_v = 10.^(PSD_hi_v/10.0);
%PSD_lo=PSD_lo_v;
%PSD_hi=PSD_hi_v;

%PSD_lo_d = 10.^(PSD_lo_d/10.0);     % (m**2/Hz)
%PSD_hi_d = 10.^(PSD_hi_d/10.0);
%PSD_lo=PSD_lo_d;
%PSD_hi=PSD_hi_d;

% plot for non log. units   seimilar to Webb 2011 BSSA, fig1
%figure
%loglog(f,PSD_hi);hold on 
%loglog(f,PSD_lo);
%xlabel('F(hz)');xlim([0.001 1]);
end

%%
function psd = disp2accel_psd(DISP,NFFT,dt)
samprate = 1/dt;
VEL= time_differentiate(DISP, samprate);
ACC= time_differentiate(VEL, samprate);
xdft = fft(ACC,NFFT);
xdft = xdft(1:NFFT/2+1);
psdx = (1/(samprate*NFFT)).*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
npts_smooth = floor(NFFT/1000)+1;
psd = smooth(psdx,npts_smooth);

end