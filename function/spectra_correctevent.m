function [spec_mat,npad0,npts,f,NFFT]=spectra_correctevent(Zraw,H1raw,H2raw,Praw,dt,taxisZ,taxis1,taxis2,taxisP,taptime);

% Spectra parameters for resampled data
npts = length(Zraw);
Ppower = nextpow2(npts);
NFFT = 2^Ppower;
samprate = 1/dt; %common named Fs
f = samprate/2*linspace(0,1,NFFT/2+1);
pts_begin = 1;
pts_end = length(Zraw);
npad0 = floor((NFFT-npts)/2);

amp_Z  = Zraw(pts_begin:pts_end);
amp_Z  = amp_Z.*flat_hanning(taxisZ(pts_begin:pts_end),taptime*dt*npts);
amp_Z  = detrend(amp_Z,0);
amp_Z  = padarray(amp_Z,[npad0 0],'both');
spectrum = fft(amp_Z,NFFT).*dt;
spec_Z = spectrum(1:NFFT/2+1);

amp_H1  = H1raw(pts_begin:pts_end);
amp_H1  = amp_H1.*flat_hanning(taxis1(pts_begin:pts_end),taptime*dt*npts);
amp_H1  = detrend(amp_H1,0);
amp_H1  = padarray(amp_H1,[npad0 0],'both');
spectrum = fft(amp_H1,NFFT).*dt;
spec_H1 = spectrum(1:NFFT/2+1);

amp_H2  = H2raw(pts_begin:pts_end);
amp_H2  = amp_H2.*flat_hanning(taxis2(pts_begin:pts_end),taptime*dt*npts);
amp_H2  = detrend(amp_H2,0);
amp_H2  = padarray(amp_H2,[npad0 0],'both');
spectrum = fft(amp_H2,NFFT).*dt;
spec_H2 = spectrum(1:NFFT/2+1);

amp_P  = Praw(pts_begin:pts_end);
amp_P  = amp_P.*flat_hanning(taxisP(pts_begin:pts_end),taptime*dt*npts);
amp_P  = detrend(amp_P,0);
amp_P  = padarray(amp_P,[npad0 0],'both');
spectrum = fft(amp_P,NFFT).*dt;
spec_P = spectrum(1:NFFT/2+1);

spec_mat = [spec_Z,spec_H1,spec_H2,spec_P];

return