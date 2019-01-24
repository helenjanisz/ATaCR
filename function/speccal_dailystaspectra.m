function[spectrum_C,cspectrum_C,is_goodwin] = speccal_dailystaspectra(nwin,f,iptwin1,iptwin2,goodwins,datraw,comp_exist,taxisC,overlap,dt,npts,npad0,NFFT)
% calculates spectral properties for a given component

if comp_exist ==1
    spectrum_C=zeros(nwin,length(f));
    cspectrum_C=zeros(nwin,length(f));
    
    for iwin = 1:length(iptwin1) % windowing the daily records in this script to start, assumes all time vectors are aligned
        pts_begin = iptwin1(iwin);
        pts_end = iptwin2(iwin);
        
        if ~isempty(intersect(iwin,goodwins))
            is_goodwin(iwin) = 1;
        else
            is_goodwin(iwin) = 0;
        end
        
        amp_C  = datraw(pts_begin:pts_end);
        amp_C  = amp_C.*flat_hanning(taxisC(pts_begin:pts_end),overlap*dt*npts);
        amp_C  = detrend(amp_C,0);
        amp_C  = padarray(amp_C,[npad0 0],'both');
        spectrum = fft(amp_C,NFFT).*dt;
        spec_C = spectrum(1:NFFT/2+1);
        cspec_C = conj(spec_C);
        spectrum_C(iwin,1:length(f)) = spec_C';
        cspectrum_C(iwin,1:length(f)) = cspec_C';
        
    end
else
    spectrum_C = [];
    cspectrum_C = [];
    is_goodwin = [];
end

return