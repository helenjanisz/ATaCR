function [max_coh_fine,max_or_fine] = spec_orient(spectrum_Z,spectrum_H1,spectrum_H2,cspectrum_Z,hangs,tiltfreq,f,isfig,dayid,isgoodwin,NFFT,dt)

c = colormap('jet');

hangint = hangs(2)-hangs(1);

days = dayofyear(str2num(dayid(1:4)),str2num(dayid(5:6)),str2num(dayid(7:8)),0,0);
cc=interp1(1:64,c,((days)/(365))*63+1);

for ih = 1:length(hangs)
    hang = hangs(ih);
    cang = cos(hang*pi/180);
    sang = sin(hang*pi/180);
    
    chh_stack=zeros(1,length(f))';
    czz_stack=zeros(1,length(f))';
    chz_stack=zeros(1,length(f))';
    nwin_stack = 0;
    for ista = 1:length(isgoodwin)
        if isgoodwin(ista)==0
            continue
        end
        spec_H = sang.*spectrum_H2(ista,:)'+cang.*spectrum_H1(ista,:)'; %rotated horizontal spectra and cross-spectra
        
        chh = abs(spec_H).^2*2/(NFFT*dt);
        czz = abs(spectrum_Z(ista,:)').^2*2/(NFFT*dt);
        chz = spec_H.*(cspectrum_Z(ista,:)')*2/(NFFT*dt);
        
        chh_stack=chh_stack+chh;
        czz_stack=czz_stack+czz;
        chz_stack=chz_stack+chz;
        nwin_stack = nwin_stack +1;
    end
    
    %Normalization
    chh_stack = chh_stack/nwin_stack;
    czz_stack = czz_stack/nwin_stack;
    chz_stack = chz_stack/nwin_stack;
    
    % Coherence
    coh_stack = abs(chz_stack).^2./(chh_stack.*czz_stack);
    % Phase
    ph_stack = 180/pi.*atan2(imag(chz_stack),real(chz_stack));
    % Admittance
    ad_stack = abs(chz_stack)./chh_stack;
    
    % plots for debugging, make sure rotation is operating correctly
%     figure(95)
%     subplot(211)
%     semilogx(f,(smooth(coh_stack,40)+ih),'-k');
%     hold on
%     semilogx(f,ih*ones(size(f)),'-r');
%     xlim([0.005 0.035])
%     
%     figure(95)
%     subplot(212)
%     semilogx(f,(ph_stack+ih*180),'.k');
%     hold on
%     semilogx(f,180*ih*ones(size(f)),'-r');
%     xlim([0.005 0.035])
    % end debugging plots
    
    % looking for average coherence value within frequency range
    
    [fmin,idx_flo] = min(abs(f-tiltfreq(1)));
    [fmin,idx_fhi] = min(abs(f-tiltfreq(2)));
    av_coh_coarse(ih) = mean(coh_stack(idx_flo:idx_fhi));
    av_ph_coarse(ih) = mean(abs(ph_stack(idx_flo:idx_fhi)));
    
end

% plotting the dependence of phase and coherence on angle
if isfig ==1
    % color scaled by day of year
    figure(105)
    subplot(1,2,1); hold on
    plot(hangs,av_coh_coarse,'LineWidth',.5,'Color',cc);
    title(sprintf('Coherence, %.3f - %.3f Hz', tiltfreq(1),tiltfreq(2)));
    xlabel('Angle from H1'); ylabel ('Coherence'); ylim([0 1]); xlim([0 360])
    subplot(1,2,2); hold on
    plot(hangs,av_ph_coarse,'LineWidth',.5,'Color',cc);
    title(sprintf('Phase, %.3f - %.3f Hz', tiltfreq(1),tiltfreq(2)));
    xlabel('Angle from H1'); ylabel ('Phase'); xlim([0 360])
end

idxph = find(abs(av_ph_coarse)<90);

[max_coh_coarse,idx]=max(av_coh_coarse(idxph),[],2);
max_or_coarse = hangs(idxph(idx));

% fine grid search to get best orientation
hangs2 = [max_or_coarse-hangint+1:max_or_coarse+hangint-1];
for ih = 1:length(hangs2)
    hang = hangs2(ih);
    if hang == max_or_coarse;
        continue
    end
    cang = cos(hang*pi/180);
    sang = sin(hang*pi/180);
    
    % Initialize output structures
    chh_stack=zeros(1,length(f))';
    czz_stack=zeros(1,length(f))';
    chz_stack=zeros(1,length(f))';
    nwin_stack = 0;
    for ista = 1:length(isgoodwin)
        if isgoodwin==0
            continue
        end
        spec_H = sang.*spectrum_H2(ista,:)'+cang.*spectrum_H1(ista,:)'; %rotated horizontal spectra and cross-spectra
        
        chh = abs(spec_H).^2*2/(NFFT*dt);
        czz = abs(spectrum_Z(ista,:)').^2*2/(NFFT*dt);
        chz = spec_H.*(cspectrum_Z(ista,:)')*2/(NFFT*dt);
        
        chh_stack=chh_stack+chh;
        czz_stack=czz_stack+czz;
        chz_stack=chz_stack+chz;
        nwin_stack = nwin_stack +1;
    end
    
    %Normalization
    chh_stack = chh_stack/nwin_stack;
    czz_stack = czz_stack/nwin_stack;
    chz_stack = chz_stack/nwin_stack;
    
    % Coherence
    coh_stack = abs(chz_stack).^2./(chh_stack.*czz_stack);
    % Phase
    ph_stack = 180/pi.*atan2(imag(chz_stack),real(chz_stack));
    % Admittance
    ad_stack = abs(chz_stack)./chh_stack;
    
    
    [fmin,idx_flo] = min(abs(f-tiltfreq(1)));
    [fmin,idx_fhi] = min(abs(f-tiltfreq(2)));
    av_coh(ih) = mean(coh_stack(idx_flo:idx_fhi));
    av_ph(ih) = mean(abs(ph_stack(idx_flo:idx_fhi)));
end

[max_coh_fine,idx]=max(av_coh,[],2);
max_or_fine = hangs2(idx);

max_coh_fine = max_coh_fine;
max_or_fine = max_or_fine;

if isfig ==1
figure(90)
    subplot(211); hold on
    plot(days,max_or_fine,'o','MarkerFaceColor',cc,'MarkerSize',5,'MarkerEdgeColor','none');
    title('Orientation of Maximum Coherence'); xlabel('Days Since January 1'); ylabel('Degrees from H1'); ylim([0 360])
    subplot(212); hold on
    plot(days,max_coh_fine,'o','MarkerFaceColor',cc,'MarkerSize',5,'MarkerEdgeColor','none');
    title('Value of Maximum Coherence'); xlabel('Days Since January 1'); ylabel('Coherence'); ylim([0 1])
end

return