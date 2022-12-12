function [specprop] =noisecal_dailystaspectra(spectrum_Z,spectrum_H1,spectrum_H2,spectrum_P,...
    cspectrum_Z,cspectrum_H1,cspectrum_H2,cspectrum_P,is_goodwin,f,comp_exist,linecol,...
    taxisZ,iptwin1,iptwin2,NFFT,dt,isfigure_powerspec,Zraw,dayid,day_deploy,hangs,tiltfreq,calrotation,isfigure_orient);

smoothlevel = floor(NFFT/1000)+1;

c11_stack=zeros(1,length(f))';
c22_stack=zeros(1,length(f))';
cpp_stack=zeros(1,length(f))';
czz_stack=zeros(1,length(f))';
c1z_stack=zeros(1,length(f))';
c2z_stack=zeros(1,length(f))';
cpz_stack=zeros(1,length(f))';
c12_stack=zeros(1,length(f))';
c1p_stack=zeros(1,length(f))';
c2p_stack=zeros(1,length(f))';


nwin_stack = 0;

for iwin = 1:length(is_goodwin)
    pts_begin = iptwin1(iwin);
    pts_end = iptwin2(iwin);
    % power spectrum for each segment
    if comp_exist(1)==1
        czz = abs(spectrum_Z(iwin,:)').^2*2/(NFFT*dt);
    else czz = nan(1,length(f))';
    end
    if comp_exist(2)==1
        c11 = abs(spectrum_H1(iwin,:)').^2*2/(NFFT*dt);
    else c11 =nan(1,length(f))';
    end
    if comp_exist(3)==1
        c22 = abs(spectrum_H2(iwin,:)').^2*2/(NFFT*dt);
    else c22= nan(1,length(f))';
    end
    if comp_exist(4)==1
        cpp = abs(spectrum_P(iwin,:)').^2*2/(NFFT*dt);
    else cpp = nan(1,length(f))';
    end

    % complex cross-spectrum for each segment
    if comp_exist(1)==1 && comp_exist(2)==1
        c1z = spectrum_H1(iwin,:)'.*cspectrum_Z(iwin,:)'*2/(NFFT*dt);
    end
    if comp_exist(1)==1 && comp_exist(3)==1
        c2z = spectrum_H2(iwin,:)'.*cspectrum_Z(iwin,:)'*2/(NFFT*dt);
    end
    if comp_exist(1)==1 && comp_exist(4)==1
        cpz = spectrum_P(iwin,:)'.*cspectrum_Z(iwin,:)'*2/(NFFT*dt);
    end
    if comp_exist(2)==1 && comp_exist(3)==1
        c12 = spectrum_H1(iwin,:)'.*cspectrum_H2(iwin,:)'*2/(NFFT*dt);
    end
    if comp_exist(2)==1 && comp_exist(4)==1
        c1p = spectrum_H1(iwin,:)'.*cspectrum_P(iwin,:)'*2/(NFFT*dt);
    end
    if comp_exist(3)==1 && comp_exist(4)==1
        c2p = spectrum_H2(iwin,:)'.*cspectrum_P(iwin,:)'*2/(NFFT*dt);
    end

    if is_goodwin(iwin) == 1;
        if comp_exist(1) == 1
            czz_stack=czz_stack+czz;
        end
        if comp_exist(2) ==1
            c11_stack=c11_stack+c11;
        end
        if comp_exist(3) ==1
            c22_stack=c22_stack+c22;
        end
        if comp_exist(4) ==1
            cpp_stack=cpp_stack+cpp;
        end

        if comp_exist(1)==1 && comp_exist(2)==1
            c1z_stack=c1z_stack+c1z;
        end
        if comp_exist(1)==1 && comp_exist(3)==1
            c2z_stack=c2z_stack+c2z;
        end
        if comp_exist(1)==1 && comp_exist(4)==1
            cpz_stack=cpz_stack+cpz;
        end
        if comp_exist(2)==1 && comp_exist(3)==1
            c12_stack=c12_stack+c12;
        end
        if comp_exist(2)==1 && comp_exist(4)==1
            c1p_stack=c1p_stack+c1p;
        end
        if comp_exist(3)==1 && comp_exist(4)==1
            c2p_stack=c2p_stack+c2p;
        end


        nwin_stack=nwin_stack+1;

        if isfigure_powerspec
            figure(1)
            subplot(421)
            loglog(f,(smooth(czz,smoothlevel)),linecol,'LineWidth',.5); hold on
            xlim([1/250,max(f)]);
            title('Z, Windows')
            ylabel('PSD (m^2/Hz)')
            subplot(423)
            loglog(f,(smooth(c11,smoothlevel)),linecol,'LineWidth',.5); hold on
            xlim([1/250,max(f)]);
            title('H1, Windows')
            ylabel('PSD (m^2/Hz)')
            subplot(425)
            loglog(f,(smooth(c22,smoothlevel)),linecol,'LineWidth',.5); hold on
            xlim([1/250,max(f)]);
            title('H2, Windows')
            ylabel('PSD (m^2/Hz)')
            subplot(427)
            loglog(f,(smooth(cpp,smoothlevel)),linecol,'LineWidth',.5); hold on
            xlim([1/250,max(f)]);
            title('P, Windows')
            xlabel('Frequency (Hz)')
            ylabel('PSD (m^2/Hz)')
        end
    end
end

if calrotation==1
    % get orientation
    [max_coh,max_or] = spec_orient(spectrum_Z,spectrum_H1,spectrum_H2,cspectrum_Z,hangs,tiltfreq,f,isfigure_orient,dayid,day_deploy,is_goodwin,NFFT,dt);
    % calculate needed spectral parameters for TFs down the line for
    % maximum orientation
    hang = max_or;
    cang = cos(hang*pi/180);
    sang = sin(hang*pi/180);

    chh_stack=zeros(1,length(f))';
    czz_stack_rot=zeros(1,length(f))';
    chz_stack=zeros(1,length(f))';
    chp_stack=zeros(1,length(f))';
    nwin_stack_rot = 0;
    for ista = 1:length(is_goodwin)
        if is_goodwin(ista)==0
            continue
        end
        spec_H = sang.*spectrum_H2(ista,:)'+cang.*spectrum_H1(ista,:)'; %rotated horizontal spectra and cross-spectra

        chh = abs(spec_H).^2*2/(NFFT*dt);
        chz = spec_H.*(cspectrum_Z(ista,:)')*2/(NFFT*dt);
        if comp_exist(4)==1
        chp = spec_H.*(cspectrum_P(ista,:)')*2/(NFFT*dt);
        else chp = nan(1,length(f))';
        end

        chh_stack=chh_stack+chh;
        chz_stack=chz_stack+chz;
        chp_stack=chp_stack+chp;
        nwin_stack_rot = nwin_stack_rot +1;
    end

    chh_stack = chh_stack/nwin_stack_rot;
    chz_stack = chz_stack/nwin_stack_rot;
    chp_stack = chp_stack/nwin_stack_rot;
end

specprop.power.c11_stack=c11_stack/nwin_stack;
specprop.power.c22_stack=c22_stack/nwin_stack;
specprop.power.cpp_stack=cpp_stack/nwin_stack;
specprop.power.czz_stack=czz_stack/nwin_stack;

if calrotation==1
    specprop.rotation.chh_stack = chh_stack;
    specprop.rotation.chz_stack = chz_stack;
    specprop.rotation.chp_stack = chp_stack;

    specprop.params.rotor = max_or;
    specprop.params.rotcoh = max_coh;
end

specprop.cross.c1z_stack=c1z_stack/nwin_stack;
specprop.cross.c2z_stack=c2z_stack/nwin_stack;
specprop.cross.cpz_stack=cpz_stack/nwin_stack;
specprop.cross.c12_stack=c12_stack/nwin_stack;
specprop.cross.c1p_stack=c1p_stack/nwin_stack;
specprop.cross.c2p_stack=c2p_stack/nwin_stack;

if isfigure_powerspec
    figure(1)
    subplot(422)
    loglog(f,(smooth(specprop.power.czz_stack,smoothlevel)),linecol,'LineWidth',1); hold on
    xlim([1/250,max(f)]);
    title('Z, Daily Average')
    subplot(424)
    loglog(f,(smooth(specprop.power.c11_stack,smoothlevel)),linecol,'LineWidth',1); hold on
    xlim([1/250,max(f)]);
    title('H1, Daily Average')
    subplot(426)
    loglog(f,(smooth(specprop.power.c22_stack,smoothlevel)),linecol,'LineWidth',1); hold on
    xlim([1/250,max(f)]);
    title('H2, Daily Average')
    subplot(428)
    loglog(f,(smooth(specprop.power.cpp_stack,smoothlevel)),linecol,'LineWidth',1); hold on
    xlim([1/250,max(f)]);
    title('P, Daily Average')
    xlabel('Frequency (Hz)')
end

return
