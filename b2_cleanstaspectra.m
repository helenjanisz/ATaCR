% cleanstasspectra
% clean noise spectra to remove spurious days and calculate deployment
% averages
% Helen Janiszewski 1/2021
% 
% The unimportant parameters '12', '1P' and '2P' have been removed, and 'HZ' 
% has been added, which is useful for analyzing the tilt characteristics.
% Updated 2022-12-12, by Yuechu Wu
%
% P1 and P2 are added back in for use with assessing shallow water
% instrument behavior. 
% NOTE: could add user choice here in future
% Smoothing removed due to issues at long periods when using high sample
% rate data. Should replace with octave averaging.
% HAJ 01-04-2023
%
% The user choice of component options have been added.
% Updated 2023-01-25, by Yuechu Wu

clear; close all


% CODE OPTIONS
isfigure = 1;
issavefigure = 1;
isoverwrite = 1; % if set to 0, will skip previously processed files

% All component options, for plotting power spectra.
% comporder = {'Z','H1','H2','P'};
comporder = {'Z','H1','H2','P'};

% All component pair options, for plotting three parameters (coherence, 
% phase, and admittance) that describe the transfer function.
% plotorder = {'1Z','2Z','HZ','PZ','12','P1','P2'};
plotorder = {'1Z','2Z','HZ','PZ','P1','P2'};

% DO NOT EDIT BELOW
setup_parameter;
for ista = 1:length(stations) % begin station loop
    clear gooddays
    clear czz_all cpp_all c11_all c22_all chh_all
    clear c12_all c1p_all c1z_all c2p_all c2z_all cpz_all chz_all chp_all
    clear spect spectsm coh_stack ph_stack ad_stack orients oriecoh
    station = stations{ista};
    netsta = [network,station];
    close all
    inpath = sprintf('%s/SPECTRA/%s/',OUTdir,netsta);
    figoutpath=sprintf('%s/STATIONS_NOISEPROP/',FIGdir);
    outpath = sprintf('%s/AVG_STA/',OUTdir);

    if ~exist(figoutpath,'dir')
        mkdir(figoutpath);
    end

    if ~exist(outpath,'dir')
        mkdir(outpath);
    end

    if ~isoverwrite && exist([outpath,netsta,'_spectraavg.mat'],'file')
        display(['Skipping ',station]);
        continue
    end

    spectra_filenames = dir(fullfile(inpath,'*.mat'));
    if isempty(spectra_filenames) == 1
        continue
    end

    for ie = 1:length(spectra_filenames)
        load(fullfile(inpath,spectra_filenames(ie).name));
        a = length(specprop.params.station);
        f = specprop.params.f;
        comp_exist = specprop.params.comp_exist;
        dayid = spectra_filenames(ie).name((a+4):(a+15));
        disp(dayid);
        nwins(ie) = length(specprop.params.goodwins);

        % Station power spectra
        czz_all(ie,:) = specprop.power.czz_stack;
        cpp_all(ie,:) = specprop.power.cpp_stack;
        c11_all(ie,:) = specprop.power.c11_stack;
        c22_all(ie,:) = specprop.power.c22_stack;
        chh_all(ie,:) = specprop.rotation.chh_stack;

        % Station cross spectra
        c12_all(ie,:) = specprop.cross.c12_stack;
        c1p_all(ie,:) = specprop.cross.c1p_stack;
        c1z_all(ie,:) = specprop.cross.c1z_stack;
        c2p_all(ie,:) = specprop.cross.c2p_stack;
        c2z_all(ie,:) = specprop.cross.c2z_stack;
        cpz_all(ie,:) = specprop.cross.cpz_stack;
        chz_all(ie,:) = specprop.rotation.chz_stack;
        chp_all(ie,:) = specprop.rotation.chp_stack;

        spect(:,ie,1) = specprop.power.czz_stack;
        spect(:,ie,2) = specprop.power.c11_stack;
        spect(:,ie,3) = specprop.power.c22_stack;
        spect(:,ie,4) = specprop.power.cpp_stack;
        
        spectsm(:,ie,1) = (specprop.power.czz_stack);
        spectsm(:,ie,2) = (specprop.power.c11_stack);
        spectsm(:,ie,3) = (specprop.power.c22_stack);
        spectsm(:,ie,4) = (specprop.power.cpp_stack);

        % orientations
        orients(ie) = specprop.params.rotor;
        oriecoh(ie) = specprop.params.rotcoh;

        % Coherence
        coh_stack(:,ie,1) = (abs(specprop.cross.c1z_stack).^2./(specprop.power.c11_stack.*specprop.power.czz_stack));
        coh_stack(:,ie,2) = (abs(specprop.cross.c2z_stack).^2./(specprop.power.c22_stack.*specprop.power.czz_stack));
        coh_stack(:,ie,3) = (abs(specprop.rotation.chz_stack).^2./(specprop.rotation.chh_stack.*specprop.power.czz_stack));
        coh_stack(:,ie,4) = (abs(specprop.cross.cpz_stack).^2./(specprop.power.cpp_stack.*specprop.power.czz_stack));
        coh_stack(:,ie,5) = (abs(specprop.cross.c12_stack).^2./(specprop.power.c11_stack.*specprop.power.c22_stack));
        coh_stack(:,ie,6) = (abs(specprop.cross.c1p_stack).^2./(specprop.power.c11_stack.*specprop.power.cpp_stack));
        coh_stack(:,ie,7) = (abs(specprop.cross.c2p_stack).^2./(specprop.power.c22_stack.*specprop.power.cpp_stack));

        % Phase
        ph_stack(:,ie,1) = 180/pi.*atan2(imag(specprop.cross.c1z_stack),real(specprop.cross.c1z_stack));
        ph_stack(:,ie,2) = 180/pi.*atan2(imag(specprop.cross.c2z_stack),real(specprop.cross.c2z_stack));
        ph_stack(:,ie,3) = 180/pi.*atan2(imag(specprop.rotation.chz_stack),real(specprop.rotation.chz_stack));
        ph_stack(:,ie,4) = 180/pi.*atan2(imag(specprop.cross.cpz_stack),real(specprop.cross.cpz_stack));
        ph_stack(:,ie,5) = 180/pi.*atan2(imag(specprop.cross.c12_stack),real(specprop.cross.c12_stack));
        ph_stack(:,ie,6) = 180/pi.*atan2(imag(specprop.cross.c1p_stack),real(specprop.cross.c1p_stack));
        ph_stack(:,ie,7) = 180/pi.*atan2(imag(specprop.cross.c1p_stack),real(specprop.cross.c2p_stack));

        % Admittance
        ad_stack(:,ie,1) = (abs(specprop.cross.c1z_stack)./specprop.power.c11_stack);
        ad_stack(:,ie,2) = (abs(specprop.cross.c2z_stack)./specprop.power.c22_stack);
        ad_stack(:,ie,3) = (abs(specprop.rotation.chz_stack)./specprop.rotation.chh_stack);
        ad_stack(:,ie,4) = (abs(specprop.cross.cpz_stack)./specprop.power.cpp_stack);
        ad_stack(:,ie,5) = (abs(specprop.cross.c12_stack)./specprop.power.c11_stack);
        ad_stack(:,ie,6) = (abs(specprop.cross.c1p_stack)./specprop.power.c11_stack);
        ad_stack(:,ie,7) = (abs(specprop.cross.c2p_stack)./specprop.power.c22_stack);
    end

    gooddays = QC_cleanstaspectra_days(spect,comp_exist,f,pb_dep,tolerance_dep,a_val_dep);

    station = specprop.params.station;
    network = specprop.params.network;

    NWINS = repmat(nwins(gooddays)',[1,length(f)]);
    ALLWINS = sum(nwins(gooddays));

    czz_all_st = czz_all(gooddays,:);
    c11_all_st = c11_all(gooddays,:);
    c22_all_st = c22_all(gooddays,:);
    cpp_all_st = cpp_all(gooddays,:);
    chh_all_st = chh_all(gooddays,:);

    czz_all = (czz_all(gooddays,:).*NWINS);
    c11_all = (c11_all(gooddays,:).*NWINS);
    c22_all = (c22_all(gooddays,:).*NWINS);
    cpp_all = (cpp_all(gooddays,:).*NWINS);
    chh_all = (chh_all(gooddays,:).*NWINS);

    c12_all = (c12_all(gooddays,:).*NWINS);
    c1p_all = (c1p_all(gooddays,:).*NWINS);
    c1z_all = (c1z_all(gooddays,:).*NWINS);
    c2p_all = (c2p_all(gooddays,:).*NWINS);
    c2z_all = (c2z_all(gooddays,:).*NWINS);
    cpz_all = (cpz_all(gooddays,:).*NWINS);
    chz_all = (chz_all(gooddays,:).*NWINS);
    chp_all = (chp_all(gooddays,:).*NWINS);

    orall = orients(gooddays);
    orcohall = oriecoh(gooddays);

    maxpow(1) = max(max(czz_all(:,1:end-1)))*100;
    maxpow(2) = max(max(c11_all(:,1:end-1)))*100;
    maxpow(3) = max(max(c22_all(:,1:end-1)))*100;
    maxpow(4) = max(max(cpp_all(:,1:end-1)))*100;

    minpow(1) = min(min(czz_all(:,1:end-1)))/10;
    minpow(2) = min(min(c11_all(:,1:end-1)))/10;
    minpow(3) = min(min(c22_all(:,1:end-1)))/10;
    minpow(4) = min(min(cpp_all(:,1:end-1)))/10;

    for  ie = 1:length(spectra_filenames)
        load(fullfile(inpath,spectra_filenames(ie).name));
        dayid = spectra_filenames(ie).name((a+4):(a+15));
        elev = specprop.params.elev;
        freqcomp = sqrt(9.8/(2*pi*elev));

        if ~isempty(find(gooddays == ie, 1))
            cc(ie,:) = [0.5 0.5 0.5];
            specprop.params.badflag = 0;
            daysused(ie).id = dayid;
        else
            cc(ie,:) = [1 0 0];
            fprintf('Killed day %s\n', dayid);
            specprop.params.badflag = 1;
            daysused(ie).id = NaN;
        end

        % Add badflag to file
        filename = [inpath,netsta,'_',dayid,'_spectra.mat'];
        save(filename,'specprop');
    end

    if isfigure
        plot_cleanstaspectra(spectsm,coh_stack,ph_stack,ad_stack,cc,netsta,f,maxpow,minpow,comporder,plotorder);
    end

    if isfigure && issavefigure
        figure(1)
        filename=sprintf('%s/%s_spectra',figoutpath,netsta);
        print(gcf,'-dpng',filename)
        figure(2)
        filename=sprintf('%s/%s_coherence',figoutpath,netsta);
        print(gcf,'-dpng',filename)
        figure(3)
        filename=sprintf('%s/%s_phase',figoutpath,netsta);
        print(gcf,'-dpng',filename)
        figure(4)
        filename=sprintf('%s/%s_admittance',figoutpath,netsta);
        print(gcf,'-dpng',filename)
    end

    % Calculate average spectra
    staavg.power.czz_mean = mean(czz_all,1).*length(gooddays)/ALLWINS;
    staavg.power.cpp_mean = mean(cpp_all,1).*length(gooddays)/ALLWINS;
    staavg.power.c11_mean = mean(c11_all,1).*length(gooddays)/ALLWINS;
    staavg.power.c22_mean = mean(c22_all,1).*length(gooddays)/ALLWINS;
    staavg.rotation.chh_mean = mean(chh_all,1).*length(gooddays)/ALLWINS;

    staavg.power.czz_std = std(czz_all_st,1);
    staavg.power.cpp_std = std(cpp_all_st,1);
    staavg.power.c11_std = std(c11_all_st,1);
    staavg.power.c22_std = std(c22_all_st,1);

    % Calculate average cross spectra
    staavg.cross.c12_mean = mean(c12_all,1).*length(gooddays)/ALLWINS;
    staavg.cross.c1p_mean = mean(c1p_all,1).*length(gooddays)/ALLWINS;
    staavg.cross.c1z_mean = mean(c1z_all,1).*length(gooddays)/ALLWINS;
    staavg.cross.c2p_mean = mean(c2p_all,1).*length(gooddays)/ALLWINS;
    staavg.cross.c2z_mean = mean(c2z_all,1).*length(gooddays)/ALLWINS;
    staavg.cross.cpz_mean = mean(cpz_all,1).*length(gooddays)/ALLWINS;
    staavg.rotation.chz_mean = mean(chz_all,1).*length(gooddays)/ALLWINS;
    staavg.rotation.chp_mean = mean(chp_all,1).*length(gooddays)/ALLWINS;

    avgprop.params.elev = elev;
    avgprop.params.freqcomp = freqcomp;
    avgprop.params.station = station;
    avgprop.params.network = network;
    avgprop.params.f = f;
    avgprop.params.ndays = length(gooddays);
    avgprop.params.nwindows = ALLWINS;
    avgprop.params.orientation = mean(orall);
    avgprop.params.oriencohere = mean(orcohall);
    avgprop.params.orientationstd = std(orall);
    avgprop.params.oriencoherestd = std(orcohall);
    avgprop.params.daysused = daysused;

    filename = [outpath,netsta,'_spectraavg.mat'];
    save(filename,'staavg','avgprop');
end % end station loop
