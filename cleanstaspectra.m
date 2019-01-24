% cleanstasspectra
% clean noise spectra to remove spurious days and calculate deployment
% averages
% Helen Janiszewski 1/2017

clear; close all


% CODE OPTIONS 
isfigure = 1;
issavefigure = 1;

% DO NOT EDIT BELOW
setup_parameter;
netsta = [network,station];
close all
inpath = sprintf('%s/SPECTRA/%s/',OUTdir,netsta);
figoutpath=sprintf('%s/STATIONS_NOISEPROP/',FIGdir);
outpath = sprintf('%s/AVG_STA/',OUTdir);

if ~exist(figoutpath)
    mkdir(figoutpath);
end

if ~exist(outpath)
    mkdir(outpath);
end

spectra_filenames = dir(fullfile(inpath,['*.mat']));

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
    
    % Station cross spectra
    c12_all(ie,:) = specprop.cross.c12_stack;
    c1p_all(ie,:) = specprop.cross.c1p_stack;
    c1z_all(ie,:) = specprop.cross.c1z_stack;
    c2p_all(ie,:) = specprop.cross.c2p_stack;
    c2z_all(ie,:) = specprop.cross.c2z_stack;
    cpz_all(ie,:) = specprop.cross.cpz_stack;
    
    spect(:,ie,1) = specprop.power.czz_stack;
    spect(:,ie,2) = specprop.power.c11_stack;
    spect(:,ie,3) = specprop.power.c22_stack;
    spect(:,ie,4) = specprop.power.cpp_stack;
    
    spectsm(:,ie,1) = smooth(specprop.power.czz_stack,40);
    spectsm(:,ie,2) = smooth(specprop.power.c11_stack,40);
    spectsm(:,ie,3) = smooth(specprop.power.c22_stack,40);
    spectsm(:,ie,4) = smooth(specprop.power.cpp_stack,40);
    
    % Coherence
    coh_stack(:,ie,1) = smooth(abs(specprop.cross.c1z_stack).^2./(specprop.power.c11_stack.*specprop.power.czz_stack),40);
    coh_stack(:,ie,2) = smooth(abs(specprop.cross.c2z_stack).^2./(specprop.power.c22_stack.*specprop.power.czz_stack),40);
    coh_stack(:,ie,3) = smooth(abs(specprop.cross.cpz_stack).^2./(specprop.power.cpp_stack.*specprop.power.czz_stack),40);
    coh_stack(:,ie,4) = smooth(abs(specprop.cross.c12_stack).^2./(specprop.power.c11_stack.*specprop.power.c22_stack),40);
    coh_stack(:,ie,5) = smooth(abs(specprop.cross.c1p_stack).^2./(specprop.power.c11_stack.*specprop.power.cpp_stack),40);
    coh_stack(:,ie,6) = smooth(abs(specprop.cross.c2p_stack).^2./(specprop.power.c22_stack.*specprop.power.cpp_stack),40);
    
    % Phase
    ph_stack(:,ie,1) = 180/pi.*atan2(imag(specprop.cross.c1z_stack),real(specprop.cross.c1z_stack));
    ph_stack(:,ie,2) = 180/pi.*atan2(imag(specprop.cross.c2z_stack),real(specprop.cross.c2z_stack));
    ph_stack(:,ie,3) = 180/pi.*atan2(imag(specprop.cross.cpz_stack),real(specprop.cross.cpz_stack));
    ph_stack(:,ie,4) = 180/pi.*atan2(imag(specprop.cross.c12_stack),real(specprop.cross.c12_stack));
    ph_stack(:,ie,5) = 180/pi.*atan2(imag(specprop.cross.c1p_stack),real(specprop.cross.c1p_stack));
    ph_stack(:,ie,6) = 180/pi.*atan2(imag(specprop.cross.c1p_stack),real(specprop.cross.c2p_stack));
    
    % Admittance
    ad_stack(:,ie,1) = smooth(abs(specprop.cross.c1z_stack)./specprop.power.c11_stack,40);
    ad_stack(:,ie,2) = smooth(abs(specprop.cross.c2z_stack)./specprop.power.c22_stack,40);
    ad_stack(:,ie,3) = smooth(abs(specprop.cross.cpz_stack)./specprop.power.cpp_stack,40);
    ad_stack(:,ie,4) = smooth(abs(specprop.cross.c12_stack)./specprop.power.c11_stack,40);
    ad_stack(:,ie,5) = smooth(abs(specprop.cross.c1p_stack)./specprop.power.c11_stack,40);
    ad_stack(:,ie,6) = smooth(abs(specprop.cross.c2p_stack)./specprop.power.c22_stack,40);
end

gooddays = QC_cleanstaspectra_days(spect,comp_exist,f,pb_dep,tolerance_dep,a_val_dep);

station = specprop.params.station;
network = specprop.params.network;

NWINS = repmat(nwins(gooddays)',[1,length(f)]);
ALLWINS = sum(nwins(gooddays));

czz_all = (czz_all(gooddays,:).*NWINS);
c11_all = (c11_all(gooddays,:).*NWINS);
c22_all = (c22_all(gooddays,:).*NWINS);
cpp_all = (cpp_all(gooddays,:).*NWINS);

c12_all = (c12_all(gooddays,:).*NWINS);
c1p_all = (c1p_all(gooddays,:).*NWINS);
c1z_all = (c1z_all(gooddays,:).*NWINS);
c2p_all = (c2p_all(gooddays,:).*NWINS);
c2z_all = (c2z_all(gooddays,:).*NWINS);
cpz_all = (cpz_all(gooddays,:).*NWINS);

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
    
    if ~isempty(find(gooddays == ie))
        cc(ie,:) = [0.5 0.5 0.5];
        specprop.params.badflag = 0;
    else
        cc(ie,:) = [1 0 0 ];
        disp(sprintf('Killed day %s', dayid));
        specprop.params.badflag = 1;
    end
    
    %     % Look for tilt orientations
    %     [specprop.directions.tilt,specprop.directions.tiltcoh,specprop.directions.micro,specprop.directions.microcoh]...
    %         = spec_orient(specprop,tiltfreq,microfreq,hangs,dayid,spectrum_Z,cspectrum_Z,spectrum_H1,spectrum_H2);
    %
    %     % Add badflag and orientations to file
%         filename = [inpath,station,'_',dayid,'_spectra.mat'];
%         save(filename,'specprop');
end

if isfigure
    plot_cleanstaspectra(spectsm,coh_stack,ph_stack,ad_stack,cc,netsta,f,maxpow,minpow);
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

staavg.power.czz_std = std(czz_all,1).*length(gooddays)/ALLWINS;
staavg.power.cpp_std = std(cpp_all,1).*length(gooddays)/ALLWINS;
staavg.power.c11_std = std(c11_all,1).*length(gooddays)/ALLWINS;
staavg.power.c22_std = std(c22_all,1).*length(gooddays)/ALLWINS;

% Calculate average cross spectra
staavg.cross.c12_mean = mean(c12_all,1).*length(gooddays)/ALLWINS;
staavg.cross.c1p_mean = mean(c1p_all,1).*length(gooddays)/ALLWINS;
staavg.cross.c1z_mean = mean(c1z_all,1).*length(gooddays)/ALLWINS;
staavg.cross.c2p_mean = mean(c2p_all,1).*length(gooddays)/ALLWINS;
staavg.cross.c2z_mean = mean(c2z_all,1).*length(gooddays)/ALLWINS;
staavg.cross.cpz_mean = mean(cpz_all,1).*length(gooddays)/ALLWINS;

staavg.cross.c12_std = std(c12_all,1).*length(gooddays)/ALLWINS;
staavg.cross.c1p_std = std(c1p_all,1).*length(gooddays)/ALLWINS;
staavg.cross.c1z_std = std(c1z_all,1).*length(gooddays)/ALLWINS;
staavg.cross.c2p_std = std(c2p_all,1).*length(gooddays)/ALLWINS;
staavg.cross.c2z_std = std(c2z_all,1).*length(gooddays)/ALLWINS;
staavg.cross.cpz_std = std(cpz_all,1).*length(gooddays)/ALLWINS;

avgprop.params.elev = elev;
avgprop.params.freqcomp = freqcomp;
avgprop.params.station = station;
avgprop.params.network = network;
avgprop.params.f = f;

filename = [outpath,netsta,'_spectraavg.mat'];
save(filename,'staavg','avgprop');



