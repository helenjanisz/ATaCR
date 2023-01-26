% dailystaspectra
% calculate noise properties for daily records at a given station
% at this point all pre-processing steps for station should be completed
% see download_data and daydata_preprocess for codes to download and
% pre-process or use your own
% Helen Janiszewski 1/2017

% updated from caxis to clim for compatability with matlab 2022. Change
% back to caxis lines 253,262,276,283 if using older matlab.

% haj 01/2023

clear; close all

% CODE OPTIONS
isfigure_spectrogram = 1; % generate spectrograms for each day
isfigure_powerspec = 1; % generate power spectra for each day
isfigure_orient = 1; % generate tilt orientation for each day
issavefigure = 1; % save output figures
isoverwrite = 1; % overwrite spectra files

% DO NOT EDIT BELOW
setup_parameter;

for ista = 1:length(stations)
    station = stations{ista};
    netsta = [network,station];
    
    outpath = sprintf('%s/SPECTRA/%s/',OUTdir,netsta);
    figoutpath=sprintf('%s/STATIONS_DAILY/',FIGdir);
    rotoutpath=sprintf('%s/STATIONS_ROTATION/',FIGdir);

    c = colormap('jet');

    if ~exist(outpath,'dir')
        mkdir(outpath);
    end

    if ~exist(figoutpath,'dir')
        mkdir(figoutpath);
    end

    if ~exist(rotoutpath,'dir')
        mkdir(rotoutpath);
    end

    disp(netsta);
    filesuff = sprintf('*_%s_%s.mat',network,station);
    day_filenames = dir(fullfile(WORKINGdir,network,station,'/',filesuff));

    hangs = [0:10:360]; 
    if isfigure_orient==1
        figure(90);clf
        figure(105);clf
    end
    if ~isempty(day_filenames)
        day_deploy = day_filenames(1).name(1:12);
    end
    for ie = 1 : length(day_filenames)
        if isfigure_spectrogram == 1
        figure(96);clf
        end
        if isfigure_powerspec==1
        figure(1);clf
        end
        day_data = load(fullfile(WORKINGdir,network, station ,'/',day_filenames(ie).name));
        traces_day = day_data.traces_day;
        dayid = day_filenames(ie).name(1:12);
        specfilename = sprintf('%s%s_%s_spectra.mat',outpath,netsta,dayid);
        if exist(specfilename,'file') && isoverwrite ==0
            disp(['Exist: ',specfilename,', Skip!']);
            continue;
        end
        disp(dayid);
        
        idxZ = find(ismember({traces_day.channel},chz_vec));
        idx1 = find(ismember({traces_day.channel},ch1_vec));
        idx2 = find(ismember({traces_day.channel},ch2_vec));
        idxP = find(ismember({traces_day.channel},chp_vec));
        
        if length(idxZ)>1 || length(idx1)>1 || length(idx2)>1 || length(idxP)>1
            disp('Skipping. Too many records for single channel.')
            continue
        end
        
        comp_exist = [~isempty(idxZ), ~isempty(idx1), ~isempty(idx2), ~isempty(idxP)];
        if ~isempty(find(comp_exist==0, 1))
            disp('Warning. At least one component is missing.')
            if allchans ==1
            continue
            end
        end
        [Zraw,elevZ,rateZ,dtZ,startZ,endZ] = varsetup_dailystaspectra(traces_day(idxZ),comp_exist(1));
        [H1raw,elev1,rate1,dt1,start1,end1] = varsetup_dailystaspectra(traces_day(idx1),comp_exist(2));
        [H2raw,elev2,rate2,dt2,start2,end2] = varsetup_dailystaspectra(traces_day(idx2),comp_exist(3));
        [Praw,elevP,rateP,dtP,startP,endP] = varsetup_dailystaspectra(traces_day(idxP),comp_exist(4));
        
        elev = [elevZ,elev1,elev2,elevP];
        elev = unique(elev(comp_exist==1));
        if length(elev)>1
            disp('Warning, mismatched elevations')
        end
        freqcomp = sqrt(9.8/(2*pi*elev)); %infragravity wave limit (theoretical estimate)
        
        starttimes = [startZ,start1,start2,startP];
        if length(unique(starttimes))>1
            if diff(unique(starttimes))>dtZ
                disp('Skipping. Not all start times are equal.')
                continue
            end
        end
        
        % check all samp rates are equal, if not resample at highest common
        % denominator
        samprate = [rateZ,rate1,rate2,rateP]; 
        if length(unique(samprate))>1
            [rates,is,ir] = unique(samprate);
            newrate = rates(1);
            for ii=2:length(rates)
                newrate = gcd(newrate,rates(ii));
            end
            dt = 1/newrate;
            [Zraw,taxisZ] = resample(Zraw,dt,dtZ);
            [H1raw,taxis1] = resample(H1raw,dt,dt1);
            [H2raw,taxis2] = resample(H2raw,dt,dt2);
            [Praw,taxisP] = resample(Praw,dt,dtP);
        else
            dt = 1/rateZ;
            taxisZ = [0:dt:(length(Zraw)-1)*dt]';
            taxis1 = [0:dt:(length(H1raw)-1)*dt]';
            taxis2 = [0:dt:(length(H2raw)-1)*dt]';
            taxisP = [0:dt:(length(Praw)-1)*dt]';
        end
        
        data_length = [length(Zraw),length(H1raw),length(H2raw),length(Praw)];
        if length(unique(data_length(comp_exist==1)))>1
            disp('Skipping. Different data lengths.')
            continue
        end
        
        % Check for all zeros
        data_zeros = [sum(abs(Zraw)), sum(abs(H1raw)), sum(abs(H2raw)), sum(abs(Praw))];
        ind_zero = data_zeros == 0;
        if sum(ind_zero) > 0
            if ind_zero(1) == 1
                disp('Z all zeros... skipping!')
            end
            if ind_zero(2) == 1
                disp('H1 all zeros... skipping!')
            end
            if ind_zero(3) == 1
                disp('H2 all zeros... skipping!')
            end
            if ind_zero(4) == 1
                disp('P all zeros... skipping!')
            end
            if allchans ==1
            continue
            end
        end
        % Spectra parameters for resampled data
        npts = T/dt;
        Ppower = nextpow2(npts);
        NFFT = 2^Ppower;
        npad0 = floor((NFFT-npts)/2);
        samprate = 1/dt; 
        f = samprate/2*linspace(0,1,NFFT/2+1);
        overn = floor(overlap*npts);
        
        if data_length/samprate<min_data_length
            disp('Skipping. Not enough data for the day.')
            continue
        end
        
        % Quality control to identify outlier windows and save spectrograms
        [goodwins,iptwin1,iptwin2] = QC_dailystaspectra_windows([Zraw,H1raw,H2raw,Praw],comp_exist,taxisZ,npts,overlap,T,isfigure_spectrogram,pb,tolerance,a_val);
        if isfigure_spectrogram && issavefigure
        figure(96)
        set(gcf,'PaperPositionMode','manual');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'PaperPosition',[.05 .05 8 10.5]);
        filename=sprintf('%s/%s_%s_spectrogram',figoutpath,netsta,dayid);
        print(gcf,'-dpng',filename)
        end
        
        % check if there are enough good windows to calculate day spectra
        nwin = length(goodwins);
        if nwin<minwin
            disp('Not enough good data segments.')
            continue
        else
            fprintf('%d good windows.Proceeding...\n',nwin);
        end
        
        % Calculating spectra for each component
        [spectrum_Z,cspectrum_Z,is_goodwin_Z] = speccal_dailystaspectra(nwin,f,iptwin1,iptwin2,goodwins,Zraw,comp_exist(1),taxisZ,overlap,dt,npts,npad0,NFFT);
        [spectrum_H1,cspectrum_H1,is_goodwin_H1] = speccal_dailystaspectra(nwin,f,iptwin1,iptwin2,goodwins,H1raw,comp_exist(2),taxis1,overlap,dt,npts,npad0,NFFT);
        [spectrum_H2,cspectrum_H2,is_goodwin_H2] = speccal_dailystaspectra(nwin,f,iptwin1,iptwin2,goodwins,H2raw,comp_exist(3),taxis2,overlap,dt,npts,npad0,NFFT);
        [spectrum_P,cspectrum_P,is_goodwin_P] = speccal_dailystaspectra(nwin,f,iptwin1,iptwin2,goodwins,Praw,comp_exist(4),taxisP,overlap,dt,npts,npad0,NFFT);
     
        is_goodwin_vec = vertcat(is_goodwin_Z,is_goodwin_H1,is_goodwin_H2,is_goodwin_P);
        first_good = find(comp_exist==1);
        first_good = first_good(1);
        is_goodwin = is_goodwin_vec(first_good(1),:);

        if isfigure_powerspec
        [specprop_all] = noisecal_dailystaspectra(spectrum_Z,spectrum_H1,spectrum_H2,spectrum_P,...
            cspectrum_Z,cspectrum_H1,cspectrum_H2,cspectrum_P,ones(size(is_goodwin)),f,comp_exist,...
            '-r',taxisZ,iptwin1,iptwin2,NFFT,dt,isfigure_powerspec,Zraw,dayid,day_deploy,hangs,tiltfreq,0,isfigure_orient);
        end
        [specprop] = noisecal_dailystaspectra(spectrum_Z,spectrum_H1,spectrum_H2,spectrum_P,...
            cspectrum_Z,cspectrum_H1,cspectrum_H2,cspectrum_P,is_goodwin,f,comp_exist,...
            '-k',taxisZ,iptwin1,iptwin2,NFFT,dt,isfigure_powerspec,Zraw,dayid,day_deploy,hangs,tiltfreq,1,isfigure_orient);
        
        % Save parameters in structure
        specprop.params.f = f;
        specprop.params.station = station;
        specprop.params.network = network;
        specprop.params.elev = elev;
        specprop.params.iptwin1 = iptwin1;
        specprop.params.iptwin2 = iptwin2;
        specprop.params.goodwins = goodwins;
        specprop.params.taxis = taxisZ;
        specprop.params.NFFT = NFFT;
        specprop.params.dt = dt;
        specprop.params.overlap = overlap;
        specprop.params.comp_exist = comp_exist;
        specprop.params.badflag = 0;
        
        filename = [outpath,netsta,'_',dayid,'_spectra.mat'];
        save(filename,'specprop');
        
        if isfigure_powerspec && issavefigure
        figure(1)
        set(gcf,'PaperPositionMode','manual');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'PaperPosition',[.05 .05 8 10.5]);
        
        figure(1)
        filename=sprintf('%s/%s_%s_dailyspectra',figoutpath,netsta,dayid);
        print(gcf,'-dpng',filename)
        end
        
        
        clear iptwin1 iptwin2
    end

    if isfigure_orient && issavefigure

        if isleapyear(str2double(day_deploy(1:4)))
            year = 366;
        else
            year = 365;
        end

        figure(105)
        subplot(121)
        clim([0 year])
        colormap jet
        colorbar
        ylabel(colorbar,'Day of year')
        subplot(122)
        clim([0 year])
        colormap jet
        colorbar
        ylabel(colorbar,'Day of year')
        set(gcf,'PaperPositionMode','manual');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'PaperPosition',[.05 .05 8 3]);
        
        figure(105)
        filename=sprintf('%s/%s_rotationcoherence',rotoutpath,netsta);
        print(gcf,'-dpng',filename)
        
        figure(90)
        subplot(211)
        clim([0 year])
        colormap jet
        colorbar
        ylabel(colorbar,'Day of year')
        box on
        grid on
        
        subplot(212)
        clim([0 year])
        colormap jet
        colorbar
        ylabel(colorbar,'Day of year')
        box on
        grid on
        set(gcf,'PaperPositionMode','manual');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'PaperPosition',[.05 .05 5 5]);
        
        figure(90)
        filename=sprintf('%s/%s_orientationdir',rotoutpath,netsta);
        print(gcf,'-dpng',filename)
    end
end % station loop