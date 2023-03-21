% tilt_characteristic
%
% The tilt characteristic of the seismometer depends on the selection of the 
% tilt frequency band, which can be obtained from Figures 2-4 plotted by 
% b2_cleanstaspectra and changed in the setup_paremeter. In the tilt frequency 
% band, the coherence between Z and H is the largest, the admittance is 
% nearly constant and the phase shift nearly zero.
%
% Reference:
% Bell S. W., D. W. Forsyth, and Y. Ruan (2015). Removing noise from the
% vertical component records of ocean‐bottom seismometers: results from
% year one of the Cascadia Initiative. Bulletin of the Seismological
% Society of America, 105(1): 300–313. https://doi.org/10.1785/0120140054
%
% Yuechu Wu
% 12131066@mail.sustech.edu.cn
% 2022-10-08
%
% Added a figure with the number of days since deployment as the horizontal axis.
% Updated 2022-12-12, by Yuechu Wu


clear; close all;

isfigure = [1 1]; % FIGURE SWITCH %
% Figure 1: The horizontal axis is the number of days since deployment,
% which requires a relatively long time to plot.
% Figure 2: The horizontal axis is the date automatically matched,
% the beginning and ending few days may be trimmed.


% DO NOT EDIT BELOW
setup_parameter;
figoutpath = sprintf('%s/STATIONS_TILT',FIGdir);

if ~exist(figoutpath,'dir')
    mkdir(figoutpath);
end

for ista = 1:length(stations) % begin station loop
    station = stations{ista};
    netsta = [network,station];
    close all;
    inpath = sprintf('%s/SPECTRA/%s',OUTdir,netsta);

    spectra_filenames = dir(fullfile(inpath,'*.mat'));

    if isempty(spectra_filenames)
        continue
    end

    splitname_deploy = strsplit(spectra_filenames(1).name,'_');
    splitname_recovery = strsplit(spectra_filenames(end).name,'_');
    day_deploy   = splitname_deploy{2};
    day_recovery = splitname_recovery{2};
    deploynum    = datenum(day_deploy(1:8),'yyyymmdd');
    recoverynum  = datenum(day_recovery(1:8),'yyyymmdd');
    x = [deploynum:recoverynum]';

    tiltangs = NaN(length(x),1);
    tiltcohs = NaN(length(x),1);

    for ie = 1:length(spectra_filenames) % begin date loop

        splitname_ie = strsplit(spectra_filenames(ie).name,'_');
        dayid = splitname_ie{2};
        fprintf('%s\n',dayid);
        load(fullfile(inpath,spectra_filenames(ie).name));

        iday = datenum(dayid(1:8),'yyyymmdd') - deploynum + 1;

        f = specprop.params.f;
        npts_smooth = floor(specprop.params.NFFT/1000)+1;

        % Coherence between Z and H
        coh = smooth(abs(specprop.rotation.chz_stack).^2./(specprop.rotation.chh_stack.*specprop.power.czz_stack),npts_smooth);

        % Admittance between Z and H
        admit = smooth(abs(specprop.rotation.chz_stack)./specprop.rotation.chh_stack,npts_smooth);

        tiltf = find(f >= tiltfreq(1) & f <= tiltfreq(2));

        % Tilt angles are the arctangent of the admittance
        tiltang = atand(mean(admit(tiltf)));
        tiltcoh = mean(coh(tiltf));

        if isfigure(1)

            plot_tilt(dayid,day_deploy,tiltang,tiltcoh);

            if isleapyear(str2double(day_deploy(1:4)))
                year = 366;
            else
                year = 365;
            end
            
            figure(1)
            subplot(211)
            title(sprintf('%s %s',network,station));
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

            figurename1 = sprintf('%s/%s_tilt_days.pdf',figoutpath,netsta);
            print('-dpdf',figurename1);

        end % end figure 1

        tiltangs(iday,:) = tiltang;
        tiltcohs(iday,:) = tiltcoh;

    end % end date loop

    if isfigure(2)

        figure(2)

        subplot(2,1,1)
        plot(x,tiltangs,'o','MarkerFaceColor','#2c6db3','MarkerSize',5,'MarkerEdgeColor','none');
        xlim([deploynum recoverynum]);
        datetick('x','mmmdd','keepticks');
        ylabel('Tilt angle (°)');
        title(sprintf('%s %s',network,station));
        box on
        grid on

        subplot(2,1,2)
        plot(x,tiltcohs,'o','MarkerFaceColor','#2c6db3','MarkerSize',5,'MarkerEdgeColor','none');
        xlim([deploynum recoverynum]);
        ylim([0 1]);
        datetick('x','mmmdd','keepticks');
        ylabel('Coherence');
        xlabel(sprintf('Date of year %s to %s',day_deploy(1:4),day_recovery(1:4)));
        set(gca,'YTick',0:0.2:1);
        box on
        grid on

        figurename2 = sprintf('%s/%s_tilt_date.pdf',figoutpath,netsta);
        print('-dpdf',figurename2);

    end % end figure 2

end % end station loop
