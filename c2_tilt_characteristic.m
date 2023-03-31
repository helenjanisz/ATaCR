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
% Added a figure with the number of days since deployment as the x-axis.
% 2022-12-12, Yuechu Wu
%
% Added write tilt information to spectra file. Added plot tilt 
% orientation, which is the same as Figure 90 of b1_dailystaspectra.
% Optimized running speed.
% 2023-02-09, Yuechu Wu



clear; close all;

% FIGURE SWITCH %
% Figures 1 & 3: The horizontal axis is the number of days since deployment.
% Figures 2 & 4: The horizontal axis is the date automatically matched,
% the beginning and ending few days may be trimmed.
isfigure_angle  = [1 1]; % TILT ANGLE, corresponding to Figures 1 & 2
isfigure_orient = [1 1]; % TILT ORIENTATION, corresponding to Figures 3 & 4

isoverwrite = 1; % overwrite tilt information to spectra file


% DO NOT EDIT BELOW
setup_parameter;
figoutpath = sprintf('%s/STATIONS_TILT',FIGdir);

if ~exist(figoutpath,'dir')
    mkdir(figoutpath);
end

for ista = 1:length(stations) % begin station loop
    station = stations{ista};
    netsta = [network,station];
    disp(netsta);
    close all;
    inpath = sprintf('%s/SPECTRA/%s',OUTdir,netsta);

    spectra_filenames = dir(fullfile(inpath,'*.mat'));

    if isempty(spectra_filenames)
        continue
    end

    splitname_deploy   = strsplit(spectra_filenames(1).name,'_');
    splitname_recovery = strsplit(spectra_filenames(end).name,'_');
    day_deploy   = splitname_deploy{2};
    day_recovery = splitname_recovery{2};
    deploynum    = datenum(day_deploy(1:8),'yyyymmdd');
    recoverynum  = datenum(day_recovery(1:8),'yyyymmdd');
    x = transpose(deploynum:recoverynum);

    deploydate = datestr(deploynum,'yyyy-mmm-dd');

    idays    = transpose(1:length(x));
    tiltangs = NaN(length(x),1);
    tiltcohs = NaN(length(x),1);
    rotors   = NaN(length(x),1);
    rotcohs  = NaN(length(x),1);

    for ie = 1:length(spectra_filenames) % begin date loop

        clear specprop
        splitname_ie = strsplit(spectra_filenames(ie).name,'_');
        dayid = splitname_ie{2};
        fprintf('%s\n',dayid);
        load(fullfile(inpath,spectra_filenames(ie).name));

        iday = datenum(dayid(1:8),'yyyymmdd') - deploynum + 1;

        rotor  = specprop.params.rotor;
        rotcoh = specprop.params.rotcoh;

        if isoverwrite || ~isfield(specprop,'tilt')

            f = specprop.params.f;
%             npts_smooth = floor(specprop.params.NFFT/1000)+1;

            % Coherence between H and Z
            coh = abs(specprop.rotation.chz_stack).^2./(specprop.rotation.chh_stack.*specprop.power.czz_stack);

            % Admittance between H and Z
            admit = abs(specprop.rotation.chz_stack)./specprop.rotation.chh_stack;

            tiltf = find(f >= tiltfreq(1) & f <= tiltfreq(2));

            % Tilt angles are the arctangent of the admittance
            tiltang = atand(mean(admit(tiltf)));
            tiltcoh = mean(coh(tiltf));

            specprop.tilt.ang = tiltang;
            specprop.tilt.coh = tiltcoh;

            % Add tilt information to file
            filename = sprintf('%s/%s_%s_spectra.mat',inpath,netsta,dayid);
            save(filename,'specprop');
        else
            tiltang = specprop.tilt.ang;
            tiltcoh = specprop.tilt.coh;
        end

        tiltangs(iday,:) = tiltang;
        tiltcohs(iday,:) = tiltcoh;
        rotors(iday,:)   = rotor;
        rotcohs(iday,:)  = rotcoh;

    end % end date loop

    % PLOTTING FIGURES
    if isfigure_angle(1)

        plot_tilt(idays,day_deploy,tiltangs,tiltcohs);

        if isleapyear(str2double(day_deploy(1:4)))
            year = 366;
        else
            year = 365;
        end

        figure(1) % tilt angle, x-axis -> days since deployment
        subplot(211)
        ylabel('Tilt angle (°)');
        title(sprintf('%s %s',network,station));
        clim([0 year]);
        colormap jet
        colorbar
        ylabel(colorbar,'Day of year');
        box on
        grid on

        subplot(212)
        ylim([0 1]);
        set(gca,'YTick',0:0.2:1);
        ylabel('Coherence');
        xlabel(sprintf('%s (%s)','Days since deployment',deploydate));
        clim([0 year]);
        colormap jet
        colorbar
        ylabel(colorbar,'Day of year');
        box on
        grid on

        figurename1 = sprintf('%s/%s_tilt_angle_days.pdf',figoutpath,netsta);
        print('-dpdf',figurename1);

    end % end figure 1


    if isfigure_angle(2)

        figure(2) % tilt angle, x-axis -> date

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

        figurename2 = sprintf('%s/%s_tilt_angle_date.pdf',figoutpath,netsta);
        print('-dpdf',figurename2);

    end % end figure 2

    if isfigure_orient(1)

        plot_tilt(idays,day_deploy,rotors,rotcohs,isfigure_orient(1));

        if isleapyear(str2double(day_deploy(1:4)))
            year = 366;
        else
            year = 365;
        end

        figure(3) % % tilt orientation, x-axis -> days since deployment
        subplot(211)
        ylim([0 360]);
        set(gca,'YTick',0:60:360);
        ylabel('Degrees from H1');
        title(sprintf('%s %s Tilt orientation',network,station));
        clim([0 year]);
        colormap jet
        colorbar
        ylabel(colorbar,'Day of year');
        box on
        grid on

        subplot(212)
        ylim([0 1]);
        set(gca,'YTick',0:0.2:1);
        ylabel('Coherence');
        xlabel(sprintf('%s (%s)','Days since deployment',deploydate));
        clim([0 year]);
        colormap jet
        colorbar
        ylabel(colorbar,'Day of year');
        box on
        grid on

        figurename3 = sprintf('%s/%s_tilt_orient_days.pdf',figoutpath,netsta);
        print('-dpdf',figurename3);

    end % end figure 3

    if isfigure_orient(2)

        figure(4) % tilt orientation, x-axis -> date

        subplot(2,1,1)
        plot(x,rotors,'o','MarkerFaceColor','#2c6db3','MarkerSize',5,'MarkerEdgeColor','none');
        xlim([deploynum recoverynum]);
        datetick('x','mmmdd','keepticks');
        ylim([0 360]);
        set(gca,'YTick',0:60:360);
        ylabel('Degrees from H1');
        title(sprintf('%s %s Tilt orientation',network,station));
        box on
        grid on

        subplot(2,1,2)
        plot(x,rotcohs,'o','MarkerFaceColor','#2c6db3','MarkerSize',5,'MarkerEdgeColor','none');
        xlim([deploynum recoverynum]);
        ylim([0 1]);
        datetick('x','mmmdd','keepticks');
        ylabel('Coherence');
        xlabel(sprintf('Date of year %s to %s',day_deploy(1:4),day_recovery(1:4)));
        set(gca,'YTick',0:0.2:1);
        box on
        grid on

        figurename4 = sprintf('%s/%s_tilt_orient_date.pdf',figoutpath,netsta);
        print('-dpdf',figurename4);

    end % end figure 4


end % end station loop
