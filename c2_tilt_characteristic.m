% tilt_characteristic
%
% Calculate and plot the tilt characteristics of seismometers during the 
% date input by the user.
% The tilt characteristic depends on the selection of the tilt frequency 
% band, which can be changed in the setup_paremeter. In the tilt frequency 
% band, the coherence is the largest, the admittance is nearly constant 
% and the phase shift nearly zero.
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

clear; close all;
startDate = '2012-03-01';
endDate   = '2012-03-30';

INPUTdir = 'NOISETC_CI/DATA/NOISETC/SPECTRA';


% DO NOT EDIT BELOW
setup_parameter;
startnum = datenum(startDate,'yyyy-mm-dd');
endnum   = datenum(endDate,'yyyy-mm-dd');
figoutpath = sprintf('%s/TILT',FIGdir);

if ~exist(figoutpath,'dir')
    mkdir(figoutpath);
end

for ista = 1:length(stations) % begin station loop
    station = stations{ista};
    netsta = [network,station];
    close all;
    inpath = [INPUTdir,'/',netsta];

    id = 0;
    for iday = startnum:endnum % begin date loop
        id = id+1;

        dayid = sprintf('%s0000',datestr(iday,'yyyymmdd'));
        spectra_file = sprintf('%s/%s_%s_spectra.mat',inpath,netsta,dayid);

        if exist(spectra_file,'file')
            fprintf('loading %s\n',spectra_file);
            load(spectra_file);

            f = specprop.params.f;
            npts_smooth = floor(specprop.params.NFFT/1000)+1;

            % Coherence between Z and H
            coh_stack(:,id) = smooth(abs(specprop.rotation.chz_stack).^2./(specprop.rotation.chh_stack.*specprop.power.czz_stack),npts_smooth);

            % Admittance between Z and H
            ad_stack(:,id) = smooth(abs(specprop.rotation.chz_stack)./specprop.rotation.chh_stack,npts_smooth);

            tiltf = find(f>=tiltfreq(1) & f<=tiltfreq(2));

            % Tilt angles are the arctangent of the admittance
            tiltang = atand(mean(ad_stack(tiltf,id)));
            tiltcoh = mean(coh_stack(tiltf,id));

        else
            fprintf('%s does not exist!\nTilt angle and coherence are assigned to NaN\n',spectra_file);
            tiltang = NaN;
            tiltcoh = NaN;
        end

        days(id,:) = iday;
        tiltangs(id,:) = tiltang;
        tiltcohs(id,:) = tiltcoh;
        
    end % end date loop


    figure(1)

    subplot(2,1,1)
    plot(days,tiltangs,'o','MarkerFaceColor','#00688B','MarkerSize',5,'MarkerEdgeColor','none');
    xlim([startnum endnum]);
    datetick('x','mmmdd','keepticks');
    ylabel('Tilt angle (°)','FontSize',12);
    grid on;

    subplot(2,1,2)
    plot(days,tiltcohs,'o','MarkerFaceColor','#00688B','MarkerSize',5,'MarkerEdgeColor','none');
    xlim([startnum endnum]);
    ylim([0 1]);
    datetick('x','mmmdd','keepticks');
    ylabel('Coherence','FontSize',12);
    xlabel(sprintf('Date of year %s to %s',startDate(1:4),endDate(1:4)),'FontSize',12);
    set(gca,'YTick',0:0.2:1,'FontSize',10);
    grid on;

    % fpos=get(gcf,'Position'); fpos(3)=fpos(3)*1.3; set(gcf,'Position',fpos);
    filename=sprintf('%s/%s_tilt.pdf',figoutpath,netsta);
    print('-dpdf',filename)

end % end station loop
