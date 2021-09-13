% eventmat2sac
%
% Convert *.mat files containing corrected seismograms to SAC files matching
% the original input SAC files.
%
% J. Russell & H. Janiszewski 
% hjaniszewski@carnegiescience.edu
% updated 11/19
%
% W. Hawley modified to use original sac start time information instead of
% matlab generated filename... there was in issue using the wrong start
% time for some files
% updated 09/20
%
% Will try to read header information from original uncorrected sac files at inpath_uncorr,
% otherwise will read header information from preprocessed mat files at inpath_eventmat.
% Updated 09/21 - J. Russell
%
 
clear;

setup_parameter;

inpath_uncorr = 'path/to/local/event/sac/files/';
inpath_eventmat = 'NOISETC_CI/DATA/datacache_prepro/'; % path to preprocessed event mat files
str_corr = 'ZP-21'; % String for correction to output
channel = 'BHZ';

%% Load data

if tf_op == 1
    corrseis_path = sprintf('%s/CORRSEIS/',OUTdir);
elseif tf_op ==2
    corrseis_path = sprintf('%s/CORRSEISAVTF/',OUTdir);
end

stadirs = dir(fullfile(corrseis_path));
for ista = 1:length(stadirs)
    station = stadirs(ista).name;
    inpath_corr = sprintf('%s/%s/',corrseis_path,station);
    if ~isdir(inpath_corr)
        continue
    end
    filenames_corr = dir(fullfile(inpath_corr,['*.mat']));
    disp(station);
    % Loop over event files
    for iev = 1:length(filenames_corr);
        if ~exist(fullfile(inpath_corr,filenames_corr(iev).name))
            continue
        end
        load(fullfile(inpath_corr,filenames_corr(iev).name))
        corr_idx = find(strcmp({corrseis.label},str_corr));
        corrdata = corrseis(corr_idx).timeseries;
        
        % Load data headers
        sacfile = fullfile(sprintf('%s/%s/%s.%s.%s.%s.sac',inpath_uncorr,corrected.params.eventid, corrected.params.eventid, corrected.params.network, corrected.params.station, channel));
        if exist(sacfile) % Check for original sac file to get header information
            sacin = rdsac(sacfile);
        else
            % If original SAC file does not exist, read required header info from the mat file
            matfile = fullfile(sprintf('%s/%s/%s_%s_%s.mat',inpath_eventmat,corrected.params.eventid, corrected.params.eventid, corrected.params.network, corrected.params.station));
            matin = load(matfile);
            icomp = find(strcmp({matin.traces.channel},channel));
            sacin.HEADER.KNETWK = matin.traces(icomp).network;
            sacin.HEADER.KSTNM = matin.traces(icomp).station;
            sacin.HEADER.KCMPNM = matin.traces(icomp).channel;
            sacin.HEADER.KHOLE = matin.traces(icomp).location;
            sacin.HEADER.STLA = matin.traces(icomp).latitude;
            sacin.HEADER.STLO = matin.traces(icomp).longitude;
            sacin.HEADER.STEL = matin.traces(icomp).elevation;
            sacin.HEADER.T0 = matin.traces(icomp).startTime;
            
            % header origin time values
            tv = datevec(sacin.HEADER.T0(1));
            sacin.HEADER.NZYEAR = tv(1);
            sacin.HEADER.NZJDAY = datenum(tv(1:3)) - datenum(tv(1),1,1) + 1;
            sacin.HEADER.NZHOUR = tv(4);
            sacin.HEADER.NZMIN = tv(5);
            sacin.HEADER.NZSEC = floor(tv(6));
            sacin.HEADER.NZMSEC = (tv(6) - sacin.HEADER.NZSEC)*1e3;
        end
        disp(corrected.params.eventid);
        for itr = 1 %1:length(traces)
                if tf_op ==1
                    opath = sprintf('%s/CORRSEIS_SAC/%s/',OUTdir,corrected.params.eventid);
                elseif tf_op ==2
                    opath = sprintf('%s/CORRSEISAVTF_SAC/%s/',OUTdir,corrected.params.eventid);
                end
                if ~exist(opath)
                    mkdir(opath);
                end
                H = sacin.HEADER;
                H.DELTA = corrected.params.dt;
                H.NPTS = length(corrected.params.taxis);
                
                data = corrseis(corr_idx).timeseries;
%                 evid =  datestr(traces(itr).startTime,'yyyymmddhhMM');
                evid = corrected.params.eventid;
                startDate = datetime(H.NZYEAR,1,H.NZJDAY); % use this to convert jday to month - day
                startYear = num2str(year(startDate), '%04d');
                startMonth = num2str(month(startDate), '%02d');
                startDay = num2str(day(startDate), '%02d');
                fullevid = [startYear,startMonth,startDay,num2str(H.NZHOUR,'%02d'),num2str(H.NZMIN,'%02d'),num2str(H.NZSEC,'%02d'),num2str(H.NZMSEC,'%03d')];
                startTime = datenum(fullevid,'yyyymmddhhMMSSFFF');
                sac_path = fullfile(sprintf('%s/%s.%s.%s.%s.sac',opath,evid, sacin.HEADER.KNETWK, sacin.HEADER.KSTNM, sacin.HEADER.KCMPNM));
                mksac(sac_path,data,startTime,H);                
        end       
        
    end
    
end
