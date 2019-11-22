% eventmat2sac
%
% Convert *.mat files containing corrected seismograms to SAC files matching
% the original input SAC files.
%
% J. Russell & H. Janiszewski 
% hjaniszewski@carnegiescience.edu
% updated 11/19

clear;

setup_parameter;

inpath_uncorr = 'path/to/local/event/sac/files/';
corr_idx = 3;
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
        corrdata = corrseis(corr_idx).timeseries;
        
        % Load data headers
        sacin = rdsac(fullfile(sprintf('%s/%s/%s.%s.%s.%s.sac',inpath_uncorr,corrected.params.eventid, corrected.params.eventid, corrected.params.network, corrected.params.station, channel)));
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
                data = corrseis(corr_idx).timeseries;
%                 evid =  datestr(traces(itr).startTime,'yyyymmddhhMM');
                evid = corrected.params.eventid;
                fullevid = [evid,num2str(H.NZSEC,'%02d'),num2str(H.NZMSEC,'%03d')];
                startTime = datenum(fullevid,'yyyymmddhhMMSSFFF');
                sac_path = fullfile(sprintf('%s/%s.%s.%s.%s.sac',opath,evid, sacin.HEADER.KNETWK, sacin.HEADER.KSTNM, sacin.HEADER.KCMPNM));
                mksac(sac_path,data,startTime,H);                
        end       
        
    end
    
end
    