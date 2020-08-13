% data_preprocess
% Collection of sample data pre-processing codes that may be applied to
% data prior to noise analysis and tilt/compliance correction. This script
% will apply the pre-processing steps and save the data in an appropriate
% file structure.
% Pre-procssing steps included here are:
%  - Downsampling
%  - Response Removal
%  - Gain Adjustment
% Each may be included or excluded depending on the needs of the user.
% Make sure preprocessing steps for a given station 
% HAJ July 2016

addpath ('function');
INPUTdir = 'NOISETC_CI/DATA/datacache/';
OUTPUTdir = 'NOISETC_CI/DATA/datacache_prepro/';

% pole_zero_dir='directory_where_SACPZ_files_are_here'; % for response removal if using sac PZ files
pole_zero_dir=''; % if not using leave blank

network = '7D';
stations = textread('./NOISETC_CI/stalist.txt','%s');
% channels = {'HHZ','HH1','HH2','HDH'};
% channels = {'HXZ','HX1','HX2','HXH'};
channels = {'BHZ','BH1','BH2','BDH'};
% channels = {'HHZ','HH1','HH2','BDH'};

% Example for case where no gain correction is applied and response is
% removed from all channels
%%%%
resprm = [1 1 1 1]; % for each channel 1 means remove response, 0 means don't. Order matches channels vector
gaincorr = [1 1 1 1]; % gain correction factor to multiply data by for each channel
samprate = 5; % new sample rate for each channel (note: for tilt/comp package must all be equal)
hp_filt = [0 0 0 0]; % apply high pass filter

% Example for case where gain correction is applied and response is only
% removed from seismometer channels
%%%%
% resprm = [1 1 1 0]; % for each channel 1 means remove response, 0 means don't. Order matches channels vector
% gaincorr = [2.37 2.37 2.37 2.37]; % gain correction factor to multiply data by for each channel
% samprate = 5; % new sample rate for each channel (note: for tilt/comp package must all be equal)
% hp_filt = [0 0 0 1]; % apply high pass filter

% parameters for high pass for response removal
lo_corner = 0.001;  % in Hz
        npoles=5;

%%%%% end user input parameters %%%%%

if ~exist(OUTPUTdir)
    mkdir(OUTPUTdir);
end

% data_filenames = dir(fullfile(INPUTdir,network,station,'/',['*.mat']));
data_filenames = dir(fullfile(INPUTdir));
for ista = 1:length(stations)
    station = stations{ista};
    for ie = 1 : length(data_filenames) % begin file loop for each station
        inpath = sprintf('%s%s/',INPUTdir,data_filenames(ie).name);
        if ~isdir(inpath)
            continue
        end
        if length(data_filenames(ie).name)~=12;
            continue
        end
        station_filenames = dir(fullfile(inpath,['*.mat']));
        if isempty(station_filenames)
            continue
        end
        for is = 1:length(station_filenames)
            if isempty(findstr([network,'_',station],station_filenames(is).name))
                continue
            end
            load(fullfile(INPUTdir,'/',data_filenames(ie).name,'/',station_filenames(is).name));
        traces_new = traces;
        eventid = data_filenames(ie).name(1:12);
        if ~exist([OUTPUTdir,eventid])
        mkdir([OUTPUTdir,eventid]);
    end
        prob=0;
        for ic = 1:length(channels) % begin channel loop
            chan = channels(ic);
            idxch = find(ismember({traces.channel},chan));
            if length(idxch)>1
                disp('Skipping. Too many records for single channel.')
                prob = 1;
                continue
            end
            if isempty(idxch)
                continue
            end
            chan_data = traces(idxch);
            rate = chan_data.sampleRate;
            dt = 1/rate;
            %%%%%%%%%%%%%%%%%
            % REMOVE RESPONSE
            %%%%%%%%%%%%%%%%%
            if resprm(ic) ==1
                chan_data = rm_resp(chan_data,eventid,lo_corner,npoles,pole_zero_dir);
                data_raw = chan_data.data_cor;
            else
                data_raw = chan_data.data;
            end
            %%%%%%%%%%%%%%%%%
            % GAIN CORRECTION
            %%%%%%%%%%%%%%%%%
            data_raw = gaincorr(ic).*data_raw;
             %%%%%%%%%%%%%%%%%
            % HIGH PASS FILTER - should link this and the bit in rm_resp
            % function to a different function, so they are the same, but the
            % user can specify parameters easily
            %%%%%%%%%%%%%%%%%
            if hp_filt(ic) ==1
            
            lo_w=2*pi*lo_corner;
            
            N = length(data_raw);
            delta = dt;
            Tr = N*delta;
            
            if mod(N,2)
                faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/Tr);
            else
                faxis = [0:N/2,-N/2+1:-1]*(1/Tr);
            end
            w = faxis.*2*pi;
            
            hpfiltfrq=( ((w./lo_w).^(2*npoles))./(1+(w./lo_w).^(2*npoles)) );
            norm_trans=hpfiltfrq;    % this is normalization transfer function
            norm_trans(find(isnan(norm_trans))) = 0;
            
            fftdata = fft(data_raw);
            fftdata = fftdata(:).*norm_trans(:);
            data_raw = real(ifft(fftdata));
            end
            %%%%%%%%%%%%%%%%%
            % DOWNSAMPLING
            %%%%%%%%%%%%%%%%%
            if rate == samprate
                taxis = [0:dt:(length(data_raw)-1)*dt]';
            else
                dt_new = 1/samprate;
                [data_raw,taxis] = resample(data_raw,dt_new,dt);
            end
            
            traces_new(idxch).data = data_raw;
            traces_new(idxch).sampleRate = samprate;
            traces_new(idxch).sampleCount = length(data_raw);
        end % end channel loop
        %save good files
        if prob ==0
            traces = traces_new;
            filename = fullfile(OUTPUTdir,eventid,'/',station_filenames(is).name);
            save(filename,'traces');
        end
        end   
    end % end file loop
end % end station loop
