% sac2mat_event
%
% Following download_event.m script but uses local data files. Loads in 
% sac data files and puts them in proper .mat structures. This assumes
% that instrument response has already been removed and is in the same
% units as the day files.
% Waveform lengths should match T in setup_parameters
%
% Assumed naming convention for local data:
% path/to/event/data/{yyyymmddhhMM}/{yyyymmddhhMM}.{network}.{station}.{component}.sac
%
% J. Russell & H. Janiszewski 
% hjaniszewski@carnegiescience.edu
% updated 2/18

clear;

addpath ('function');
startlist = 'NOISETC_CI/eventtimes_CItest.txt'; % list of start times for data download, this will be the beginning of the waveform so specify appropriately, signal of interest should not be at edges of time series
datalength = 7200; % length of time series after each start time in seconds (default 86400, code not thoroughly tested for other values)
saceventdata = 'path/to/local/event/sac/files/'; % path to local event sac files

download_networks = '7D'; % list of networks to download
download_stations = textread('./NOISETC_CI/stalist.txt','%s');  % list of stations to download (* for all)
% Channel Names
chz_vec = 'BHZ'; %'BHZ'; % list of acceptable names for Z component
ch1_vec = 'BH1'; %'BH1'; % list of acceptable names for H1 component
ch2_vec = 'BH2'; %'BH2'; % list of acceptable names for H2 component
chp_vec = 'BDH'; %'BDH'; % list of acceptable names for P component

datacache = 'NOISETC_CI/DATA/datacache_day'; % output folder for data

%%%%% end user input parameters %%%%%

if ~exist(datacache,'dir')
    mkdir(datacache)
end

startlist = textread(startlist,'%s');

for id = 1:length(startlist)
   eventid = cell2mat(startlist(id));
   disp(sprintf('Start Time: %s',eventid));
   if length(eventid) == 12
        otime = datenum(eventid,'yyyymmddHHMM'); % look at this when I get back
   elseif length(eventid) == 14
       otime = datenum(eventid,'yyyymmddHHMMSS');
       eventid = eventid(1:12); % for naming purposes only, start time will still be saved to the second in traces file
   end
   starttime = datestr(otime,'yyyy-mm-dd HH:MM:SS');
   endtime = datestr(otime+datalength/3600/24,'yyyy-mm-dd HH:MM:SS');   
   
   for ista =1:length(download_stations)
       clear traces
       error = 0;
       stnm = download_stations{ista};
       network = download_networks;
       if ~exist(fullfile(datacache,eventid),'dir')
           mkdir(fullfile(datacache,eventid));
       end
       sta_filename = fullfile(datacache,eventid,[eventid,'_',network,'_',stnm,'.mat']);
       if exist(sta_filename,'file')
           disp(['Exist: ',sta_filename,', Skip!']);
           continue;
       end
       disp(['SAC to MAT station: ',stnm,' From:',starttime,' To:',endtime]);
		try
            ich = 0;
            for ch = {chp_vec ch1_vec ch2_vec chz_vec}
                ich = ich + 1;
                sac_filename = [eventid,'.',network,'.',stnm,'.',ch{:},'.sac'];
                sac = rdsac(fullfile(saceventdata,eventid,sac_filename));
                traces(ich) = sac2mat( sac );
            end
            save(sta_filename,'traces');
		catch e
            e.message;
            error = 1;
            display('Could not load data file');
        end
    end
   
end