% sac2mat_data
%
% Following download_data.m script but uses local data files. Loads in 24
% hour sac data files and puts them in proper .mat structures. This assumes
% that instrument response has already been removed and is in the same
% units as the event files.
%
% Assumed naming convention for local data:
% path/to/day/data/{station}/{station}.{yyyy}.{jday}.{hh}.{MM}.{ss}.{component}.sac
%
% J. Russell & H. Janiszewski 
% hajanisz@hawaii.edu
% updated 11/19

clear;

addpath ('function');

startlist = 'NOISETC_CI/starttimes_CItest.txt'; % list of start times for data download
datalength = 86400; % length of time series after each start time in seconds (default 86400, code not thoroughly tested for other values)
sacdaydata = '/path/to/local/day/sac/files/'; % path to local day sac files

download_networks = '7D'; %'2D'; % list of networks to download
download_stations = textread('./NOISETC_CI/stalist.txt','%s'); % list of stations to download (* for all)

% Channel Names
chz_vec = 'BHZ'; % list of acceptable names for Z component
ch1_vec = 'BH1'; % list of acceptable names for H1 component
ch2_vec = 'BH2'; % list of acceptable names for H2 component
chp_vec = 'BDH'; % list of acceptable names for P component

datacache = 'NOISETC_CI/DATA/datacache_day'; % output folder for data

%%%%% end user input parameters %%%%%

if ~exist(datacache,'dir')
    mkdir(datacache)
end

startlist = textread(startlist,'%s');

for id = 1:length(startlist)
   eventid = cell2mat(startlist(id));
   disp(sprintf('Start Time: %s',eventid));
   otime = datenum(eventid,'yyyymmddHHMM');
   starttime = datestr(otime,'yyyy-mm-dd HH:MM:SS');
   endtime = datestr(otime+datalength/3600/24,'yyyy-mm-dd HH:MM:SS');   
   jday = floor(otime - datenum(year(otime),1,1) + 1);
   
   for ista =1:length(download_stations)
       clear traces_day
       error = 0;
       stnm = download_stations{ista};
       network = download_networks;
       if ~exist(fullfile(datacache,network),'dir')
           mkdir(fullfile(datacache,network));
       end
       if ~exist(fullfile(datacache,network,stnm),'dir')
           mkdir(fullfile(datacache,network,stnm));
       end
       sta_filename = fullfile(datacache,network,stnm,[eventid,'_',network,'_',stnm,'.mat']);
       if exist(sta_filename,'file')
           disp(['Exist: ',sta_filename,', Skip!']);
           continue;
       end
       disp(['SAC to MAT station: ',stnm,' From:',starttime,' To:',endtime]);
		try
            ich = 0;
            for ch = {chp_vec ch1_vec ch2_vec chz_vec}
                ich = ich + 1;
                sac_filename = [stnm,'.',num2str(year(otime)),'.',num2str(jday,'%03d'),'.00.00.00.',ch{:},'.sac'];
                sac = rdsac(fullfile(sacdaydata,stnm,sac_filename));
                traces_day(ich) = sac2mat( sac );
            end
            save(sta_filename,'traces_day');
		catch e
            e.message;
            display('Missing data file');
            error = 1;
        end
    end
   
end