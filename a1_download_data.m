% download_data

% downloads the data files used in calculating the noise spectra and the
% transfer functions for tilt and compliance corrections and saves them as
% matfiles (default is 24 hours of data in each file).

% H. Janiszewski 
% hajanisz@hawaii.edu
% updated 2/18

clear;

javaaddpath('IRIS-WS-2.0.19.jar');

startlist = 'starttimes.txt'; % list of start times for data download
datalength = 86400; % length of time series after each start time in seconds (default 86400, code not thoroughly tested for other values)

download_networks = '7D'; % list of networks to download
download_stations = textread('stalist.txt','%s'); % list of stations to download (* for all)

% Channel Names
chz_vec = 'HHZ,BHZ'; % list of acceptable names for Z component
ch1_vec = 'HH1,BH1'; % list of acceptable names for H1 component
ch2_vec = 'HH2,BH2'; % list of acceptable names for H2 component
chp_vec = 'HDH,BDH'; % list of acceptable names for P component

datacache = 'DATA/datacache_day'; % output folder for data

%%%%% end user input parameters %%%%%

if ~exist(datacache,'dir')
    mkdir(datacache)
end

startlist = textread(startlist,'%s');
chanlist = sprintf('%s,%s,%s,%s',chz_vec,ch1_vec,ch2_vec,chp_vec);

for id = 1:length(startlist)
   eventid = cell2mat(startlist(id));
   disp(sprintf('Start Time: %s',eventid));
   otime = datenum(eventid,'yyyymmddHHMM');
   starttime = datestr(otime,'yyyy-mm-dd HH:MM:SS');
   endtime = datestr(otime+datalength/3600/24,'yyyy-mm-dd HH:MM:SS');
   
   stations_info = irisFetch.Stations('channel',download_networks,download_stations,'*',chz_vec,'startTime',starttime,'endTime',endtime);
   
   
   for ista =1:length(stations_info)
       error = 0;
       stnm = stations_info(ista).StationCode;
       network = stations_info(ista).NetworkCode;
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
       disp(['Downloading station: ',stnm,' From:',starttime,' To:',endtime]);
		try
            traces_day = irisFetch.Traces(network,stnm,'*',chanlist,starttime,endtime,'includePZ');
			save(sta_filename,'traces_day');
		catch e
            e.message;
            error = 1;
        end
        if error ==1
            try % to try and get around the missing zeros for some pressure components                
                traces_day = irisFetch.Traces(network,stnm,'*',chanlist,starttime,endtime);
                save(sta_filename,'traces_day');
            catch e
                e.message;
                continue;
            end
        end
    end
   
end