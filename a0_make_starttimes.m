% make_starttimes
%
% Make start time files for daily data downloads from list of event times in 
% format YYYYmmddhhMM
%
% J. Russell & H. Janiszewski 
% hajanisz@hawaii.edu
% updated 11/19

% Input eventfile
evfile = 'NOISETC_CI/eventtimes_CItest.txt';

% Output dayfile
dayfile = 'NOISETC_CI/starttimes_CItest.txt';

% Number of days prior to event for calculation of spectra
Ndays = 4;

% Load event list
evlist = textread(evfile,'%s');

fid = fopen(dayfile,'w');
for iev = 1:length(evlist)
    evnum = datenum(evlist(iev),'yyyymmddHHMM');
    daynums = flip(evnum - [1:Ndays]);
    for iday = 1:length(daynums)
        day = datestr(daynums(iday),'yyyymmddHHMM');
        fprintf(fid,'%s\n',day);
    end    
end
fclose(fid);

