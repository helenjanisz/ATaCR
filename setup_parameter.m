% parameters for the NoiseTiltComp_pck

% path for matlab codes and functions
addpath ('function');
% location of the continuous matlab files for spectral properties
WORKINGdir = 'NOISETC_SAMPLE_CI/DATA/datacache_day_preproc/';
% output directory for spectra
OUTdir = 'NOISETC_SAMPLE_CI/DATA/NOISETC';
 % directory for figure output
FIGdir = 'NOISETC_SAMPLE_CI/FIGURES/NOISETC';


% information for station
network = '7D'; 
station = 'M08A'; 

% Channel Names
chz_vec = {['HHZ'], ['BHZ'], ['LHZ']}; % list of acceptable names for Z component
ch1_vec = {['HH1'], ['BH1'], ['LHN'], ['LH1']}; % list of acceptable names for H1 component
ch2_vec = {['HH2'], ['BH2'], ['LHE'], ['LH2']}; % list of acceptable names for H2 component
chp_vec = {['HDH'], ['BDH']}; % list of acceptable names for P component

% Spectral Properties Windowing
T    = 7200;  % the legnth of each time window, in sec, should match the event data length for corrections
overlap = 0.3; %fraction of window overlap for spectra calculation

% Quality Control Parameters for Daily Windows
pb = [0.004 .2]; %pass-band, in Hz
tolerance = 1.5;
a_val = 0.05;
minwin = 10;    % minimum numbers of time window for segment to be accepted

% Tilt orientation - only matters if using transfer functions with the 'H'
% option, but package needs variables specified to run; leave as default if
% not using
tiltfreq = [.005, .035]; % specifying frequency ranges for maximum coherence search
% tiltfreq = [1/20, 1/5]; % specifying frequency ranges for maximum coherence search

% Quality Control Parameters for Deployment Days
pb_dep = [0.004 .2]; %pass-band, in Hz
tolerance_dep = 1.5;
a_val_dep = 0.05;

% Transfer Function Options
TF_list = {'ZP','Z1','Z2-1','ZP-21','ZH','ZP-H'};
% convention: Z1 = Z/1, Z - signal, 1 - noise
%             Z2-1 = Z - signal, 2 - noise, with Z1 already taken into account
%             ZH - Z - signal, H - rotated max horizontal noise direction,
% Typical options for tilt compliance removal: 'ZP','1P','2P','Z1','Z2-1','ZP-21','ZH','ZP-H'
% Note: should be able to put in other options (e.g. Z2-1P, removing
% compliance then component wise tilt, but those not thoroughly tested yet)
% Currently rotation only implemeted for ZH

% Correction Options
taptime = 0.075; % taper for seismogram correction step

tf_op = 1; %option for using either daily (1) or station average (2) TF in correction

filop = 2; %how to filter TF
% 1 - user specified constant high pass and low pass
% 2 - %lowpass - 0.005+freqcomp, adaptive to the infragravity wave cutoff, no high pass;
tap_width = 0.01; %width in percent of frequency vector of cosine taper is this actually used?????
taper_lim = [0 1]; % upper and lower freuqncy in Hz limit for taper if user specified (option 1); zero means not applied
