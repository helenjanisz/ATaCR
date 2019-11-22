function [ mat ] = sac2mat( sac )
%SAC2MAT convert sac file from rdsac to matlab structure mimicking
% irisFetch.Traces
%
% JBR 9/24/19

mat.network = sac.HEADER.KNETWK;
mat.station = sac.HEADER.KSTNM;
mat.location = '';
mat.channel= sac.HEADER.KCMPNM;
mat.quality = '';
mat.latitude = sac.HEADER.STLA;
mat.longitude = sac.HEADER.STLO;
mat.elevation = sac.HEADER.STEL;
try
    mat.depth = sac.HEADER.EVDP;
catch
    mat.depth = 0;
end
try
    mat.azimuth = sac.HEADER.AZ;
catch
    mat.azimuth = 0;
end
mat.dip = 0;
mat.sensitivity =  0;
mat.sensitivityFrequency = 0;
mat.instrument = '';
mat.sensitivityUnits = '';
mat.data = sac.d;
mat.sampleCount = sac.HEADER.NPTS;
mat.sampleRate = sac.HEADER.DELTA;
mat.startTime = datenum(sac.HEADER.NZYEAR,0,sac.HEADER.NZJDAY);
mat.endTime = mat.startTime + sac.HEADER.DELTA*sac.HEADER.NPTS/60/60/24;
mat.sacpz = [];

end

