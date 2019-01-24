function [Zraw,H1raw,H2raw,Praw,taxisZ,taxis1,taxis2,taxisP,dt] = varsetup_correctevent(sta,chz_vec,ch1_vec,ch2_vec,chp_vec)

idxZ = find(ismember({sta.traces.channel},chz_vec));
idx1 = find(ismember({sta.traces.channel},ch1_vec));
idx2 = find(ismember({sta.traces.channel},ch2_vec));
idxP = find(ismember({sta.traces.channel},chp_vec));

Z = sta.traces(idxZ);
H1 = sta.traces(idx1);
H2 = sta.traces(idx2);
P = sta.traces(idxP);
rateZ = Z.sampleRate;
rate1 = H1.sampleRate;
rate2 = H2.sampleRate;
rateP = P.sampleRate;
dtZ = 1/rateZ;
dt1 = 1/rate1;
dt2 = 1/rate2;
dtP = 1/rateP;

Zraw = Z.data;
H1raw = H1.data;
H2raw = H2.data;
Praw = P.data;

% check all samp rates are equal, if not resample at highest common
% denominator
samprate = [rateZ,rate1,rate2,rateP];
if length(unique(samprate))>1
    [rates,is,ir] = unique(samprate);
    newrate = rates(1);
    for ii=2:length(rates)
        newrate = gcd(newrate,rates(ii));
    end
    dt = 1/newrate;
    [Zraw,taxisZ] = resample(Zraw,dt,dtZ);
    [H1raw,taxis1] = resample(H1raw,dt,dt1);
    [H2raw,taxis2] = resample(H2raw,dt,dt2);
    [Praw,taxisP] = resample(Praw,dt,dtP);
else dt = 1/rateZ;
    taxisZ = [0:dt:(length(Zraw)-1)*dt]';
    taxis1 = [0:dt:(length(H1raw)-1)*dt]';
    taxis2 = [0:dt:(length(H2raw)-1)*dt]';
    taxisP = [0:dt:(length(Praw)-1)*dt]';
end

data_length = [length(Zraw),length(H1raw),length(H2raw),length(Praw)];
if length(unique(data_length))>1
    disp('Problem. Different data lengths.')
    Zraw=nan; H1raw=nan; H2raw=nan; Praw=nan; taxisZ=nan; taxis1=nan; taxis2=nan; taxisP=nan; dt=nan;
    return
end

if data_length*dt<7200
    disp('Not enough data. Skipping.')
    Zraw=nan; H1raw=nan; H2raw=nan; Praw=nan; taxisZ=nan; taxis1=nan; taxis2=nan; taxisP=nan; dt=nan;
    return
end

return