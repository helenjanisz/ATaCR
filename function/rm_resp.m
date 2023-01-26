function outtrace = rm_resp(intrace,eventid,lo_corner,npoles,pole_zero_dir)
%% function to remove instrument response of irisFetch data structure or SACPZ files
% written by Ge Jin, 2014/02/27, adapted by H Janiszewski

intrace = intrace(1);
if isempty(intrace.sacpz.poles) % if sacpz does not exist in trace data, look for SACPZ file
    if isempty(pole_zero_dir)
        pzexist=0;
        disp('ERROR MISSING SACPZ FILE');
    else
        % find SACPZ file for station (CAREFUL, ASSUMES ONE PER STATION FOR
        % NOW)
        filesuff = sprintf('SACPZ.%s.%s.--.%s',intrace.network,intrace.station,intrace.channel);
        sacpz_filenames = dir(fullfile(pole_zero_dir,'/',filesuff));
        if length(sacpz_filenames) == 1
            pzfn_good=[pole_zero_dir '/' sacpz_filenames(1).name];
            pzexist = 1;
        else
            disp('ERROR BAD SACPZ FILE');
            pzexist = 0;
        end
    end
else
    pzexist = 0;
end

% add back in reading a SAC PZ file, but clean it up

A(1,1) = cellstr(intrace.station);
A(1,2) = cellstr(intrace.channel);
A(1,3) = cellstr(eventid);
A(1,4) = cellstr(intrace.network);

isfigure = 0;


data = intrace.data;
data = data - mean(data); % demean 
data = detrend(data);
data = flat_hanning_win(1:length(data),data,1,length(data),50);

N = intrace.sampleCount;
delta = 1/intrace.sampleRate;
T = N*delta;

if mod(N,2)
     faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/T);
else
     faxis = [0:N/2,-N/2+1:-1]*(1/T);
end

poles = intrace.sacpz.poles;
zeros = intrace.sacpz.zeros;
gain = intrace.sacpz.constant;

%%%%%%%%%%
% My edits haj
%%%%%%%%%%

if pzexist == 1 % fix this... trust the matlab one if it exists, if not look for SAC pz file
    [zz,pp,constant] = read_sac_pole_zero(pzfn_good);
    
    poleseq=isequal(poles,pp);
    zeroseq=isequal(zeros,zz);
    gaineq=isequal(gain,constant);
    
    if poleseq+zeroseq+gaineq~=3
        zeros=zz;
        poles=pp;
        gain=constant;
        
        % not sure what this B variable does
        B(1,1) = cellstr(intrace.station);
            B(1,2) = cellstr(intrace.channel);
            B(1,3) = cellstr(eventid);
            B(1,4) = cellstr(intrace.network);
            if poleseq~=1
                B(1,5) = cellstr('poles');
            elseif zeroseq~=1
                B(1,5) = cellstr('zeros');
            elseif gaineq~=1
                B(1,5) = cellstr('gain');
            end
    end
end

w = faxis.*2*pi;
resp = ones(size(w));
if isempty(zeros) == 1
    zeros = 0;
end
for ip = 1:length(poles)
	resp = resp./(i*w - poles(ip));
end
for ip = 1:length(zeros)
	resp = resp.*(i*w - zeros(ip));
end
resp = resp*gain;
if isfigure
	figure(35)
	clf
	set(gcf,'position',[360   514   900   400]);
	hold on
	subplot(1,2,1)
	set(gca,'fontsize',18)
	semilogy(faxis,abs(resp),'rx');
	subplot(1,2,2)
	set(gca,'fontsize',18)
	plot(faxis,angle(resp),'rx');
end

lo_w=2*pi*lo_corner;
hpfiltfrq=( ((w./lo_w).^(2*npoles))./(1+(w./lo_w).^(2*npoles)) );
norm_trans=hpfiltfrq./resp;    % this is normalization transfer function
norm_trans(find(isnan(norm_trans))) = 0;

fftdata = fft(data);
fftdata = fftdata(:).*norm_trans(:);
data_cor = real(ifft(fftdata));

outtrace = intrace;
outtrace.data_cor = data_cor;

disp(['Station: ',intrace.station,'.',intrace.channel,' deconv to ',intrace.sacpz.units]);

return

