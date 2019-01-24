function [Craw,elev,rateC,dtC,startC,endC] = varsetup_dailystaspectra(traces_day_C,comp_exist);

if comp_exist ==1
C = traces_day_C;
Craw = C.data;
elev = -C.elevation;
rateC = C.sampleRate;
dtC = 1/rateC;
%filter goes here if needed
    %%%%%%%%%%% Processing to data to filter, without removing resp
%     if APG ==1
%     lo_corner = 0.005;  % in Hz
%     npoles=5;
%     lo_w=2*pi*lo_corner;
%     
%     N = length(Praw);
%     delta = dtP;
%     Tr = N*delta;
%     
%     if mod(N,2)
%         faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/Tr);
%     else
%         faxis = [0:N/2,-N/2+1:-1]*(1/Tr);
%     end
%     w = faxis.*2*pi;
%     
%     hpfiltfrq=( ((w./lo_w).^(2*npoles))./(1+(w./lo_w).^(2*npoles)) );
%     norm_trans=hpfiltfrq;    % this is normalization transfer function
%     norm_trans(find(isnan(norm_trans))) = 0;
%     
%     fftdata = fft(Praw);
%     fftdata = fftdata(:).*norm_trans(:);
%     Praw = real(ifft(fftdata));
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startC = C.startTime;
endC = C.endTime;

else
    Craw= [];
    elev=[];
    freqcomp=[];
    rateC=[];
    dtC=[];
    startC=[];
    endC=[];
    
end

return