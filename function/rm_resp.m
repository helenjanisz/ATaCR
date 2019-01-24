function outtrace = rm_resp(intrace,eventid,lo_corner,npoles,pole_zero_dir)
%% function to remove instrument response of irisFetch data structure or SACPZ files
% written by Ge Jin, 2014/02/27, adapted by H Janiszewski
intrace = intrace(1);

pole_zero_file=['SAC_PZs_',intrace.network,'_',intrace.station,'_',intrace.channel,'_*'];
pole_zero_file_names=dir(fullfile(pole_zero_dir,pole_zero_file));
eqjdy=num2str(dayofyear(str2num(eventid(1:4)),str2num(eventid(5:6)),str2num(eventid(7:8)),0,0));
eqyr=eventid(1:4);
eqhr=eventid(9:10);
eqmn=eventid(11:12);
if str2num(eqjdy)>=100 % possible problem in this.
eqjdt=str2num([eqyr,eqjdy,eqhr,eqmn]);
elseif str2num(eqjdy)>=10 & str2num(eqjdy)<=99
    eqjdt=str2num([eqyr,'0',eqjdy,eqhr,eqmn]);
elseif str2num(eqjdy)<=9
    eqjdt=str2num([eqyr,'00',eqjdy,eqhr,eqmn]);
end
pzexist=0;
for nfile=1:length(pole_zero_file_names)
    pole_zero_file_name=pole_zero_file_names(nfile).name;
    if strcmp(intrace.station,'ELW')==1
        datesindex(1)=strfind(pole_zero_file_name,'_1');
        datesindex(2)=strfind(pole_zero_file_name,'_2');
    else
    datesindex=strfind(pole_zero_file_name,'_2');
    end
    yrst=pole_zero_file_name(datesindex(1)+1:datesindex(1)+4);
    dyst=pole_zero_file_name(datesindex(1)+6:datesindex(1)+8);
    hrst=pole_zero_file_name(datesindex(1)+10:datesindex(1)+11);
    mnst=pole_zero_file_name(datesindex(1)+13:datesindex(1)+14);
    yred=pole_zero_file_name(datesindex(2)+1:datesindex(2)+4);
    dyed=pole_zero_file_name(datesindex(2)+6:datesindex(2)+8);
    hred=pole_zero_file_name(datesindex(2)+10:datesindex(2)+11);
    mned=pole_zero_file_name(datesindex(2)+13:datesindex(2)+14);
    julst=str2num([yrst,dyst,hrst,mnst]);
    juled=str2num([yred,dyed,hred,mned]);
    if eqjdt>=julst && eqjdt<=juled
        %may need to check sample rate - check what chan for COR reads
        %as... might not be a problem for now, COR may not get included
        pzfn_good=[pole_zero_dir '/' pole_zero_file_name];
        pzexist=1;
        break
    else
        if nfile == length(pole_zero_file_names)
            pzexist=0;
            A(1,1) = cellstr(intrace.station);
            A(1,2) = cellstr(intrace.channel);
            A(1,3) = cellstr(eventid);
            A(1,4) = cellstr(intrace.network);
%             if ~exist(sprintf('%s%s',outdir,nopz_file),'file') ;
%                 a=A;
%                 save(sprintf('%s%s',outdir,nopz_file),'a');
%                 clear A
%             elseif exist(sprintf('%s%s',outdir,nopz_file),'file') ;
%                load(sprintf('%s%s',outdir,nopz_file))
%                a = vertcat(a,A);
%                save(sprintf('%s%s',outdir,nopz_file),'a');
%                clear A
%             end
%             disp(['Missing SAC PZ file, ',' ',intrace.network,' ', intrace.station,' ',intrace.channel,'. Using downloaded, proceed with caution.' ])
        end
    end
end

isfigure = 0;

data = intrace.data;
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

if pzexist == 1
    [zz,pp,constant] = read_sac_pole_zero(pzfn_good);
    
    poleseq=isequal(poles,pp);
    zeroseq=isequal(zeros,zz);
    gaineq=isequal(gain,constant);
    
    if poleseq+zeroseq+gaineq~=3
        zeros=zz;
        poles=pp;
        gain=constant;
        
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
%             if ~exist(sprintf('%s%s',outdir,erpz_file),'file') ;
%                 b=B;
%                 save(sprintf('%s%s',outdir,erpz_file),'b');
%                 clear B
%             elseif exist(sprintf('%s%s',outdir,erpz_file),'file') ;
%                load(sprintf('%s%s',outdir,erpz_file))
%                b = vertcat(b,B);
%                save(sprintf('%s%s',outdir,erpz_file),'b');
%                clear B
%             end
            
%         disp(['PZ file mismatch ', intrace.station,' ',intrace.channel,'.' ])
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

