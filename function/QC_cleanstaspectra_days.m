function [ goodwins ] = QC_staspectra_days(spect,comp_exist,f,pb,tolerance,a_val)
% [ goodwins ] = QC_staspectra_windows( Zraw,H1raw,H2raw,Praw,tt,Nwin,ifplot )
%  Function to go through the different noise windows and do a QC on which
%  ones to keep based on if any of the windows particularly increase the
%  standard deviation of the average spectrum.
%  
% HAJ from Z. Eilon 03/2016
% edit: HAJ 1/2017 for component number flexibility
% spect must be vector with format [spect_Z,spect_1,spect_2,spect_P]. Even if
% component doesn't exist should just be an empty vector

comporder = {'Z', 'H1', 'H2', 'P'};

for ic = 1:4
    if comp_exist(ic)==0
        continue
    end

[b,Nwin] = size(spect(:,:,ic));

% fwt
f_wt = ones(size(spect(:,:,ic),1),1);
f_wt(f<pb(1)) = 0;
f_wt(f>pb(2)) = 0;

a_D = abs(spect(:,:,ic));
la_D = log10(a_D).*(f_wt*ones(1,Nwin));
sla_D = zeros(size(a_D));
npts_smooth = floor(length(f)*2/1000)+1;
for ii = 1:Nwin
    sla_D(:,ii) = smooth(la_D(:,ii),npts_smooth);
    % sla_D(:,ii) = la_D(:,ii);
end
dsla_D(:,:,ic) = sla_D-mean(sla_D,2)*ones(1,Nwin);

% debugging only 
%     figure(97)
%     subplot(4,1,ic)
%     semilogx(FD,sla_D)
%     hold on
% %     
%     figure(98)
%     subplot(4,1,ic)
%     semilogx(FD,dsla_D(:,:,ic))
%     hold on
%
end
%% cycle through trying to kill high-std-norm windows
moveon = false;
goodwins = [1:Nwin]';
while moveon==false
    clear('normvar', 'ubernorm');
    in = 1;
    for ic = 1:4
        if comp_exist(ic)==1
            for ii = 1:length(goodwins)
                ind = goodwins; ind(ii) = [];
                normvar(ii) = norm(std(dsla_D(:,ind,ic),0,2));
            end
            ubernorm(:,in) = [normvar];
            in = in+1;
        end
        
    end
    penalty = ones(length(goodwins),1)*median(ubernorm) - ubernorm; % large penalty if norm is v. different from the median 
    penalty = sum(penalty,2); % sum penalty across components
    
%     figure(2), clf
%     plot([1:Nwin],detrend(ubernorm,'constant'),'o-')
%     figure(3), clf
%     plot([1:Nwin],sum(ubernorm,2),'o-')

    kill = penalty>tolerance*std(penalty);
    if isempty(kill); break; end
    trypenalty = penalty(~kill,:);
    
    if ftest(penalty,1,trypenalty,1) < a_val
        goodwins = goodwins(~kill);
        Nwin = length(goodwins);
        moveon = false; % loop back
    else
        moveon=true; % no reason to remove more windows
    end
end

%debugging only
% figure(99),clf
% subplot(411), semilogx(FZ,dsla_Z(:,goodwins))
% subplot(412), semilogx(FZ,dsla_1(:,goodwins))
% subplot(413), semilogx(FZ,dsla_2(:,goodwins))
% subplot(414), semilogx(FZ,dsla_P(:,goodwins))
%

end

