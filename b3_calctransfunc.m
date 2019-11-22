% calctransfunc
% calculate daily and deployment averaged transfer functions
% operates on all stations for which spectra have been calculated
% Helen Janiszewski 1/2017

clear all; close all

isfigure = 1;
issavefigure = 1;
isoverwrite =1;

%%%%%% DO NOT EDIT BELOW %%%%%%%
setup_parameter;

inpath_dir = sprintf('%s/SPECTRA/',OUTdir);
station_dir = dir(fullfile(inpath_dir));

for itf = 1:length(TF_list)
    TF_size(itf) = length(TF_list{itf});
end
[vec,TF_list_sort] = sort(TF_size,'descend');

TF_list_orig = TF_list;
TF_list_sort_orig = TF_list_sort;

figoutpath=sprintf('%s/STATIONS_TRANSFUNC/',FIGdir);
if ~exist(figoutpath)
    mkdir(figoutpath);
end

for is = 1 : length(station_dir)
    clear TF_matrix
    if isfigure==1
    figure(1);clf
    end
    station = station_dir(is).name;
    inpath = sprintf('%s/SPECTRA/%s/',OUTdir,station);
    inpathav = sprintf('%s/AVG_STA/',OUTdir);
    if ~isdir(inpath)
        continue
    end
    spectra_filenames = dir(fullfile(inpath,['*.mat']));
    if isempty(spectra_filenames)
        continue
    end
    a = length(station);
    disp(station);
    
    
    outpath = sprintf('%s/TRANSFUN/%s/',OUTdir,station);
    if ~exist(outpath)
        mkdir(outpath);
    end
    
for ie = 1 : length(spectra_filenames)+1
    clear cnn_stack cnm_stack
    TF_check = zeros(size(TF_size));
    if ie == length(spectra_filenames)+1
        dayid = 'AVERAGE';
        transfilename = sprintf('%s%s_%s_transfun.mat',outpath,station,dayid);
        if isoverwrite ==0
            if exist(transfilename,'file')
                disp(['Exist: ',transfilename,', Skip!']);
                continue;
            end
        end
    load(sprintf('%s%s_spectraavg.mat',inpathav,station));
    cnn_stack(:,1) = staavg.power.czz_mean;
    cnn_stack(:,2) = staavg.power.c11_mean;
    cnn_stack(:,3) = staavg.power.c22_mean;
    cnn_stack(:,4) = staavg.power.cpp_mean;
    
    cnm_stack(:,1) = staavg.cross.c1z_mean;
    cnm_stack(:,2) = staavg.cross.c2z_mean;
    cnm_stack(:,3) = staavg.cross.cpz_mean;
    cnm_stack(:,4) = staavg.cross.c12_mean;
    cnm_stack(:,5) = staavg.cross.c1p_mean;
    cnm_stack(:,6) = staavg.cross.c2p_mean;
    
    transprop.params.f = avgprop.params.f;
    transprop.params.station = station;
    transprop.params.freqcomp = avgprop.params.freqcomp;
    transprop.params.dayid = dayid;
    transprop.params.NFFT = NFFT;
    
    lineprop = '-';
    cc(ie,:) = [0 0 0];
    
    else
        dayid = spectra_filenames(ie).name((a+2):(a+13));
        transfilename = sprintf('%s%s_%s_transfun.mat',outpath,station,dayid);
        if isoverwrite == 0
            if exist(transfilename,'file')
                disp(['Exist: ',transfilename,', Skip!']);
                continue;
            end
        end
        
    disp(dayid);
    load(fullfile(inpath,spectra_filenames(ie).name));

    if specprop.params.badflag == 1
        disp('Bad day. Skipping')
        continue
    end
    
    f = specprop.params.f;
    elev = specprop.params.elev;
    freqcomp = sqrt(9.8/(2*pi*elev));
        
    cnn_stack(:,1) = specprop.power.czz_stack';
    cnn_stack(:,2) = specprop.power.c11_stack';
    cnn_stack(:,3) = specprop.power.c22_stack';
    cnn_stack(:,4) = specprop.power.cpp_stack';
    
    % the horizontally rotated components
    crr_stack(:,1) = specprop.power.czz_stack';
    crr_stack(:,2) = specprop.rotation.chh_stack';
    crr_stack(:,3) = specprop.rotation.chz_stack';
    crr_stack(:,4) = specprop.rotation.chp_stack';
    crr_stack(:,5) = specprop.power.cpp_stack';
    crr_stack(:,6) = specprop.cross.cpz_stack';
    
    cnm_stack(:,1) = specprop.cross.c1z_stack';
    cnm_stack(:,2) = specprop.cross.c2z_stack';
    cnm_stack(:,3) = specprop.cross.cpz_stack';
    cnm_stack(:,4) = specprop.cross.c12_stack';
    cnm_stack(:,5) = specprop.cross.c1p_stack';
    cnm_stack(:,6) = specprop.cross.c2p_stack';
    
    NFFT = specprop.params.NFFT;
    dt = specprop.params.dt;
    transprop.params.f = f;
    transprop.params.station = specprop.params.station;
    transprop.params.freqcomp = freqcomp;
    transprop.params.dayid = dayid;
    transprop.params.NFFT = specprop.params.NFFT;
    transprop.params.hang = specprop.params.rotor;
    
    lineprop = '-';
    cc(ie,:) = [0.5 0.5 0.5];
    end
    
    ii = 1;
    for itf = 1:length(TF_list_sort)
        TF_cal = TF_list(TF_list_sort(itf));
        if TF_check(TF_list_sort(itf)) == 0 % check if TF calculated already
            % first check if rotational or not
            
            % if rotational
            if ~isempty(strfind(cell2mat(TF_cal),'H'))
                if  ie == length(spectra_filenames)+1 % average doesn't make sense for this
                    continue
                end
                if length(TF_cal{1}) == 2 % 1 component rotational TF, i.e. ZH
                    [lc2c1,label_list] = comp1rotate_calctransfunc(TF_cal,crr_stack);
                    TF_check(TF_list_sort(itf)) = 1;
                    TFs(ii).label = label_list{1};
                    TFs(ii).transfunc = lc2c1;
                    ii = ii+1;
                    TF_matrix(:,ie,(TF_list_sort(itf))) = lc2c1; %for plotting
                elseif length(TF_cal{1}) == 4 % ZP-H
                    [lc2c1,lc2c3,lc3c1_c2,label_list] = comp2rotate_calctransfunc(TF_cal,crr_stack,f);
                    TF_check(TF_list_sort(itf)) = 1;
                    TFs(ii).label = label_list{1};
                    TFs(ii).transfunc = lc2c1;
                    if ~isempty(find(strcmp(label_list{1},TF_list)==1))
                        TFidx = find(strcmp(label_list{1},TF_list)==1);
                        TF_check(TFidx)=1;
                        TF_matrix(:,ie,TFidx) = lc2c1; %for plotting
                    end
                    TFs(ii+1).label = label_list{2};
                    TFs(ii+1).transfunc = lc2c3;
                    if ~isempty(find(strcmp(label_list{2},TF_list)==1))
                        TFidx = find(strcmp(label_list{2},TF_list)==1);
                        TF_check(TFidx)=1;
                        TF_matrix(:,ie,TFidx) = lc2c3; %for plotting
                    end
                    TFs(ii+2).label = label_list{3};
                    TFs(ii+2).transfunc = lc3c1_c2;
                    if ~isempty(find(strcmp(label_list{3},TF_list)==1))
                        TFidx = find(strcmp(label_list{3},TF_list)==1);
                        TF_check(TFidx)=1;
                        TF_matrix(:,ie,TFidx) = lc3c1_c2; %for plotting
                    end
                    ii = ii+3;
                end
                
            else % if component-wise
            if length(TF_cal{1}) == 2 % 1 component TF
                [lc2c1,label_list] = comp1_calctransfunc(TF_cal,cnn_stack,cnm_stack,f);
                TF_check(TF_list_sort(itf)) = 1;
                TFs(ii).label = label_list{1};
                TFs(ii).transfunc = lc2c1;
                ii = ii+1;
                TF_matrix(:,ie,(TF_list_sort(itf))) = lc2c1; %for plotting
            elseif length(TF_cal{1}) == 5 % 3 component TF
                [lc2c1,lc2c3,lc2c4,lc3c4_c2,lc3c1_c2,lc4c1_c3c2,label_list] = comp3_calctransfunc(TF_cal,cnn_stack,cnm_stack,f);
                TF_check(TF_list_sort(itf)) = 1;
                TFs(ii).label = label_list{1};
                TFs(ii).transfunc = lc2c1;
                if ~isempty(find(strcmp(label_list{1},TF_list)==1))
                    TFidx = find(strcmp(label_list{1},TF_list)==1);
                    TF_check(TFidx)=1;
                    TF_matrix(:,ie,TFidx) = lc2c1; %for plotting
                end
                TFs(ii+1).label = label_list{2};
                TFs(ii+1).transfunc = lc2c3;
                if ~isempty(find(strcmp(label_list{2},TF_list)==1))
                    TFidx = find(strcmp(label_list{2},TF_list)==1);
                    TF_check(TFidx)=1;
                    TF_matrix(:,ie,TFidx) = lc2c3; %for plotting
                end
                TFs(ii+2).label = label_list{3};
                TFs(ii+2).transfunc = lc2c4;
                if ~isempty(find(strcmp(label_list{3},TF_list)==1))
                    TFidx = find(strcmp(label_list{3},TF_list)==1);
                    TF_check(TFidx)=1;
                    TF_matrix(:,ie,TFidx) = lc2c4; %for plotting
                end
                TFs(ii+3).label = label_list{4};
                TFs(ii+3).transfunc = lc3c4_c2;
                if ~isempty(find(strcmp(label_list{4},TF_list)==1))
                    TFidx = find(strcmp(label_list{4},TF_list)==1);
                    TF_check(TFidx)=1;
                    TF_matrix(:,ie,TFidx) = lc3c4_c2; %for plotting
                end
                TFs(ii+4).label = label_list{5};
                TFs(ii+4).transfunc = lc3c1_c2;
                if ~isempty(find(strcmp(label_list{5},TF_list)==1))
                    TFidx = find(strcmp(label_list{5},TF_list)==1);
                    TF_check(TFidx)=1;
                    TF_matrix(:,ie,TFidx) = lc3c1_c2; %for plotting
                end
                TFs(ii+5).label = label_list{6};
                TFs(ii+5).transfunc = lc4c1_c3c2;
                if ~isempty(find(strcmp(label_list{6},TF_list)==1))
                    TFidx = find(strcmp(label_list{6},TF_list)==1);
                    TF_check(TFidx)=1;
                    TF_matrix(:,ie,TFidx) = lc4c1_c3c2; %for plotting
                end
                ii = ii+6;
            end
            end
            
        end
        
        
    end
    
  filename = [outpath,station,'_',dayid,'_transfun.mat'];
        save(filename,'TFs','transprop')
        clear TFs
end
if isfigure==1 && isoverwrite==1
figure(1)
for itf = 1:length(TF_list_sort)
subplot(length(TF_list_sort),1,itf)
set(gca,'ColorOrder',cc,'NextPlot','replacechildren');
loglog(f,abs(TF_matrix(:,:,itf)),'-');
title([station,'Transfer Function',TF_list{itf}]); xlim([10^-4 max(f)]); %ylim([y1/10, y2*10])
set(gca,'yscale','log','xscale','log');
end
figure(1)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);
end

if isfigure==1 && issavefigure==1 && isoverwrite==1
filename=sprintf('%s/%s_tfs',figoutpath,station);
print(gcf,'-dpng',filename)
end

end
