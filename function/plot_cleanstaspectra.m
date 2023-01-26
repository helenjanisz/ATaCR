function plot_cleanstaspectra(spect,coh_stack,ph_stack,ad_stack,cc,station,f,maxpow,minpow,comporder,plotorder)

all_comporder = {'Z','H1','H2','P'};
all_plotorder = {'1Z','2Z','HZ','PZ','12','P1','P2'};

for ico = 1:length(comporder)
    co_id(ico) = find(ismember(all_comporder,comporder(ico)));
end

for ipo = 1:length(plotorder)
    po_id(ipo) = find(ismember(all_plotorder,plotorder(ipo)));
end


% Plotting Power Spectra
figure(1)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);

for ip = 1:length(comporder)
    subplot(length(comporder),1,ip)
    set(gca,'ColorOrder',cc,'NextPlot','replacechildren');
    loglog(f,spect(:,:,co_id(ip)),'-','LineWidth',.5);
%     hold on
    title(sprintf('%s-component, Station: %s',comporder{ip},station));
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    set(gca,'yscale','log','xscale','log');
    xlim([10e-4 max(f)]); 
    if isnan(minpow(co_id(ip)))
        ylim([0 1])
    elseif minpow(co_id(ip))==maxpow(co_id(ip))
        ylim([0 1])
    else
        ylim([minpow(co_id(ip)) maxpow(co_id(ip))]);
    end
    box on
end

figure(2)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);

figure(3)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);

figure(4)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);

for ip = 1:length(plotorder)
    figure(2)
    subplot(length(plotorder),1,ip)
    set(gca,'ColorOrder',cc,'NextPlot','replacechildren');
    semilogx(f',coh_stack(:,:,po_id(ip)),'-','LineWidth',.5); %,('Color',cc);
    hold on
    title(sprintf('%s Coherence: %s',station, plotorder{ip}))
    xlabel('Frequency (Hz)')
    ylabel('Coherence')
    ylim([0 1]); xlim([10^-4 max(f)]);
    set(gca,'xscale','log');
    box on
    
    figure(3)
    subplot(length(plotorder),1,ip)
    set(gca,'ColorOrder',cc,'NextPlot','replacechildren');
    semilogx(f',ph_stack(:,:,po_id(ip)),'o','MarkerSize',1);
    hold on
    title(sprintf('%s Phase: %s',station,plotorder{ip}))
    xlabel('Frequency (Hz)')
    ylabel('Phase')
    xlim([10^-4 max(f)]);
    set(gca,'xscale','log');
    box on
    
    figure(4)
    subplot(length(plotorder),1,ip)
    set(gca,'ColorOrder',cc,'NextPlot','replacechildren');
    loglog(f',ad_stack(:,:,po_id(ip)),'-','LineWidth',.5)
    hold on
    title(sprintf('%s Admittance: %s',station,plotorder{ip}))
    xlabel('Frequency (Hz)')
    ylabel('Admittance')
    xlim([10^-4 max(f)]);
    set(gca,'yscale','log','xscale','log');
    box on
    
end

return