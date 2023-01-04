function plot_tilt(dayid,day_deploy,ang,coh)
%
% Plot angle and coherence versus days since deployment.
% Yuechu Wu
% 12131066@mail.sustech.edu.cn
% 2022-12-12

c = colormap('jet');

if isleapyear(str2double(day_deploy(1:4)))
    year = 366;
else
    year = 365;
end

days = dayofyear(str2double(dayid(1:4)),str2double(dayid(5:6)),str2double(dayid(7:8)),0,0);
cc = interp1(1:length(c),c,((days)/(year))*(length(c)-1)+1);

deploynum = datenum(day_deploy(1:8),'yyyymmdd');
deploydate = datestr(deploynum,'yyyy-mmm-dd');
iday = datenum(dayid(1:8),'yyyymmdd') - deploynum + 1;

figure(1)
subplot(211);hold on
plot(iday,ang,'o','MarkerFaceColor',cc,'MarkerEdgeColor','none','MarkerSize',5);
ylabel('Tilt angle (°)','FontSize',10);

subplot(212);hold on
plot(iday,coh,'o','MarkerFaceColor',cc,'MarkerEdgeColor','none','MarkerSize',5);
ylim([0 1]);
set(gca,'YTick',0:0.2:1,'FontSize',10);
ylabel('Coherence','FontSize',10);
xlabel(sprintf('%s (%s)','Days since deployment',deploydate),'FontSize',10);

return