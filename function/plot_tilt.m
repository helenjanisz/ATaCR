function plot_tilt(idays,day_deploy,angs,cohs,isfigure_orient)

% Plot tilt angle and coherence versus days since deployment.
% Yuechu Wu
% 12131066@mail.sustech.edu.cn
% 2022-12-12

if nargin == 4
    isfigure_orient = 0;
end

c = colormap('jet');
ccs = NaN(length(idays),1);

deploynum = datenum(day_deploy(1:8),'yyyymmdd');

if isleapyear(str2double(day_deploy(1:4)))
    year = 366;
else
    year = 365;
end

for i = 1:length(idays)
    iday = idays(i);
    dayid = datestr(deploynum + iday - 1,'yyyymmdd');
    jday = dayofyear(str2double(dayid(1:4)),str2double(dayid(5:6)),str2double(dayid(7:8)),0,0);
    cc = interp1(1:length(c),c,((jday)/(year))*(length(c)-1)+1);
    ccs(i,1:3) = cc;
end

if isfigure_orient
    figure(3)
else
    figure(1)
end

for id = 1:length(idays)
    subplot(211)
    plot(idays(id),angs(id),'o','MarkerFaceColor',ccs(id,:),'MarkerEdgeColor','none','MarkerSize',5);
    hold on

    subplot(212)
    plot(idays(id),cohs(id),'o','MarkerFaceColor',ccs(id,:),'MarkerEdgeColor','none','MarkerSize',5);
    hold on
end


return