filename = 'D:\KYBF\TrailPPT\Datas\ÄÆÎ²¼£ÎÂ¶È»ã×Ü.xlsx';
Data = readtable(filename, 'sheet',2);
Alt = table2array(Data(:,2));
Ttrail = table2array(Data(:,3));
Tevmt = table2array(Data(:,7));
Tdlt = table2array(Data(:,8));

FitTrailTemp = fit( Ttrail,Alt, 'exp1');
% FitRatioAVG = polyfit(Ttrail,Alt,2);
% FitTrailTemp = polyval(FitRatioAVG, Ttrail);

scatter(Tevmt,Alt, 50, 'b', 'o','LineWidth', 2)
box on;
grid on;
hold on;
scatter(Ttrail,Alt, 50, 'r', 'o','LineWidth', 2)
hold on;
box on;
grid on;

% for js = 1:size(Tevmt,1)
%     plot([Tevmt(js) Ttrail(js)],[Alt(js) Alt(js)],'--k','linewidth',1.5);
%     hold on;
% end
% 
% plot(FitTrailTemp, Ttrail, Alt);
% % plot(Ttrail,FitTrailTemp,'--m','linewidth',1.5)
ylim([80 100])
xlabel('Temperature (K)');
ylabel('Altitude (km)');
title('Na Meteor Trail Temperature');
set(gca,'fontname','arial','fontsize',12,'fontweight','bold');
legend('T_{env.}','T_{Trail}')


