filename = 'D:\KYBF\TrailPPT\Datas\CaCapRatio.xlsx';
Data = readtable(filename, 'sheet',2);
Alt = table2array(Data(:,2));
RatioAVG = table2array(Data(:,3))./100;
RatioSUM = table2array(Data(:,4))./100;

FitRatioAVG = polyfit(Alt,RatioAVG,1);
FitAVG = polyval(FitRatioAVG, Alt);
FitRatioSUM = polyfit(Alt,RatioSUM,1);
FitSUM = polyval(FitRatioSUM, Alt);

figure('name','AVG')
plot(FitAVG,Alt,'--m','linewidth',1.5)
hold on;
scatter(RatioAVG,Alt, 50, 'b', 'o','LineWidth', 2)
box on;
grid on;
hold on;
for js = 1:size(RatioAVG,1)
    plot([RatioAVG(js) 0],[Alt(js) Alt(js)],'--k','linewidth',1.5);
    hold on;
end
plot([0 0],[80 96],'-k','linewidth',2)
title('Mean Ratio Variation of Ca')
xlabel('Ratio')
ylabel('Altitude (km)')
legend('fitted curve')
set(gca,'fontname','arial','fontsize',12,'fontweight','bold')

figure('name','SUM')
plot(FitSUM,Alt,'--m','linewidth',1.5)
hold on;
scatter(RatioSUM,Alt, 50, 'b', 'o','LineWidth', 2)
box on;
grid on;
hold on;
for js = 1:size(RatioSUM,1)
    plot([RatioSUM(js) 0],[Alt(js) Alt(js)],'--k','linewidth',1.5);
    hold on;
end
plot([0 0],[80 96],'-k','linewidth',2)
title('Sum Ratio Variation of Ca')
xlabel('Ratio ')
ylabel('Altitude (km)')
legend('fitted curve')
set(gca,'fontname','arial','fontsize',12,'fontweight','bold')

