%% 尾迹变化图
TimeCap = TimeXCap(:,2347:2349)';
TimeCa = TimeXCa(:,2354:2356)';
TimeNumCap = datenum(TimeCap);
TimeNumCa = datenum(TimeCa);
DenTrailCap = [111 106 111 111 90 26 13];
DenTrailCa = [160 366 594 766 773 419 74];
% DenTrailCap = [40 60 34];
% DenTrailCa = [78 141 100];
Profx = 2079;
SEC = 1/(24*60*60);
TimeIntCap = (TimeNumCap(1,:):SEC:TimeNumCap(end,:))';
TimeIntCa = (TimeNumCa(1,:):SEC:TimeNumCa(end,:))';
FitCap = interp1(TimeNumCap,DenTrailCap,TimeIntCap,'linear');
FitCa = interp1(TimeNumCa,DenTrailCa,TimeIntCa,'linear');

% 检查图像
figure('Name','Check')
plot(TimeIntCa,FitCa,'-b','linewidth',1.5)
hold on;
plot(TimeIntCap,FitCap,'--b','linewidth',1.5)
datetick('x','HH:MM:SS')
% ylim([0 600])
ylabel('Density of Trail (cm^{-3})')
grid on
TDstr = 'starttime: '+string(Profx)+' | '+TimeList(:,Profx);
title(TDstr)
xlabel('Time')
datetick('x','HH:MM:SS','keeplimits','keepticks');
legend('Original Ca','Original Ca^+','location','northwest')

%% 修改与裁剪

% Ca 左

% StdT = round((TimeIntCap(1,:)-TimeIntCa(1,:))./SEC);
% EndT = round((TimeIntCap(end,:)-TimeIntCa(end,:))./SEC);
% 
% IntCap = FitCap(1:end-EndT,:);
% IntCa = FitCa(StdT+1:end,:);
% TimeCut = TimeIntCap(1:end-EndT,:);

% Cap 左

% StdT = round((TimeIntCa(1,:)-TimeIntCap(1,:))./SEC);
% EndT = round((TimeIntCa(end,:)-TimeIntCap(end,:))./SEC);
% 
% IntCa = FitCa(1:end-EndT,:);
% IntCap = FitCap(StdT+1:end,:);
% TimeCut = TimeIntCap(StdT+1:end,:);

% 定制

StdT = round((TimeIntCa(1,:)-TimeIntCap(1,:))./SEC);
EndT = round((TimeIntCap(end,:)-TimeIntCa(end,:))./SEC);

IntCa = FitCa(1:end,:);
IntCap = FitCap(StdT+1:end-EndT,:);
TimeCut = TimeIntCap(StdT+1:end-EndT,:);

SunCaCap = IntCap+IntCa;
RatioCap = IntCap./SunCaCap;
RatioCa = IntCa./SunCaCap;
RatioCaCap = IntCa./IntCap;

RatioSum = zeros(size(RatioCaCap,1)-1,1);
for jr = 2:size(RatioCaCap,1)
    RatioSum(jr-1,:) = RatioCaCap(jr,:)-RatioCaCap(jr-1,:);
end
RatioCaAVG = mean(RatioSum)
RatioCaSUM = sum(RatioSum)


%% 绘图
figure('Name','Origin','position',[50,100,1450,735])
subplot(2,2,1)
plot(TimeIntCa,FitCa,'-b','linewidth',1.5)
hold on;
plot(TimeIntCap,FitCap,'--b','linewidth',1.5)
datetick('x','HH:MM:SS')
% ylim([0 600])
ylabel('Density of Trail (cm^{-3})')
grid on
TDstr = 'starttime: '+string(Profx)+' | '+TimeList(:,Profx);
title(TDstr)
xlabel('Time')
datetick('x','HH:MM:SS','keeplimits','keepticks');
legend('Original Ca','Original Ca^+','location','northwest')

% figure('Name','Cut')
subplot(2,2,2)
plot(TimeCut,IntCa,'-b','linewidth',1.5)
hold on;
plot(TimeCut,IntCap,'--b','linewidth',1.5)
datetick('x','HH:MM:SS')
% ylim([0 600])
ylabel('Density of Trail (cm^{-3})')
grid on
TDstr = 'starttime: '+string(Profx)+' | '+TimeList(:,Profx);
title(TDstr)
xlabel('Time')
datetick('x','HH:MM:SS','keeplimits','keepticks');
legend('Fit Ca','Fit Ca^+','location','northwest')

% figure('Name','Add')
subplot(2,2,3)
plot(TimeCut,SunCaCap,'-r','linewidth',1.5)
datetick('x','HH:MM:SS')
% ylim([0 600])
ylabel('Density of Trail (cm^{-3})')
grid on
TDstr = 'starttime: '+string(Profx)+' | '+TimeList(:,Profx);
title(TDstr)
xlabel('Time')
datetick('x','HH:MM:SS','keeplimits','keepticks');
legend('Sum Density','location','northwest')

% figure('Name','Ratio')
subplot(2,2,4)
plot(TimeCut,RatioCaCap,'-b','linewidth',1.5)
% hold on;
% plot(TimeCut,RatioCaCap.*100,'--b','linewidth',1.5)
datetick('x','HH:MM:SS')
% ylim([0 600])
ylabel('Ratio of Trail')
grid on
TDstr = 'starttime: '+string(Profx)+' | '+TimeList(:,Profx);
title(TDstr)
xlabel('Time')
datetick('x','HH:MM:SS','keeplimits','keepticks');
legend('Ratio of Ca/Ca ion','location','northwest')


%% 直接对齐

Xiexls = 'E:\TrailIMG\CaCapRatioXie.xlsx';
Table = readtable(Xiexls,'sheet','Calcu');
TabDat = table2array(Table);
TimeS = TabDat(:,1);
DenCa = TabDat(:,2);
DenCap = TabDat(:,3);

TimeCap2 = TimeXCap(:,423:425)';
% TimeCa2 = TimeXCa(:,2085:2901)';
% TimeNumCap = datenum(TimeCap);
% TimeNumCa = datenum(TimeCa);
DenTrailCap2 = DenCap';
DenTrailCa2 = DenCa';
Profx = 423;

RatioCaCap2 = DenTrailCa2./DenTrailCap2;
SunCaCap2 = DenTrailCa2 + DenTrailCap2;

% RatioCaAVG2 = (RatioCaCap2(end)-RatioCaCap2(1))./((size(RatioCaCap2,2)-1).*66)
RatioCaSUM2 = RatioCaCap2(end)-RatioCaCap2(1)
RatioCaAVG2 = RatioCaSUM2./(11)

figure('Name','Origin','position',[50,100,1450,735])
% subplot(2,2,1)
% plot(TimeIntCa,FitCa,'-b','linewidth',1.5)
% hold on;
% plot(TimeIntCap,FitCap,'--b','linewidth',1.5)
% datetick('x','HH:MM:SS')
% % ylim([0 600])
% ylabel('Density of Trail (cm^{-3})')
% grid on
% TDstr = 'starttime: '+string(Profx)+' | '+TimeList(:,Profx);
% title(TDstr)
% xlabel('Time (UT)')
% datetick('x','HH:MM:SS','keeplimits','keepticks');
% legend('Original Ca','Original Ca^+','location','northwest')

% figure('Name','Cut')
subplot(2,2,1)
plot(TimeCap2,DenTrailCa2,'-b','linewidth',1.5)
hold on;
plot(TimeCap2,DenTrailCap2,'--b','linewidth',1.5)
hold on;
plot(TimeCap2,SunCaCap2,'-r','linewidth',1.5)
datetick('x','HH:MM:SS')
% ylim([0 600])
ylabel('Density of Trail (cm^{-3})')
grid on
TDstr = "(a) "+TimeList(:,Profx);
title(TDstr)
xlabel('Time (UT)')
datetick('x','HH:MM:SS','keeplimits','keepticks');
legend('Ca','Ca^+','Sum Density','location','northwest')
set(gca,'fontsize',11);

% figure('Name','Add')

subplot(2,2,3)
imgtrail = "D:\File\20241118ppt文章\ProfCaCap-1.jpg";
imshow(imgtrail);
axis tight; % 确保图像铺满子图范围

% subplot(2,2,2)
% 
% plot(TimeCap2,RatioCaCap2,'-r','linewidth',1.5)
% box on; grid on;
% datetick('x','HH:MM:SS')
% ylabel('Ratio of Trail')
% grid on
% TDstr = '(b) starttime: '+string(Profx)+' | '+TimeList(:,Profx);
% title(TDstr)
% xlabel('Time (UT)')
% datetick('x','HH:MM:SS','keeplimits','keepticks');
% legend('Sum Density','Ratio of Ca/Ca ion','location','northwest')

subplot(2,2,2)

RatiAB = DenCa./DenCap;
RatiABB = DenCa./(DenCa+DenCap);
scatter(TimeS,RatiAB,50,'r','linewidth',2)
hold on; box on; grid on;
customModel = fittype('a*x+b',...
    'independent', 'x', 'coefficients', {'a', 'b'});
XieFit = fit(TimeS,RatiAB,customModel);
coeffs1 = coeffvalues(XieFit);
Slope1 = coeffs1(1);Slope1 = round(Slope1,4);
Xmix1 = XieFit(TimeS);
plot(TimeS,Xmix1,'--k','linewidth',2)
TDstr = '(b) Slope: '+string(Slope1)+' | '+TimeList(:,Profx);
title(TDstr)
ylabel('Ratio');
set(gca,'fontsize',11);
legend('Ratio of Ca/Ca^+','Fit curve');

subplot(2,2,4)
scatter(TimeS,RatiABB,50,'r','linewidth',2)
hold on; box on; grid on;
XieFit2 = fit(TimeS,RatiABB,customModel);
coeffs2 = coeffvalues(XieFit2);
Slope2 = coeffs2(1);Slope2 = round(Slope2,4);
Xmix2 = XieFit2(TimeS);
plot(TimeS,Xmix2,'--k','linewidth',2)
TDstr = '(d) Slope: '+string(Slope2)+' | '+TimeList(:,Profx);
title(TDstr)
ylabel('Ratio');
set(gca,'fontsize',11);
legend('Ratio of Ca+/Ca+Ca^+','Fit curve');







