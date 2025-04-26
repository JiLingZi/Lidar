%{
AltCaCap = [92.34,   86.14,   86.78,   89.55,   87.71,   ...
            87.64,   89.95,   89.92,   85.98,   84.82,   ...
            93.17,   95.02,   87.64,   86.23,   89.27,   81.56];
DenCa =    [3152,    529,     1584,    604,     863,     ...
            202,     319,     450,     532,     55,      ...
            42,      59,      1310,    394,     813,     342];
DenCap =   [650,     120,     433,     310,     61,      ...
            294,     134,     48,      13,      70,      ...
            133,     183,     1148,    76,      191,     59];
%}

XlsName ='D:\KYBF\TrailPPT\Datas\CaCapSumDenAlt.xlsx';
Table = readtable(XlsName,'sheet','Sheet1');
TabDat = table2array(Table);
AltCaCap = TabDat(:,1);
DenCa = TabDat(:,2);
DenCap = TabDat(:,3);

RatioCaCap = (DenCa ./ DenCap);
AddCaCap = DenCa + DenCap;

FitMix = [RatioCaCap,AltCaCap];
ShunMix = sortrows(FitMix,2);

Fexp = fit(ShunMix(:,2),ShunMix(:,1), 'exp1');

customModel = fittype('a/x+b',...
    'independent', 'x', 'coefficients', {'a', 'b'});

fr = fit(ShunMix(:,1),ShunMix(:,2), customModel);

yFit = feval(fr, ShunMix(:,1));
% fr = polyfit(ShunMix(:,2),ShunMix(:,1), 3);

% 卡方
x = ShunMix(:,1);
y = ShunMix(:,2);
yFit = feval(fr, x);% 计算预测值
residuals = y - yFit;% 计算残差
chiSquare = sum((residuals.^2) ./ yFit);% 计算卡方
fprintf('卡方值: %f\n', chiSquare);% 显示卡方值

% R方
yFit = feval(fr, x);
% 计算残差平方和 (RSS)
residuals = y - yFit;
RSS = sum(residuals.^2);
TSS = sum((y - mean(y)).^2);% 计算总平方和 (TSS)
R_square = 1 - (RSS / TSS);% 计算R平方值
fprintf('R平方值: %f\n', R_square);% 显示R平方值


% frmix = 1140*log10(ShunMix(:,2))./ShunMix(:,2)-22.16;
fitx = (0.1:0.01:50)';
frmix = 1.083./fitx + 85.92;

% cftool(ShunMix(:,2),ShunMix(:,1))

figure('name','合比')
scatter(RatioCaCap,AltCaCap,50,'b','linewidth',2);
grid on;box on;hold on;
plot(fitx,frmix,'--k','linewidth',2)
title('Sum Density Ratio (Ca/Ca ion)')
xlabel('Ratio');
ylabel('Altitude (km)');
xlim([-5 50]);
% xlim([0.01 100])
% set(gca, 'XScale', 'log');
hold on;
plot([1 1 ],[75 105],'-k','linewidth',2)
% xticks([0.01 0.1 1 10 100]);
% xticklabels({'0.01', '0.1', '1', '10', '100'});
legend('Ratio','Fit Curve')
set(gca,'fontsize',12)



%%
FitX = 80:0.1:98.5;
px =7.12;
frmix = -0.001447*(FitX-px).^3+...
    0.3891*(FitX-px).^2-35.05*(FitX-px)+1058;
% frmix = 1140*(log10(0.9.*FitX))./(1.0.*FitX)-22.16;
log10(0.9)
figure('name','合比')
scatter(RatioCaCap,AltCaCap,50,'b','linewidth',2);
grid on;box on;hold on;
plot(frmix,FitX,'-r','linewidth',2)
title('Sum Density Ratio (Ca/Ca ion)')
xlabel('Ratio');
ylabel('Altitude (km)');
xlim([-5 50]);
% xlim([0.01 100])
% set(gca, 'XScale', 'log');
hold on;
plot([1 1 ],[75 105],'-k','linewidth',2)
% xticks([0.01 0.1 1 10 100]);
% xticklabels({'0.01', '0.1', '1', '10', '100'});
legend('Ratio of Ca/Ca^+','Fit Curve')

%%

figure('name','分颜色')
RxCa = DenCa./(DenCa+DenCap);
RxCap = DenCap./(DenCa+DenCap);
FitRxCa = fit(RxCa,AltCaCap, 'exp1');
FitRxCap = fit(RxCap,AltCaCap, 'exp1');
MixX = 1e-2:0.001:1;
MixRxCa = FitRxCa(MixX);
MixRxCap = FitRxCap(MixX);
scatter(RxCa,AltCaCap,50,'b','linewidth',2);
grid on;box on;hold on;
plot(MixX,MixRxCa,'--k','linewidth',1.5);
% hold on;
% scatter(RxCap,AltCaCap,50,'r','linewidth',2);
% hold on;grid on;box on;
% plot(MixX,MixRxCap,'--k','linewidth',1.5);
title('The altitude distribution of Ca/(Ca+Ca^+) Density')
xlabel('Ratio');
ylabel('Altitude (km)');
xticks([0.1 1]);
xticklabels({'0.1', '1'});
xlim([0.08 1])
ylim([80 100])
legend('Ca/(Ca+Ca^+)','Fit curve','location','northwest')
set(gca, 'XScale', 'log','fontsize',12);

%% 月份分布

startDate = datetime(2023, 10, 1);
endDate = datetime(2024, 9, 1);
stlim = datetime(2023, 9, 15);
edlim = datetime(2024, 9, 15);
Yue = startDate:calmonths(1):endDate;

XlsNameN ='D:\KYBF\TrailPPT\Datas\时长数量柱状图.xlsx';
TableN = readtable(XlsNameN,'sheet','CaCap');
AryN = table2array(TableN);
Num = AryN(:,2);%个数
bar(Yue,Num,0.7,'b','linewidth',1)
grid on;
box on;
xlim([stlim,edlim]);
datetick('x','mmm','keeplimits')
xlabel('Time');
ylabel('Simultaneous Meteor Trail Number')
title('Simultaneous occurrences of Ca and Ca^{+} trails')
set(gca, 'fontsize',12);

%% 斜率先增后减
AltXx = [92.34   86.11   86.78   89.95   89.95   87.74    ];
% Slope = [0.14755 0.17635 0.01679 0.04956 0.01803 -0.01817 ];
Slope = [0.0042  0.0033  -0.002 0.0043  0.0003   -0.0035  ];

TableSl = readtable('F:\TrailIMG\CaCapRatioXie.xlsx','sheet','Xie');
AltXx = table2array(TableSl(:,1));
Slope = table2array(TableSl(:,2));

scatter(Slope,AltXx,50,'r','linewidth',2);
grid on; box on; hold on;

ftsp = polyfit(AltXx,Slope, 1);
x_fit = linspace(min(AltXx), max(AltXx), 100);
y_fit = polyval(ftsp, x_fit);
plot(y_fit,x_fit,'--k','linewidth',2)

hold on;
plot([0 0],[80 100],'-k','linewidth',2);
% xlim([-5e-3 5e-3]);
% xticks([-5e-3,-4e-3,-3e-3,-2e-3,-1e-3,0e-3,1e-3,2e-3,3e-3,4e-3,5e-3,])
xlabel('Slope');
ylabel('Altitude (km)');
set(gca,'fontsize',10)
title('Slope of Ca^+/(Ca+Ca^+) ratio variation in complete meteor trails')
legend('Slope','Fit curve')
set(gca,'fontsize',10)
% cftool(AltXx,Slope)


%% 钠尾迹温度高度关系

AltNas = [105.25, 84.20, 85.45, 83.06];
TemNas = [365,    272,   317,   250  ];
scatter(TemNas,AltNas,50,'r','linewidth',2);
hold on; box on; grid on;
ylim([80 110]);
xlabel('Temperature (K)');
ylabel('Altitude (km)');
title('Temperature of Na Meteor Trail');

%% 三维图

xls3D = "D:\KYBF\TrailPPT\Datas\CaCap高度密度时长三维.xlsx";
Tab3D = readtable(xls3D,'sheet','Sheet1');
Alt3D = table2array(Tab3D(:,1));
DTm3D = table2array(Tab3D(:,2));
Den3D = table2array(Tab3D(:,3));
scatter3(DTm3D, Den3D, Alt3D, 36, 'b', 'filled');
xlabel('Duration (s)');
ylabel('Density (cm^{-3})');
zlabel('Altitude (km)');
grid on; box on; hold on;


%% 
T = 150:1:350;
R = 3.8.*(10e-12).*(T/200).^(-2.5);
plot(T,R)

%% 四成分中位高度与CAMOD对比
NaAlt = [87.48, 94.0, 93.5];
K_Alt = [86.41, 93.8, 93.2];
FeAlt = [86.63, 85.5, 82.4];
CaAlt = [86.56, 82.3, 81.1];

Xtk = [-5; 0; 5];

plot(Xtk,NaAlt,'-r','linewidth',2);
hold on; grid on; box on;
plot(Xtk,K_Alt,'-b','linewidth',2);
hold on; grid on; box on;
plot(Xtk,FeAlt,'-g','linewidth',2);
hold on; grid on; box on; 
plot(Xtk,CaAlt,'-m','linewidth',2);
hold on; grid on; box on;

xticks([-8 -5 0 5 8]);
xticklabels({' ', 'This work', 'CAMOD 37°','CAMOD 0°', ' '})
xlim([-8 8]); ylim([80 100])
ylabel('Altitude (km)');









