%% 创建示例数据
x = [75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105];
y = [0  0  0  0  0  0  0  1  7  9  13  16  8  8  8  8  6  6  4  5  4  2  0  1  1  0   0   0   0   0   0];
% y = [0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0];

% 绘制横向柱状图
barh(x, y, 'b'); % 'b' 表示蓝色
xlabel('流星尾迹数');
ylabel('Altitude (km)');
title('Meteoric Trails in Ca Channel');
xlim([0 20]);
% hold on;

%% 创建示例数据
x = [75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 125];
y = [0   0   0    0   0    0   0   99   41   426   653    992   1668   1934   2744    2342   1728   1792   766   505    0   61   0   2391   291    88     0      52     0      0     0     0];
% Cay = [0  0  0  0  0  0  0  19  319  531  886  2182  515  1158  811  717  462  652  223  519  288  123  0  10  34  0   0   0   0   0   0  0];
% y = [0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0];

% 绘制横向柱状图
yyaxis left;
plot(x, y , '-b', 'linewidth', 2); % 'b' 表示蓝色
xlabel('Altitude (km)');
ylabel('Ca ion Trail Density (cm^{-3})');
title('Meteoric Trails Density in Ca ion Channel');
xlim([75 125]);

% 平均密度

xls1 = xlsread('C:\Users\18333\Desktop\近期汇报\漠河流星尾迹\CaCap.xlsx');
Alt = xls1(2:end, 1);
DenAll = xls1(2:end, 2:end);
SumAll = sum(DenAll,2);
MeanAll = SumAll./ (613+478+431+581-4);
yyaxis right;
plot(Alt, MeanAll , '-r', 'linewidth', 2);
xlabel('Altitude (km)');
ylabel('Ca ion Density (cm^{-3})');
xlim([75 125]);

%% 分布特征：钙原子离子化
xls1 = xlsread('C:\Users\18333\Desktop\四月\漠河流星尾迹\CaCap.xlsx');
Alt = xls1(:,2);
CapTiDen = xls1(:,3);
CaTiDen = xls1(:,4);
SumTiDen = xls1(:,5);
RCP2C = xls1(:,6);
Rion = xls1(:,7);
scatter(Alt, Rion,SumTiDen, 'linewidth', 2);
xlim([80 105])
box on;
grid on;

% title('Ca^{+} / Ca Trail Density Ratio');
title('Calcium Ionization Ratio');

xlabel('Altitude (km)')
ylabel('Ratio (%)')

% hold on;
% plot([80 105],[1 1],'-r','linewidth',2)

%% 相关性计算

dataMatrix = [Alt, SumTiDen, Rion];
correlation_matrix = corr(dataMatrix, 'Type', 'Kendall');
disp(correlation_matrix);















