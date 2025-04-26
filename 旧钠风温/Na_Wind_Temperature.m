%-------------------------------------------------------------------------
%=================中高层大气钠风温激光雷达风温反演程序======================
%-------------------------------------------------------------------------
% 项目名称：漠河钠原子风温反演程序
% 作者名称：王可鑫
% 作者单位：江苏科技大学
% 作者导师：王泽龙
% 作者邮箱：wang.kexin@stu.just.edu.cn
% 创建时间：2023-07-09 17:30(UTC+8)
% 修改时间：2023-07-10 17:30(UTC+8)
% 数据来源：子午工程二期漠河激光雷达台站 
%===============================注意事项==================================
% 切记不要使用 clear all；浪费时间
% 高距离门：bin_num = 8912
% 时间分辨率：2 min
% 空间分辨率：30 m
% 激光线宽：80 MHz
% 三频切换频移量:585 MHz
% V：垂直方向
% N：北方向，天顶倾角30°
% W：东方向（标记错误），天顶倾角30°
% 大气模型日期快速计算：
%                     闰年：03月+060；04月+091；05月+121；06月+152；
% % % % 01月+000；02月+031；03月+059；04月+090；05月+120；06月+151；
% 闰年：07月+182；08月+213；09月+244；10月+274；11月+305；12月+335。
% % % % 07月+181；08月+212；09月+243；10月+273；11月+304；12月+334。
%-------------------------------------------------------------------------

% 加载三频有效后向散射截面矩阵
load Data_eff_0;
load Data_eff_R;
load Data_eff_L;

% 温度数据
T = [80 100 120 140 160 180 200 220 240 260 280 300 320]';
T = repmat(T,1,size(T,1));

% 风速数据
V = [-120 -100 -80 -60 -40 -20 0 20 40 60 80 100 120];
V = repmat(V,size(V,2),1);

% 风温比率矩阵
R_T = (Data_eff_R + Data_eff_L) ./ (2 * Data_eff_0);
R_V = (Data_eff_R - Data_eff_L) ./ (Data_eff_0);

% 对温度进行二维曲面拟合
ft = fittype( 'poly55' );
[fit_T, gof_T] = fit( [R_V(:), R_T(:)], T(:), ft );
coe_fit_t=coeffvalues(fit_T);
f_T = fit_T;

% 画出温度拟合曲面
figure(1);
plot(f_T);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Fitting Surface of Temperature');
xlabel('R_{V}');
ylabel('R_{T}');
zlabel('Temperature (K)');

% 对风速进行二维曲面拟合
ft = fittype( 'poly55' );
[fit_V, gof_V] = fit( [R_V(:), R_T(:)], V(:), ft );
coe_fit_v=coeffvalues(fit_V);
f_V = fit_V;

% 画出风速拟合曲面
figure(2);
plot(f_V);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Fitting Surface of Wind Velocity');
xlabel('R_{V}');
ylabel('R_{T}');
zlabel('Wind Velocity (m/s)');

%%
% 读取钠风温雷达采集数据矩阵
% Na_table = readtable('D:\ZWDATA\Na20230707\20230707\Na\Na20230707171542258.dat');
% Na_table = readtable('D:\ZWDATA\20230709Na\Na\Na20230709060144901.dat');
Na_table = readtable('G:\Mohe_Data\589Na\Na20230709150515310.dat');%首次反演出的日期
Na_data = table2array(Na_table);
height_num = Na_data(:,1);
Na_photon = Na_data(:,2:10);
Na_photon = Na_photon - mean(Na_photon(4232:4558,:),1);
Na_photon = Na_photon ./ max(Na_photon);

% 每1km合并光子数行
Na_i_read = [];
for i = 1:256
    Na_i_read(i,:) = sum(Na_photon(1+32*(i-1):32*i,:),1);
end
Na_photon = Na_i_read;
height = 1:256;
height = height * 0.9830390625;
height = height(1:221);

% 对不同方向不同频率进行拆分
V_F_0 = Na_photon(:,1);
V_F_R = Na_photon(:,2);
V_F_L = Na_photon(:,3);
N_F_0 = (Na_photon(:,4))*(sqrt(3)/2);
N_F_R = (Na_photon(:,5))*(sqrt(3)/2);
N_F_L = (Na_photon(:,6))*(sqrt(3)/2);
W_F_0 = (Na_photon(:,7))*(sqrt(3)/2);
W_F_R = (Na_photon(:,8))*(sqrt(3)/2);
W_F_L = (Na_photon(:,9))*(sqrt(3)/2);

% 计算RV垂直方向、北向、东向
R_V_real_V = (V_F_R - V_F_L) ./ (V_F_0);  % 铅垂RV
R_V_real_N = (N_F_R - N_F_L) ./ (N_F_0);  % 北向RV
R_V_real_W = (W_F_R - W_F_L) ./ (W_F_0);  % 东向RV

% 计算RT
R_T_real_V = (V_F_R + V_F_L) ./ (2 * V_F_0);
R_T_real_N = (N_F_R + N_F_L) ./ (2 * N_F_0);
R_T_real_W = (W_F_R + W_F_L) ./ (2 * W_F_0);

% 计算垂直风速
V_real_V = f_V(R_V_real_V,R_T_real_V);
V_real_V = V_real_V(1:221,:);

% 计算北向风速,并乘以sin（30°）
V_real_N = zeros(221,1);
for i = 1:221
    V_real_N(i,:) = f_V(R_V_real_N(round(i*(2/sqrt(3))),:),R_T_real_N(round(i*(2/sqrt(3))),:));
end
V_real_N = V_real_N * 0.5;

% 计算东向风速,并乘以sin（30°）
V_real_W = zeros(221,1);
for i = 1:221
    V_real_W(i,:) = f_V(R_V_real_W(round(i*(2/sqrt(3))),:),R_T_real_W(round(i*(2/sqrt(3))),:));
end
V_real_W = V_real_W * 0.5;

% 计算风速
V_real_Array = [V_real_N,V_real_W,V_real_V];
V_real_length = zeros(221,1);
V_real_direction = zeros(221,3);
for i = 1:221
    V_real_length(i,:) = norm(V_real_Array(i,:));
    V_real_direction(i,:) = V_real_Array(i,:) / norm(V_real_Array(i,:));
end

% 计算温度
T_real = f_T(R_V_real_V,R_T_real_V);
T_real = T_real(1:221,:);

% % 东向与北向高度修正因子，铅垂方向必须关闭
% height = height * (sqrt(3)/2);

% 画温度图
yyyy = 2023;
dd = 181+9;
hh = 15;
[T_model,Den_model] = atmosnrlmsise00([height*1e3],53.5,122.3,yyyy,dd,hh);
T_model = T_model(:,2);
figure('Name','钠原子温度反演图');
plot(T_real,height,'-b','linewidth',2);
hold on
plot(T_model,height,'-r','linewidth',2);
legend('T_{Mohe}', 'T_{Msise-00}', 'Location', 'northwest');
set(legend, 'EdgeColor', 'w');
xlim([120 320]);
ylim([80 105]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Temperature 20230709 23:05');
xlabel('Temperature (K)');
ylabel('Altitude (km)');

% 画风速图,V_real_length是三向合风速，只有正值
Meridional_wind = V_real_N;     % 经向风
Zonal_wind = V_real_W;          % 纬向风
Vertical_wind = V_real_V;       %垂直风
figure('Name','钠原子风速反演图');
plot(V_real_length,height,'-b','linewidth',2);
xlim([0 120]);
ylim([80 105]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Wind Velocity 20230709 23:05');
xlabel('Wind Velocity (m/s)');
ylabel('Altitude (km)');

% 画风场图,53.49N,122.34E
figure('Name','钠原子风场反演图');
X = [53.49];
X = repmat(X,221,1);
Y = [122.34];
Y = repmat(Y,221,1);
quiver3(X, Y, height', V_real_direction(:,1), V_real_direction(:,2), V_real_direction(:,3), 'AutoScale','off','MaxHeadSize', 1 , 'LineWidth', 2);
zlim([85 100]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Wind Direction 20230709 23:05');
xlabel('Latitude (N)');
ylabel('Longitude (E)');
zlabel('Altitude (km)');








