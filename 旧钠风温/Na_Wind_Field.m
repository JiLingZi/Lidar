%-------------------------------------------------------------------------
%======================中高层大气钠风场反演程序============================
%-------------------------------------------------------------------------
% 项目名称：漠河钠原子风场反演程序
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

%% 常数赋值区

% 通用常数
c = 2.99792458e8;                           %光速(m/s)
k_B = 1.3806505e-23;                        %玻尔兹曼常数(J/K)
e = 1.60217662e-19;                         %电子电量(C)
m_e = 9.10938215e-31;                       %电子质量(kg)
epsilon_0 = 8.854187817e-12;                %真空电容率(F/m)

% 激光常数
lambda_L = 589.1583e-9;                     %激光波长(m)
nu_L = c / lambda_L;                        %激光中心频率(Hz)
sigma_FWHM = 80e6;                          %激光线宽
sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));  %高斯线型的RMS宽度

% 钠原子常数
M = 3.82e-26;                               %Na原子绝对质量(kg)
fD2 = 0.6408;                               %Na原子D2线的振子强度
%% 风温数据区

% 温度数据
T = [80 100 120 140 160 180 200 220 240 260 280 300 320];

% 风速数据
V = [-120 -100 -80 -60 -40 -20 0 20 40 60 80 100 120];

%% 有效后向散射截面矩阵循环

Data_eff_0 = zeros(size(T,2),size(V,2));
Data_eff_R = zeros(size(T,2),size(V,2));
Data_eff_L = zeros(size(T,2),size(V,2));

% 进入风速循环
for j = 1:length(V)    
    %% 谱线参数赋值区
    
    % Na的6条跃迁谱线的中心频率偏差值(Hz)
    nu01 = (nu_L - 0.7328e9 + 0.6302e9) * ((c+V(j))/c);
    nu02 = (nu_L - 0.6962e9 + 0.6302e9) * ((c+V(j))/c);
    nu03 = (nu_L - 0.6302e9 + 0.6302e9) * ((c+V(j))/c);
    nu04 = (nu_L + 1.0333e9 + 0.6302e9) * ((c+V(j))/c);
    nu05 = (nu_L + 1.0552e9 + 0.6302e9) * ((c+V(j))/c);
    nu06 = (nu_L + 1.0919e9 + 0.6302e9) * ((c+V(j))/c);
    
    % Na的超精细跃迁强度因子
    f1 = 1/32;
    f2 = 5/32;
    f3 = 14/32;
    f4 = 2/32;
    f5 = 5/32;
    f6 = 5/32;
    
    %% 归一化线型
    nu = nu_L-3e9:1e6:nu_L+3e9;
    g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));
    
    % 温度循环
    for i = 1:length(T)
        %% 多普勒增宽的RMS宽度
        sigma_D_1 = nu01 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_2 = nu02 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_3 = nu03 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_4 = nu04 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_5 = nu05 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_6 = nu06 * sqrt(k_B*T(i) / (M*c^2));
       
        %% 峰值处的吸收截面
        sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_2 = (1 / (sqrt(2*pi)*sigma_D_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_3 = (1 / (sqrt(2*pi)*sigma_D_3)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_4 = (1 / (sqrt(2*pi)*sigma_D_4)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_5 = (1 / (sqrt(2*pi)*sigma_D_5)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_6 = (1 / (sqrt(2*pi)*sigma_D_6)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        
        %% 六条跃迁线的单个吸收截面
        sigma_abs_1 = sigma_0_1 * exp(-(nu-nu01).^2 / (2*sigma_D_1^2));
        sigma_abs_2 = sigma_0_2 * exp(-(nu-nu02).^2 / (2*sigma_D_2^2));
        sigma_abs_3 = sigma_0_3 * exp(-(nu-nu03).^2 / (2*sigma_D_3^2));
        sigma_abs_4 = sigma_0_4 * exp(-(nu-nu04).^2 / (2*sigma_D_4^2));
        sigma_abs_5 = sigma_0_5 * exp(-(nu-nu05).^2 / (2*sigma_D_5^2));
        sigma_abs_6 = sigma_0_6 * exp(-(nu-nu06).^2 / (2*sigma_D_6^2));
       
        %% 六条吸收截面乘以超精细跃迁强度因子
        sigma_abs_1_fi = sigma_abs_1*f1;
        sigma_abs_2_fi = sigma_abs_2*f2;
        sigma_abs_3_fi = sigma_abs_3*f3;
        sigma_abs_4_fi = sigma_abs_4*f4;
        sigma_abs_5_fi = sigma_abs_5*f5;
        sigma_abs_6_fi = sigma_abs_6*f6;
       
        %% 吸收截面叠加
        sigma_abs = sigma_abs_1_fi + sigma_abs_2_fi + sigma_abs_3_fi + sigma_abs_4_fi + sigma_abs_5_fi + sigma_abs_6_fi;
        
        %% 有效后向散射截面
        sigma_eff_fx = conv(sigma_abs,g_L);
        sigma_eff_fx = sigma_eff_fx * 1e6;
        Data_eff_0(i,j) = sigma_eff_fx(6001)/(4*pi);
        Data_eff_R(i,j) = sigma_eff_fx(6001+585)/(4*pi);
        Data_eff_L(i,j) = sigma_eff_fx(6001-585)/(4*pi);
    end
end

% 定义拟合风速和温度数据
T = T';
T = repmat(T,1,size(V,2));
V = repmat(V,size(T,1),1);

% 风温比率矩阵
R_T = (Data_eff_R + Data_eff_L) ./ (2 * Data_eff_0);
R_V = (Data_eff_R - Data_eff_L) ./ (Data_eff_0);

% 对温度进行二维曲面拟合
ft = fittype( 'poly55' );
[fit_T, gof_T] = fit( [R_V(:), R_T(:)], T(:), ft );
coe_fit_t=coeffvalues(fit_T);
f_T = fit_T;

% 对风速进行二维曲面拟合
ft = fittype( 'poly55' );
[fit_V, gof_V] = fit( [R_V(:), R_T(:)], V(:), ft );
coe_fit_v=coeffvalues(fit_V);
f_V = fit_V;

% 读取钠风温雷达采集数据矩阵
folder = 'D:\ZWDATA\Mohe_20230625_0712\Field\';
% folder = 'D:\ZWDATA\Mohe_20230625_0712\Na07112300_0030\';
files = dir(fullfile(folder, '*.dat'));
num_files = length(files);  
V_F_0 = zeros(8192, num_files);
V_F_R = zeros(8192, num_files);
V_F_L = zeros(8192, num_files);
N_F_0 = zeros(8192, num_files);
N_F_R = zeros(8192, num_files);
N_F_L = zeros(8192, num_files);
W_F_0 = zeros(8192, num_files);
W_F_R = zeros(8192, num_files);
W_F_L = zeros(8192, num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    Na_table = readtable(filename);
    Na_data = table2array(Na_table);
    V_F_0(:,j) = Na_data(:,2);
    V_F_R(:,j) = Na_data(:,3);
    V_F_L(:,j) = Na_data(:,4);
    N_F_0(:,j) = Na_data(:,5);
    N_F_R(:,j) = Na_data(:,6);
    N_F_L(:,j) = Na_data(:,7);
    W_F_0(:,j) = Na_data(:,8);
    W_F_R(:,j) = Na_data(:,9);
    W_F_L(:,j) = Na_data(:,10);
end
height_num = Na_data(:,1);
H_FUCK = height_num*(sqrt(3)/2);

% 扣除130 to 150 km 噪声,并将经纬向光子数乘以cos(30)°后以36km瑞利信号归一化*(sqrt(3)/2))
V_F_0 = (V_F_0 - mean(V_F_0(4232:4558,:),1))./mean(V_F_0(977:1140,:),1);
V_F_R = (V_F_R - mean(V_F_R(4232:4558,:),1))./mean(V_F_R(977:1140,:),1);
V_F_L = (V_F_L - mean(V_F_L(4232:4558,:),1))./mean(V_F_L(977:1140,:),1);
N_F_0 = (N_F_0 - mean(N_F_0(4232:4558,:),1))./mean(N_F_0(977:1140,:),1);
N_F_R = (N_F_R - mean(N_F_R(4232:4558,:),1))./mean(N_F_R(977:1140,:),1);
N_F_L = (N_F_L - mean(N_F_L(4232:4558,:),1))./mean(N_F_L(977:1140,:),1);
W_F_0 = (W_F_0 - mean(W_F_0(4232:4558,:),1))./mean(W_F_0(977:1140,:),1);
W_F_R = (W_F_R - mean(W_F_R(4232:4558,:),1))./mean(W_F_R(977:1140,:),1);
W_F_L = (W_F_L - mean(W_F_L(4232:4558,:),1))./mean(W_F_L(977:1140,:),1);

% 每1km合并光子数行
V_F_0_READ = zeros(256, size(V_F_0,2));
V_F_R_READ = zeros(256, size(V_F_0,2));
V_F_L_READ = zeros(256, size(V_F_0,2));
N_F_0_READ = zeros(256, size(V_F_0,2));
N_F_R_READ = zeros(256, size(V_F_0,2));
N_F_L_READ = zeros(256, size(V_F_0,2));
W_F_0_READ = zeros(256, size(V_F_0,2));
W_F_R_READ = zeros(256, size(V_F_0,2));
W_F_L_READ = zeros(256, size(V_F_0,2));
for i = 1:256
    V_F_0_READ(i,:) = sum(V_F_0(1+32*(i-1):32*i,:),1);
    V_F_R_READ(i,:) = sum(V_F_R(1+32*(i-1):32*i,:),1);
    V_F_L_READ(i,:) = sum(V_F_L(1+32*(i-1):32*i,:),1);
    N_F_0_READ(i,:) = sum(N_F_0(1+32*(i-1):32*i,:),1);
    N_F_R_READ(i,:) = sum(N_F_R(1+32*(i-1):32*i,:),1);
    N_F_L_READ(i,:) = sum(N_F_L(1+32*(i-1):32*i,:),1);
    W_F_0_READ(i,:) = sum(W_F_0(1+32*(i-1):32*i,:),1);
    W_F_R_READ(i,:) = sum(W_F_R(1+32*(i-1):32*i,:),1);
    W_F_L_READ(i,:) = sum(W_F_L(1+32*(i-1):32*i,:),1);    
end
V_F_0 = V_F_0_READ;
V_F_R = V_F_R_READ;
V_F_L = V_F_L_READ;
N_F_0 = N_F_0_READ;
N_F_R = N_F_R_READ;
N_F_L = N_F_L_READ;
W_F_0 = W_F_0_READ;
W_F_R = W_F_R_READ;
W_F_L = W_F_L_READ;
height = 1:256;
height = height * 0.9830390625;
height = height(1:221);

% 计算RV垂直方向、北向、东向
R_V_real_V = (V_F_R - V_F_L) ./ (V_F_0);        % 铅垂RV
R_V_real_N = (N_F_R - N_F_L) ./ (N_F_0);        % 北向(经向)RV
R_V_real_W = (W_F_R - W_F_L) ./ (W_F_0);        % 东向(纬向)RV

% 计算RT
R_T_real_V = (V_F_R + V_F_L) ./ (2 * V_F_0);    % 铅垂RT
R_T_real_N = (N_F_R + N_F_L) ./ (2 * N_F_0);    % 北向(经向)RT
R_T_real_W = (W_F_R + W_F_L) ./ (2 * W_F_0);    % 东向(纬向)RT

% 计算垂直风速
V_real_V = f_V(R_V_real_V,R_T_real_V);
V_real_V = V_real_V(1:221,:);

% 画垂直风时变图
figure('Name','钠原子垂直风');
X = 1:num_files;
Y = height';
window = fspecial('average', [2,5]);
V_V_Window = imfilter(V_real_V, window, 'symmetric', 'same');
[Map, Line]=contourf(X,Y(87:102,:),V_V_Window(87:102,:),7);
set(Line,'LineColor','k');
clabel(Map, Line);
% xlim([120 320]);
ylim([85 100]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Vertical Wind 20230709/10');
xlabel('Local Time');
ylabel('Altitude (km)');
colorbar;
ylabel(colorbar,'Wind Velocity (m/s)');
colormap(jet);

% 计算北向（经向风）风速,并乘以sin（30°）
V_real_N = zeros(221,size(V_F_0,2));
for i = 1:221
    V_real_N(i,:) = f_V(R_V_real_N(round(i*(2/sqrt(3))),:),R_T_real_N(round(i*(2/sqrt(3))),:));
end
V_real_N = V_real_N * 0.5;

% 画经向风时变图
figure('Name','钠原子经向风');
X = 1:num_files;
Y = height';
window = fspecial('average', [1,5]);
V_N_Window = imfilter(V_real_N, window, 'symmetric', 'same');
[Map, Line]=contourf(X,Y(87:102,:),V_N_Window(87:102,:),7);
set(Line,'LineColor','k');
% clabel(Map, Line);
% xlim([120 320]);
ylim([85 100]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Meridional Wind 20230709/10');
xlabel('Local Time');
ylabel('Altitude (km)');
colorbar;
ylabel(colorbar,'Wind Velocity (m/s)');
colormap(jet);

% 计算东向（纬向风）风速,并乘以sin（30°）
V_real_W = zeros(221,size(V_F_0,2));
for i = 1:221
    V_real_W(i,:) = f_V(R_V_real_W(round(i*(2/sqrt(3))),:),R_T_real_W(round(i*(2/sqrt(3))),:));
end
V_real_W = V_real_W * 0.5;

% 画纬向风时变图
figure('Name','钠原子纬向风');
X = 1:num_files;
Y = height';
window = fspecial('average', [1,5]);
V_W_Window = imfilter(V_real_W, window, 'symmetric', 'same');
[Map, Line]=contourf(X,Y(87:102,:),V_W_Window(87:102,:),7);
set(Line,'LineColor','k');
% clabel(Map, Line);
% xlim([120 320]);
ylim([85 100]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Zonal Wind 20230709/10');
xlabel('Local Time');
ylabel('Altitude (km)');
colorbar;
ylabel(colorbar,'Wind Velocity (m/s)');
colormap(jet);

% 给出风速大小和方向0矩阵
% V_real_length = zeros(221,1);
V_real_direction = zeros(221,3);

% 画风场图,53.49N,122.34E[1:size(V_F_0,2)(:,j)]第一行X和quiver3内的X已修改
fig = figure('Name','钠原子风场反演图','position',[500 50 400 580]);
X = zeros(221,1);
Y = zeros(221,1);
V_Vertical = zeros(221,1);
Z = height';
for j = 1:size(V_F_0,2)
    V_real_Array = [V_real_N(:,j),V_real_W(:,j),V_real_V(:,j)];
    for i = 1:221
        V_real_direction(i,:) = V_real_Array(i,:) / norm(V_real_Array(i,:));
    end
    clf;
    % 经向风向量要变成负的才对
    quiver3(X, Y, Z, V_real_direction(:,2), V_real_direction(:,1), V_Vertical, 'AutoScale','off','MaxHeadSize', 0.3, 'LineWidth', 2);
    xlim([-3 3]);
    ylim([-3 2]);
    zlim([85 100]);
    view([-20 15]);
    set(gca,'FontSize',12,'FontName','Times New Roman');
    title('Na Wind Field 20230709/10');
    xlabel('Zonal');
    ylabel('Meridional');
    zlabel('Altitude (km)');
    hold on;
    pause(0.3);
    frames(j) = getframe(fig);% 录制动画
end

figure('Name','播放动画','Position',[500 0 500 680]);
axis off;
% implay(frames,10); 
movie(frames);

% 保存为视频
writerObj = VideoWriter('Na_Wind_Field.mp4');
writerObj.FrameRate = 5; % 设置帧率为每秒10帧
open(writerObj);
writeVideo(writerObj, frames);
close(writerObj);

% 画风场线图,53.49N,122.34E
figure('Name','钠原子风场反演图');
for j = 1:size(V_F_0,2)
    plot(height,0.1.*V_real_N(:,j)+j-1)
    hold on;
end
% 26(43) 28(44) 30(45) 32(46) 34(47) 36(48) 38(49) 40(50) 42(51) 44(52)
xlim([86 100]);
ylim([40 55]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Wind Field 20230709/10');
xlabel('Local Time');
ylabel('53.49N,122.34E');
zlabel('Altitude (km)');

% 画单个风场
plot(height,V_real_N(:,45))
height_pis = height';

