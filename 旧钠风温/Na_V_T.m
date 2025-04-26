%-------------------------------------------------------------------------
%=================中高层大气钠风温激光雷达风温反演程序======================
%-------------------------------------------------------------------------
% 项目名称：漠河钠原子风温反演程序
% 作者名称：王可鑫
% 作者单位：江苏科技大学
% 作者导师：王泽龙
% 作者邮箱：wang.kexin@stu.just.edu.cn
% 创建时间：2023-07-09 17:30(UTC+8)
% 修改时间：2023-07-11 00:30(UTC+8)
% 数据来源：子午工程二期漠河激光雷达台站 
%===============================注意事项==================================
% 本程序仅可反演单个钠风温三频三向文件，如需反演风温随时间变化须加循环
% 本程序作为应急使用，平时不建议使用；
% 本程序跑完一次需要5748秒，请合理安排电脑工作时间，预分配内存
% 平时建议使用‘Na_Wind_Temperature’，以‘Na_Wind_Temperature_EFF’为辅助
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

%% F_0有效后向散射截面矩阵循环

Data_eff_0 = [];

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
    syms nu;
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
        sigma_eff = int(sigma_abs*g_L,nu,-Inf,Inf);
        sigma_eff = sigma_eff / (4*pi);
        sigma_eff = double(sigma_eff);
        Data_eff_0(i,j) = sigma_eff;
    end
end

%% F_+有效后向散射截面矩阵循环

Data_eff_R = [];

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
    syms nu;
    g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-(nu_L+585e6)).^2 ./ (2.*sigma_L.^2));
    
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
        sigma_eff = int(sigma_abs*g_L,nu,-Inf,Inf);
        sigma_eff = sigma_eff / (4*pi);
        sigma_eff = double(sigma_eff);
        Data_eff_R(i,j) = sigma_eff;
    end
end

%% F_-有效后向散射截面矩阵循环

Data_eff_L = [];

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
    syms nu;
    g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-(nu_L-585e6)).^2 ./ (2.*sigma_L.^2));
    
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
        sigma_eff = int(sigma_abs*g_L,nu,-Inf,Inf);
        sigma_eff = sigma_eff / (4*pi);
        sigma_eff = double(sigma_eff);
        Data_eff_L(i,j) = sigma_eff;
    end
end

% 定义拟合风速和温度数据
T = T';
T = repmat(T,1,size(V,2));
V = repmat(V,size(T,1),1);

% 拟合风温比率矩阵
R_T = (Data_eff_R + Data_eff_L) ./ (2 * Data_eff_0);
R_V = (Data_eff_R - Data_eff_L) ./ (Data_eff_0);

% 对温度进行二维曲面拟合
ft = fittype( 'poly55' );
[fit_T, gof_T] = fit( [R_V(:), R_T(:)], T(:), ft );
coe_fit_t=coeffvalues(fit_T);
f_T = fit_T;

% 画出温度拟合曲面
figure('Name','温度拟合曲面');
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
figure('Name','风速拟合曲面');
plot(f_V);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Fitting Surface of Wind Velocity');
xlabel('R_{V}');
ylabel('R_{T}');
zlabel('Wind Velocity (m/s)');

% 读取钠风温雷达采集数据矩阵
Na_table = readtable('G:\Mohe_Data\589Na\Na20230709150515310.dat');
Na_data = table2array(Na_table);
height_num = Na_data(:,1);
Na_photon = Na_data(:,2:10);
Na_photon = Na_photon - mean(Na_photon(4232:4558,:),1);
Na_photon = Na_photon ./ max(Na_photon);

% 以1km为空间分辨率进行光子数行合并
Na_i_read = [];
for i = 1:256
    Na_i_read(i,:) = sum(Na_photon(1+32*(i-1):32*i,:),1);
end
Na_photon = Na_i_read;

% 矩阵个数实际高度修正
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

% 计算RT三向
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

% 画风速图
figure('Name','钠原子风速反演图');
plot(V_real_length,height,'-b','linewidth',2);
xlim([0 120]);
ylim([80 105]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Wind Velocity 20230709 23:05');
xlabel('Wind Velocity (m/s)');
ylabel('Altitude (km)');

% 画风场图
figure('Name','钠原子风场反演图');
X = [53.49];
X = repmat(X,221,1);
Y = [122.34];
Y = repmat(Y,221,1);
quiver3(X, Y, height', V_real_direction(:,1), V_real_direction(:,2), V_real_direction(:,3), 'AutoScale','off', 'MaxHeadSize', 1 , 'LineWidth', 2);
zlim([80 105]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Wind Direction 20230709 23:05');
xlabel('Latitude (N)');
ylabel('Longitude (E)');
zlabel('Altitude (km)');


