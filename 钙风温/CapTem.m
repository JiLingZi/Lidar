%% 钙离子风温反演

%% 低空散射截面-风温曲面拟合

% 通用常数
c = 2.99792458e8;                           %光速(m/s)
k_B = 1.3806505e-23;                        %玻尔兹曼常数(J/K)
e = 1.60217662e-19;                         %电子电量(C)
m_e = 9.10938215e-31;                       %电子质量(kg)
epsilon_0 = 8.854187817e-12;                %真空电容率(F/m)

% 激光常数
lambda_L = 393.477469e-9;                     %激光波长(m)
nu_L = c / lambda_L;                        %激光中心频率(Hz)
sigma_FWHM = 150e6;                         %激光线宽
sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));  %高斯线型的RMS宽度
% OffsFreq = 650;                             %左右频移量(MHz)

% 偏移量
lambda0 = 393.477469e-9;
lambdaR = 393.477048e-9;
lambdaL = 393.477845e-9;
DlambdaR = lambdaR-lambda0;
DlambdaL = lambdaL-lambda0;
DnuR = round((1e-6)*(-1*c*DlambdaR)/(lambda0^2));
DnuL = round((1e-6)*(-1*c*DlambdaL)/(lambda0^2));

% 钙离子常数
M = 6.665e-26;                              %钙离子绝对质量(kg)
fD2 = 0.69;                                 %钙离子振子强度

% 低空风温数据
% T = [80 100 120 140 160 180 200 220 240 260 280 300 320];
T = 80:10:350;
% V = [-120 -100 -80 -60 -40 -20 0 20 40 60 80 100 120];
V = -120:10:120;

% 有效后向散射截面循环矩阵
Data_eff_0 = zeros(size(T,2),size(V,2));
Data_eff_R = zeros(size(T,2),size(V,2));
Data_eff_L = zeros(size(T,2),size(V,2));

% 进入风速循环
for j = 1:length(V)
    
    nu_0_1 = nu_L * ((c+V(j))/c);
    f1 = 1;
    nu = nu_L-3e9:1e6:nu_L+3e9;
    g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));
    
    % 进入温度循环
    for i = 1:length(T)
        sigma_D_1 = nu_0_1 * sqrt(k_B*T(i) / (M*c^2));
        sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_abs_1 = sigma_0_1 * exp(-(nu-nu_0_1).^2 / (2*sigma_D_1^2));
        sigma_abs_1_fi = sigma_abs_1*f1;
        sigma_abs = sigma_abs_1_fi;
        sigma_eff_fx = conv(sigma_abs,g_L);
        sigma_eff_fx = sigma_eff_fx * 1e6;
        Data_eff_0(i,j) = sigma_eff_fx(6001)/(4*pi);
        Data_eff_R(i,j) = sigma_eff_fx(6001+round(1*DnuR))/(4*pi);
        Data_eff_L(i,j) = sigma_eff_fx(6001+round(1*DnuL))/(4*pi);
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
f_TL = fit_T;

% % 画出温度拟合曲面
% figure('Name','温度拟合曲面');
% plot(f_TL);
% set(gca,'FontSize',12,'FontName','Times New Roman');
% title('Fitting Surface of Temperature');
% xlabel('R_{V}');
% ylabel('R_{T}');
% zlabel('Temperature (K)');

% 对风速进行二维曲面拟合
ft = fittype( 'poly55' );
[fit_V, gof_V] = fit( [R_V(:), R_T(:)], V(:), ft );
coe_fit_v=coeffvalues(fit_V);
f_VL = fit_V;

% % 画出风速拟合曲面
% figure('Name','风速拟合曲面');
% plot(f_VL);
% set(gca,'FontSize',12,'FontName','Times New Roman');
% title('Fitting Surface of Wind Velocity');
% xlabel('R_{V}');
% ylabel('R_{T}');
% zlabel('Wind Velocity (m/s)');


%% 高空散射截面-风温曲面拟合

% 通用常数
c = 2.99792458e8;                           %光速(m/s)
k_B = 1.3806505e-23;                        %玻尔兹曼常数(J/K)
e = 1.60217662e-19;                         %电子电量(C)
m_e = 9.10938215e-31;                       %电子质量(kg)
epsilon_0 = 8.854187817e-12;                %真空电容率(F/m)

% 激光常数
lambda_L = 393.477469e-9;                     %激光波长(m)
nu_L = c / lambda_L;                        %激光中心频率(Hz)
sigma_FWHM = 150e6;                         %激光线宽
sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));  %高斯线型的RMS宽度
OffsFreq = 650;                             %左右频移量(MHz)

% 钙离子常数
M = 6.665e-26;                              %钙离子绝对质量(kg)
fD2 = 0.69;                                 %钙离子振子强度

% 低空风温数据
% T = [80 100 120 140 160 180 200 220 240 260 280 300 320]+600;
T = 600:20:1100;
% V = [-120 -100 -80 -60 -40 -20 0 20 40 60 80 100 120];
V = -120:20:120;

% 有效后向散射截面循环矩阵
Data_eff_0 = zeros(size(T,2),size(V,2));
Data_eff_R = zeros(size(T,2),size(V,2));
Data_eff_L = zeros(size(T,2),size(V,2));

% 进入风速循环
for j = 1:length(V)
    
    nu_0_1 = nu_L * ((c+V(j))/c);
    f1 = 1;
    nu = nu_L-3e9:1e6:nu_L+3e9;
    g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));
    
    % 进入温度循环
    for i = 1:length(T)
        sigma_D_1 = nu_0_1 * sqrt(k_B*T(i) / (M*c^2));
        sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_abs_1 = sigma_0_1 * exp(-(nu-nu_0_1).^2 / (2*sigma_D_1^2));
        sigma_abs_1_fi = sigma_abs_1*f1;
        sigma_abs = sigma_abs_1_fi;
        sigma_eff_fx = conv(sigma_abs,g_L);
        sigma_eff_fx = sigma_eff_fx * 1e6;
        Data_eff_0(i,j) = sigma_eff_fx(6001)/(4*pi);
        Data_eff_R(i,j) = sigma_eff_fx(6001+DnuR)/(4*pi);
        Data_eff_L(i,j) = sigma_eff_fx(6001+DnuL)/(4*pi);
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
f_TH = fit_T;

% 画出温度拟合曲面
figure('Name','温度拟合曲面');
plot(f_TH);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Fitting Surface of Temperature');
xlabel('R_{V}');
ylabel('R_{T}');
zlabel('Temperature (K)');

% 对风速进行二维曲面拟合
ft = fittype( 'poly55' );
[fit_V, gof_V] = fit( [R_V(:), R_T(:)], V(:), ft );
coe_fit_v=coeffvalues(fit_V);
f_VH = fit_V;

% % 画出风速拟合曲面
% figure('Name','风速拟合曲面');
% plot(f_VH);
% set(gca,'FontSize',12,'FontName','Times New Roman');
% title('Fitting Surface of Wind Velocity');
% xlabel('R_{V}');
% ylabel('R_{T}');
% zlabel('Wind Velocity (m/s)');


%% 风温反演

% 读取钙离子风温雷达采集数据矩阵
Cap_table = readtable('H:\ZWDATA\FregionCap\20231223_5_photon_add_0028.txt');
Cap_data = table2array(Cap_table);
height_num = Cap_data(:,1);
Cap_photon = Cap_data(:,2:4);

% 读取钠风温雷达采集数据矩阵 23:00 to 00:00
folder = 'H:\ZWDATA\FregionCap\20231223_5\';
files = dir(fullfile(folder, '*.dat'));
num_files = length(files); 
Cap_photon = zeros(8192,3);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    Na_table = readtable(filename);
    Na_data = table2array(Na_table);
    Cap_photon = Cap_photon+Na_data(:,2:4);
end
height_num = Na_data(:,1);

% 以1km为空间分辨率进行光子数行合并
SumiNum = 9;
Cap_i_read = zeros(floor((size(Cap_photon,1)/SumiNum)),size(Cap_photon,2));
Alt = zeros(floor((size(Cap_photon,1)/SumiNum)),1);
for i = 1:floor((size(Cap_photon,1)/SumiNum))
    Cap_i_read(i,:) = sum(Cap_photon(1+SumiNum*(i-1):SumiNum*i,:),1);
    Alt(i,:) = sum(height_num(1+SumiNum*(i-1):SumiNum*i,:),1)/SumiNum;
end
Noise = mean(Cap_i_read(118:136,:),1);  % 310-330
Cap_phSum = Cap_i_read-Noise;

% 对不同方向不同频率进行拆分,30KM处归一化
V_F_0 = Cap_phSum(:,1)./Cap_phSum(28,1);
V_F_R = Cap_phSum(:,2)./Cap_phSum(28,2);
V_F_L = Cap_phSum(:,3)./Cap_phSum(28,3);
NorBack = (Cap_phSum(4,1)+Cap_phSum(4,2)+Cap_phSum(4,3))./3;

% RTRV
R_V_real_V = (V_F_R - V_F_L) ./ (V_F_0);  % 铅垂RV
R_T_real_V = (V_F_R + V_F_L) ./ (2 * V_F_0);% 铅垂RT

% R_V_real_V = (V_F_L - V_F_0) ./ (V_F_R);  % 铅垂RV
% R_T_real_V = (V_F_L + V_F_0) ./ (2 * V_F_R);% 铅垂RT

% 计算低空温度
T_realL = f_TL(R_V_real_V,R_T_real_V);
% 计算高空温度
T_realH = f_TH(R_V_real_V,R_T_real_V);

% 计算风速
V_realL = f_VL(R_V_real_V,R_T_real_V);
V_realH = f_VH(R_V_real_V,R_T_real_V);


%% 计算误差
N0 = V_F_0.*NorBack;
NR = V_F_R.*NorBack;
NL = V_F_L.*NorBack;
RT_real = (NR+NL)./(2*N0);
RV_real = (NR-NL)./N0;

% 定义 RT RV 变量
syms RT RV;

dvalueTRTL = zeros(size(RT_real,1),1);
dvalueTRVL = zeros(size(RT_real,1),1);
dvalueTRTH = zeros(size(RT_real,1),1);
dvalueTRVH = zeros(size(RT_real,1),1);
dvalueVRTL = zeros(size(RT_real,1),1);
dvalueVRVL = zeros(size(RT_real,1),1);
dvalueVRTH = zeros(size(RT_real,1),1);
dvalueVRVH = zeros(size(RT_real,1),1);

for i = 1:size(RT_real,1)
    % 定义 fTL 函数
    tL00 =        29.9;
    tL10 =   4.622e-05;
    tL01 =         393;
    tL20 =      -428.8;
    tL11 =   -0.001585;
    tL02 =      -775.7;
    tL30 =    0.001282;
    tL21 =        1567;
    tL12 =    0.004152;
    tL03 =        2682;
    tL40 =      -269.5;
    tL31 =   -0.004658;
    tL22 =       -2159;
    tL13 =   -0.003733;
    tL04 =       -3849;
    tL50 =  -1.957e-05;
    tL41 =       541.5;
    tL32 =    0.005158;
    tL23 =       260.6;
    tL14 =   -0.001483;
    tL05 =        3056;
    fTL(RT,RV) = tL00+tL10*RT+tL01*RV+tL20*RT^2+tL11*RT*RV+tL02*RV^2+tL30*RT^3 ...
             +tL21*RT^2*RV+tL12*RT*RV^2+tL03*RV^3+tL40*RT^4+tL31*RT^3*RV ...
             +tL22*RT^2*RV^2+tL13*RT*RV^3+tL04*RV^4+tL50*RT^5+tL41*RT^4*RV ... 
             +tL32*RT^3*RV^2+tL23*RT^2*RV^3+tL14*RT*RV^4+tL05*RV^5;

     % 对 fTL 的 RT 和 RV 分别求偏导,并根据实测所得RT RV求偏导常数
    dT_dRTL = diff(fTL(RT,RV),RT);
    dT_dRVL = diff(fTL(RT,RV),RV);
    dvalueTRTL(i,:) = subs(dT_dRTL,[RT RV],[RV_real(i,:) RT_real(i,:)]);
    dvalueTRTL(i,:) = double(dvalueTRTL(i,:));
    dvalueTRVL(i,:) = subs(dT_dRVL,[RT RV],[RV_real(i,:) RT_real(i,:)]);
    dvalueTRVL(i,:) = double(dvalueTRVL(i,:));

    % 定义 fTH 函数
    tH00 =  -8.791e+05;
    tH10 =      -2.052;
    tH01 =   5.675e+06;
    tH20 =   2.466e+05;
    tH11 =       10.72;
    tH02 =  -1.468e+07;
    tH30 =      0.3725;
    tH21 =  -9.372e+05;
    tH12 =      -21.06;
    tH03 =   1.902e+07;
    tH40 =       -5687;
    tH31 =     -0.9523;
    tH22 =   1.191e+06;
    tH13 =       18.41;
    tH04 =  -1.236e+07;
    tH50 =   -0.003054;
    tH41 =        7596;
    tH32 =      0.6126;
    tH23 =  -5.075e+05;
    tH14 =      -6.058;
    tH05 =   3.223e+06;
    fTH(RT,RV) = tH00+tH10*RT+tH01*RV+tH20*RT^2+tH11*RT*RV+tH02*RV^2+tH30*RT^3 ...
             +tH21*RT^2*RV+tH12*RT*RV^2+tH03*RV^3+tH40*RT^4+tH31*RT^3*RV ...
             +tH22*RT^2*RV^2+tH13*RT*RV^3+tH04*RV^4+tH50*RT^5+tH41*RT^4*RV ... 
             +tH32*RT^3*RV^2+tH23*RT^2*RV^3+tH14*RT*RV^4+tH05*RV^5;

    % 对 fTH 的 RT 和 RV 分别求偏导,并根据实测所得RT RV求偏导常数
    dT_dRTH = diff(fTH(RT,RV),RT);
    dT_dRVH = diff(fTH(RT,RV),RV);
    dvalueTRTH(i,:) = subs(dT_dRTH,[RT RV],[RV_real(i,:) RT_real(i,:)]);
    dvalueTRTH(i,:) = double(dvalueTRTH(i,:));
    dvalueTRVH(i,:) = subs(dT_dRVH,[RT RV],[RV_real(i,:) RT_real(i,:)]);
    dvalueTRVH(i,:) = double(dvalueTRVH(i,:));

    % 定义 fVL 函数
    vL00 =  -4.435e-05;
    vL10 =       361.7;
    vL01 =   0.0006769;
    vL20 =  -0.0001372;
    vL11 =       -1311;
    vL02 =   -0.003888;
    vL30 =       36.46;
    vL21 =   0.0006721;
    vL12 =        3187;
    vL03 =     0.01076;
    vL40 =  -2.887e-05;
    vL31 =      -267.3;
    vL22 =   -0.000916;
    vL13 =       -3379;
    vL04 =    -0.01448;
    vL50 =       15.53;
    vL41 =   4.997e-05;
    vL32 =         106;
    vL23 =    0.000276;
    vL14 =        1662;
    vL05 =    0.007595;
    fVL(RT,RV) = vL00+vL10*RT+vL01*RV+vL20*RT^2+vL11*RT*RV+vL02*RV^2+vL30*RT^3 ...
             +vL21*RT^2*RV+vL12*RT*RV^2+vL03*RV^3+vL40*RT^4+vL31*RT^3*RV ...
             +vL22*RT^2*RV^2+vL13*RT*RV^3+vL04*RV^4+vL50*RT^5+vL41*RT^4*RV ... 
             +vL32*RT^3*RV^2+vL23*RT^2*RV^3+vL14*RT*RV^4+vL05*RV^5;

    % 对 fVL 的 RT 和 RV 分别求偏导
    dVL_dRT = diff(fVL(RT,RV),RT);
    dVL_dRV = diff(fVL(RT,RV),RV);
    dvalueVRTL(i,:) = subs(dVL_dRT,[RT RV],[RV_real(i,:) RT_real(i,:)]);
    dvalueVRTL(i,:) = double(dvalueVRTL(i,:));
    dvalueVRVL(i,:) = subs(dVL_dRV,[RT RV],[RV_real(i,:) RT_real(i,:)]);
    dvalueVRVL(i,:) = double(dvalueVRVL(i,:));

    % 定义 fVH 函数
    vH00 =       -1.35;
    vH10 =   1.258e+05;
    vH01 =       8.304;
    vH20 =    -0.08086;
    vH11 =  -6.372e+05;
    vH02 =      -20.42;
    vH30 =  -1.646e+04;
    vH21 =      0.2889;
    vH12 =   1.216e+06;
    vH03 =        25.1;
    vH40 =  -0.0004014;
    vH31 =    4.23e+04;
    vH22 =     -0.3436;
    vH13 =  -1.036e+06;
    vH04 =      -15.42;
    vH50 =       146.7;
    vH41 =   0.0004509;
    vH32 =  -2.757e+04;
    vH23 =      0.1361;
    vH14 =   3.329e+05;
    vH05 =       3.786;
    fVH(RT,RV) = vH00+vH10*RT+vH01*RV+vH20*RT^2+vH11*RT*RV+vH02*RV^2+vH30*RT^3 ...
             +vH21*RT^2*RV+vH12*RT*RV^2+vH03*RV^3+vH40*RT^4+vH31*RT^3*RV ...
             +vH22*RT^2*RV^2+vH13*RT*RV^3+vH04*RV^4+vH50*RT^5+vH41*RT^4*RV ... 
             +vH32*RT^3*RV^2+vH23*RT^2*RV^3+vH14*RT*RV^4+vH05*RV^5;

    % 对 fVH 的 RT 和 RV 分别求偏导
    dVH_dRT = diff(fVH(RT,RV),RT);
    dVH_dRV = diff(fVH(RT,RV),RV);
    dvalueVRTH(i,:) = subs(dVH_dRT,[RT RV],[RV_real(i,:) RT_real(i,:)]);
    dvalueVRTH(i,:) = double(dvalueVRTH(i,:));
    dvalueVRVH(i,:) = subs(dVH_dRV,[RT RV],[RV_real(i,:) RT_real(i,:)]);
    dvalueVRVH(i,:) = double(dvalueVRVH(i,:));
end
% 定义并计算Pre_TL
Pre_TL = sqrt((NR.*(0.5.*dvalueTRTL+dvalueTRVL).^2+NL.*(0.5.*dvalueTRTL-dvalueTRVL).^2+N0.*(RT_real.*dvalueTRTL+RV_real.*dvalueTRVL).^2)./(N0.^2));
Pre_TL = real(double(Pre_TL));

% 定义并计算Pre_TH
Pre_TH = sqrt((NR.*(0.5.*dvalueTRTH+dvalueTRVH).^2+NL.*(0.5.*dvalueTRTH-dvalueTRVH).^2+N0.*(RT_real.*dvalueTRTH+RV_real.*dvalueTRVH).^2)./(N0.^2));
Pre_TH = real(double(Pre_TH));

%% 画图区，不计算

% 画低空温度图
yyyy = 2023;
dd = 358;
hh = 15;
ModelAlt = 1:1:300;
[T_model,Den_model] = atmosnrlmsise00(ModelAlt*1e3,40.5,116.0,yyyy,dd,hh);
T_model = T_model(:,2);
figure('Name','钙离子低空温度反演图');
plot(T_realL,Alt,'-b','linewidth',2);
hold on
plot(T_model,ModelAlt,'-k','linewidth',2);
hold on;
% errorbar(T_realL, Alt, Pre_TL, 'horizontal', 'b.',  'LineWidth', 2);
legend('T_{Wang}', 'T_{Msise-00}', 'Location', 'southeast');
set(legend, 'EdgeColor', 'w');
xlim([160 300]);
ylim([85 105]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Cap Tem [85 105] km');
xlabel('Temperature (K)');
ylabel('Altitude (km)');

% 读取Before温度文件
Before_table = readtable('H:\ZWDATA\FregionCap\20220605_2000_2700_temperature_wind_H.txt', 'ReadVariableNames', true, 'HeaderLines', 8);
Before_data = table2array(Before_table);
height_before = Before_data(:,1);
Before_Tem = Before_data(:,4);
Before_Err = Before_data(:,5);
figure('Name','钙离子低空温度反演图');
plot(Before_Tem,height_before,'-r','linewidth',2);
hold on
plot(T_realH,Alt,'-b','linewidth',2);
hold on
plot(T_model,ModelAlt,'-k','linewidth',2);
hold on;
errorbar(T_realH, Alt, Pre_TH, 'horizontal', 'b.',  'LineWidth', 2);
hold on;
errorbar(Before_Tem, height_before, Before_Err, 'horizontal', 'r.',  'LineWidth', 2);
legend('T_{Wu}', 'T_{Wang}', 'T_{Msise-00}', 'Location', 'southeast');
set(legend, 'EdgeColor', 'w');
xlim([500 1200]);
ylim([220 290]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Cap Tem [220 290] km');
xlabel('Temperature (K)');
ylabel('Altitude (km)');





