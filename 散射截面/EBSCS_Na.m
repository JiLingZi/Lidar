%{
项目名称：钠原子D2线的有效后向散射截面计算程序
作者名称：王可鑫
作者单位：江苏科技大学理学院
作者邮箱：wangkexin1998@qq.com
导师名称：王泽龙
版本编号：MATLAB-3
创建时间：2022/09/15/14:31:51
修改时间：2022/09/18/10:58:07
%}

%% 常数赋值区
%通用常数
c = 2.99792458e8;%光速(m/s)
V = 150;
k_B = 1.3806505e-23;%玻尔兹曼常数(J/K)
e = 1.60217662e-19;%电子电量(C)
m_e = 9.10938215e-31;%电子质量(kg)
epsilon_0 = 8.854187817e-12;%真空电容率(F/m)

%激光常数  589.1583e-9
lambda_L = 589.1576e-9;%激光波长(m),【铁、钠流星尾迹的激光雷达观测研究及钠层全天时观测技术,P12】
nu_L = (c) / lambda_L;%激光中心频率(Hz)
%专用常数
T = 200;%温度为200K条件下【铁、钠流星尾迹的激光雷达观测研究及钠层全天时观测技术，P14】
M = 3.82e-26;%Na原子绝对质量(kg)【根据相对原子质量计算得出】
fD2 = 0.6408;%Na原子D2线的振子强度【高空钠层、钾层同时探测的激光雷达】
IsoAbdn = 1;%大气Na没有同位素，丰度取1
sigma_FWHM = 80e6;%激光测量线宽1.5GHz；半高全宽【高空钠层、钾层同时探测的激光雷达】
sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));%高斯线型的RMS宽度

% 谱线参数赋值区
%Na的6条跃迁谱线的中心频率偏差值【正号，因为激光线型不动，动的是光谱】 + 0.6302e9
nu01 = (nu_L - 0.7328e9) * ((c+V)/c);%Na第1条跃迁谱线的中心频率【铁、钠流星尾迹的激光雷达观测研究及钠层全天时观测技术,P12】
nu02 = (nu_L - 0.6962e9) * ((c+V)/c);%Na-2
nu03 = (nu_L - 0.6302e9) * ((c+V)/c);%Na-3
nu04 = (nu_L + 1.0333e9) * ((c+V)/c);%Na-4
nu05 = (nu_L + 1.0552e9) * ((c+V)/c);%Na-5
nu06 = (nu_L + 1.0919e9) * ((c+V)/c);%Na-6
%Na的超精细跃迁强度因子【铁、钠流星尾迹的激光雷达观测研究及钠层全天时观测技术,P12】
f1 = 1/32;
f2 = 5/32;
f3 = 14/32;
f4 = 2/32;
f5 = 5/32;
f6 = 5/32;

% 归一化线型，1.11

syms nu;%定义变量频率
nu = nu_L-2.5e9:1e7:nu_L+2.5e9;
%归一化激光线型
g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));%公式1.11

% 多普勒增宽的RMS宽度，1.14
%Na的6条线的多普勒增宽的RMS宽度
sigma_D_1 = nu01 * sqrt(k_B*T / (M*c^2));
sigma_D_2 = nu02 * sqrt(k_B*T / (M*c^2));
sigma_D_3 = nu03 * sqrt(k_B*T / (M*c^2));
sigma_D_4 = nu04 * sqrt(k_B*T / (M*c^2));
sigma_D_5 = nu05 * sqrt(k_B*T / (M*c^2));
sigma_D_6 = nu06 * sqrt(k_B*T / (M*c^2));

% 峰值处的吸收截面，1.15
%Na的6条线的峰值处的吸收截面
sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
sigma_0_2 = (1 / (sqrt(2*pi)*sigma_D_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
sigma_0_3 = (1 / (sqrt(2*pi)*sigma_D_3)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
sigma_0_4 = (1 / (sqrt(2*pi)*sigma_D_4)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
sigma_0_5 = (1 / (sqrt(2*pi)*sigma_D_5)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
sigma_0_6 = (1 / (sqrt(2*pi)*sigma_D_6)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;

% 单个吸收截面，1.13
%Na的6条线的吸收截面
sigma_abs_1 = sigma_0_1 * exp(-(nu-nu01).^2 / (2*sigma_D_1^2));
sigma_abs_2 = sigma_0_2 * exp(-(nu-nu02).^2 / (2*sigma_D_2^2));
sigma_abs_3 = sigma_0_3 * exp(-(nu-nu03).^2 / (2*sigma_D_3^2));
sigma_abs_4 = sigma_0_4 * exp(-(nu-nu04).^2 / (2*sigma_D_4^2));
sigma_abs_5 = sigma_0_5 * exp(-(nu-nu05).^2 / (2*sigma_D_5^2));
sigma_abs_6 = sigma_0_6 * exp(-(nu-nu06).^2 / (2*sigma_D_6^2));

% Na的6条吸收截面乘以超精细跃迁强度因子绘图区
sigma_abs_1_fi = sigma_abs_1*f1;
sigma_abs_2_fi = sigma_abs_2*f2;
sigma_abs_3_fi = sigma_abs_3*f3;
sigma_abs_4_fi = sigma_abs_4*f4;
sigma_abs_5_fi = sigma_abs_5*f5;
sigma_abs_6_fi = sigma_abs_6*f6;

% Na的吸收截面叠加绘图区
sigma_abs = IsoAbdn .* (sigma_abs_1_fi + sigma_abs_2_fi + sigma_abs_3_fi + sigma_abs_4_fi + sigma_abs_5_fi + sigma_abs_6_fi);

figure('Name','Na Scatter')
plot(nu-nu_L,sigma_abs,'-k','linewidth',2)
hold on
plot([-6.5e8,-6.5e8],[0,1e-15],'--m','linewidth',1.5);
hold on
plot([-6.5e8-585e6,-6.5e8-585e6],[0,1e-15],'--b','linewidth',1.5);
hold on
plot([-6.5e8+585e6,-6.5e8+585e6],[0,1e-15],'--r','linewidth',1.5);
ylim([0,1e-15]);
title('Na Scattering Cross Section');
xlabel('Offset Frequency (GHz)');
ylabel('Cross Section (m^2/sr)');
set(gca,'FontName','Times New Roman','FontSize',15);

%% 有效后向散射截面，画图时注释掉收起来,计算时记得关闭绘图精度

sigma_eff_fx = sigma_abs.*g_L;
sigma_eff = int(sigma_eff_fx,nu,-inf,inf);
sigma_K = sigma_eff / (4*pi);
sprintf('Na原子的D2线的有效后向散射截面为\n%e m^2/sr',sigma_K)
