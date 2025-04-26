%{
项目名称：钙原子的有效后向散射截面计算程序
作者名称：王可鑫
作者单位：江苏科技大学理学院
作者邮箱：wangkexin1998@qq.com
导师名称：王泽龙
版本编号：MATLAB-3
创建时间：2022/09/14/14:31:51
修改时间：2022/09/15/10:58:07
%}

%% 常数赋值区
%通用常数
c = 2.99792458e8;%光速(m/s)
k_B = 1.3806505e-23;%玻尔兹曼常数(J/K)
e = 1.60217662e-19;%电子电量(C)
m_e = 9.10938215e-31;%电子质量(kg)
epsilon_0 = 8.854187817e-12;%真空电容率(F/m)

%激光常数
lambda_L = 422.673e-9;%激光波长(m)
nu_L = c / lambda_L;%激光中心频率(Hz)

%专用常数
T = 200;%温度为200K条件下
M = 6.665e-26;%原子绝对质量(kg)
fD1 = 1.75;%振子强度
IsoAbdn = 1;%Fe的同位素丰度
sigma_FWHM = 0.3e9;%激光测量线宽；半高全宽
sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));%高斯线型的RMS宽度

%% 谱线参数赋值区
%跃迁谱线的中心频率偏差值
nu_0_1 =c/422.673e-9;
%超精细跃迁强度因子
f1 = 1;

%% 归一化线型，1.11

syms nu;%定义变量频率

%归一化激光线型绘图精度，绘图精度在计算时必须关闭
%nu = nu_L-2e9:1e6:nu_L+2e9;%绘图时，将nu变量换成等距向量，下公式要点乘


g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));%公式1.11
%{
int_g_L = int(g_L,nu,-inf,inf);%对g_L从-oo到+oo积分验证归一化，得出值为1
sprintf('归一化验证\n%e',int_g_L)
%}
%{
plot(nu-nu_L,g_L,'.-k');%绘制归一化的激光线型
max_g_L = max(g_L);
hold on
plot([0 0], [0 max_g_L],'-k')
set(gca,'FontSize',20);
xlabel('Frequency offset(Hz)');
ylabel('Normalized intensity');
%}


%% 多普勒增宽的RMS宽度，1.14
%Fe的1条线的多普勒增宽的RMS宽度
sigma_D_1 = nu_0_1 * sqrt(k_B*T / (M*c^2));
%sigma_D_2 = nu_0_2 * sqrt(k_B*T / (M*c^2));

%% 峰值处的吸收截面，1.15
%Fe的1条线的峰值处的吸收截面
sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%sigma_0_2 = (1 / (sqrt(2*pi)*sigma_D_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;


%% 绘图精度定义区，不绘图时记得注释掉

%nu = nu_L-2e9:1e5:nu_L+2e9;%绘图时，将nu变量换成等距向量，下公式要点乘

%% 单个吸收截面，1.13
%Fe的1条线的吸收截面
sigma_abs_1 = sigma_0_1 * exp(-(nu-nu_0_1).^2 / (2*sigma_D_1^2));
%sigma_abs_2 = sigma_0_2 * exp(-(nu-nu_0_1).^2 / (2*sigma_D_2^2));

%% 吸收截面叠加绘图区，不绘图时注释掉收起来
%{
plot(nu-nu_L,sigma_abs_1,'.-b');%绘制K39第1条线的吸收截面
hold on

max_sigma_abs_1 = max(sigma_abs_1);
hold on

%记住画直线的起终点是[x1 x2], [y1 y2],也就是先定义两个x，再定义两个y，而不是[x1 y1], [x2 y2]
plot([c / 372.0993,c / 372.0993],[0,max_sigma_abs_1],'.-b')
hold on

set(gca,'FontSize',20);
xlabel('Frequency offset(Hz)');
ylabel('Absorption cross-section(m^2)');
%}

%% 吸收截面乘以超精细跃迁强度因子绘图区
sigma_abs_1_fi = sigma_abs_1*f1;
%sigma_abs_2_fi = sigma_abs_2*f2;

%{
plot(nu-nu_L,sigma_abs_1_fi,'.-b');%绘制K39第1条线的吸收截面
hold on

max_sigma_abs_1_fi = max(sigma_abs_1_fi);
hold on

%记住画直线的起终点是[x1 x2], [y1 y2],也就是先定义两个x，再定义两个y，而不是[x1 y1], [x2 y2]
plot([c / 372.0993,c / 372.0993],[0,max_sigma_abs_1_fi],'.-b')
hold on

set(gca,'FontSize',20);
xlabel('Frequency offset(Hz)');
ylabel('Absorption cross-section(m^2)');
%}


%% K39与K41的吸收截面叠加绘图区
sigma_abs = IsoAbdn * sigma_abs_1_fi;
%sigma_abs_1 = IsoAbdn * (sigma_abs_1_fi);
%sigma_abs_2 = IsoAbdn * (sigma_abs_2_fi);

%fx = max(sigma_abs)/(4*pi)

%{
plot(nu-nu_L,sigma_abs,'.-b');%绘制K39的吸收截面
hold on
plot(nu-nu_L,sigma_abs_41,'.-r');%绘制K39的吸收截面

max_sigma_abs_39 = max(sigma_abs_39);
max_sigma_abs_41 = max(sigma_abs_41);

hold on
%记住画直线的起终点是[x1 x2], [y1 y2],也就是先定义两个x，再定义两个y，而不是[x1 y1], [x2 y2]
plot([-0.096e9,-0.096e9],[0,max_sigma_abs_39],'.-b')
hold on
plot([0.225e9,0.225e9],[0,max_sigma_abs_41],'.-r')

set(gca,'FontSize',20);
xlabel('频率偏差(Hz)');
ylabel('k39与K41的的吸收截面(m^2)');
legend('K39的吸收截面','K41的吸收截面');
%}

%% 钾原子的吸收截面与散射截面叠加绘图区
sigma_abs = sigma_abs;

%fx = max(sigma_abs)/(4*pi)
%{
sigma_eff_fx = conv(sigma_abs,g_L);
sigma_eff_fx = sigma_eff_fx * 1e6;
X = nu_L-2e9:0.5e6:nu_L+2e9;
Y = sigma_eff_fx;
%sigma_eff_fx = sigma_abs.*g_L;
%
plot(nu-nu_L,sigma_abs,'-b','Linewidth',2);%绘制钾原子的吸收截面
hold on
box on
plot((X-nu_L)*2,Y,'.-r','Linewidth',2);%绘制钾原子的散射截面
xlim([-2e9,2e9]);
set(gca,'FontSize',19);
xlabel('Frequency offset(Hz)');
ylabel('Cross-section(m^2)');
legend('Absorption','Effective scattering');
%}

%% 有效后向散射截面，画图时注释掉收起来,计算时记得关闭绘图精度

sigma_eff_fx = sigma_abs.*g_L;
sigma_eff = int(sigma_eff_fx,nu,-inf,inf);
sigma_K = sigma_eff / (4*pi);
sprintf('Ca原子的有效后向散射截面为\n%e m^2/sr',sigma_K)
%{
plot(sigma_FWHM,sigma_K,'-m','Linewidth',2)
set(gca,'FontSize',18);
xlabel('Line width(Hz)');
ylabel('Cross-section(m^2/sr)');

%}