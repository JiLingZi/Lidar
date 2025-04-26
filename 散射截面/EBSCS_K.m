%{
项目名称：钾原子D1线的有效后向散射截面计算程序
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
lambda_L = 769.898e-9;%激光波长(m)
nu_L = c / lambda_L;%激光中心频率(Hz)
%专用常数
T = 200;%温度为200K条件下
M_39 = 6.4679e-26;%K39原子绝对质量(kg)
M_41 = 6.7996e-26;%K41原子绝对质量(kg)
fD1 = 0.3327;%钾原子D1线的振子强度
IsoAbdn_39 = 0.933;%K39的同位素丰度
IsoAbdn_41 = 0.067;%K39的同位素丰度
sigma_FWHM =1.5e9;%激光测量线宽1.44GHz；半高全宽
sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));%高斯线型的RMS宽度

%% 谱线参数赋值区
%K39的4条跃迁谱线的中心频率偏差值
nu_0_39_1 = nu_L+0.310e9;%K39第1条跃迁谱线的中心频率
nu_0_39_2 = nu_L+0.254e9;%K39第2条跃迁谱线的中心频率
nu_0_39_3 = nu_L-0.152e9;%K39第3条跃迁谱线的中心频率
nu_0_39_4 = nu_L-0.208e9;%K39第4条跃迁谱线的中心频率
%K41的4条跃迁谱线的中心频率偏差值
nu_0_41_1 = nu_L+0.405e9;%K41第1条跃迁谱线的中心频率
nu_0_41_2 = nu_L+0.375e9;%K41第2条跃迁谱线的中心频率
nu_0_41_3 = nu_L+0.151e9;%K41第3条跃迁谱线的中心频率
nu_0_41_4 = nu_L+0.121e9;%K41第4条跃迁谱线的中心频率
%K的超精细跃迁强度因子
f1 = 5/16;
f2 = 1/16;
f3 = 5/16;
f4 = 5/16;

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
plot([-0.75e9 0.75e9],[3.1e-10 3.1e-10],'.-k')
set(gca,'FontSize',22);
xlabel('Frequency offset(Hz)');
ylabel('Normalized intensity');
%}


%% 多普勒增宽的RMS宽度，1.14
%K39的4条线的多普勒增宽的RMS宽度
sigma_D_39_1 = nu_0_39_1 * sqrt(k_B*T / (M_39*c^2));
sigma_D_39_2 = nu_0_39_2 * sqrt(k_B*T / (M_39*c^2));
sigma_D_39_3 = nu_0_39_3 * sqrt(k_B*T / (M_39*c^2));
sigma_D_39_4 = nu_0_39_4 * sqrt(k_B*T / (M_39*c^2));

%K41的4条线的多普勒增宽的RMS宽度
sigma_D_41_1 = nu_0_41_1 * sqrt(k_B*T / (M_41*c^2));
sigma_D_41_2 = nu_0_41_2 * sqrt(k_B*T / (M_41*c^2));
sigma_D_41_3 = nu_0_41_3 * sqrt(k_B*T / (M_41*c^2));
sigma_D_41_4 = nu_0_41_4 * sqrt(k_B*T / (M_41*c^2));

%% 峰值处的吸收截面，1.15
%K39的4条线的峰值处的吸收截面
sigma_0_39_1 = (1 / (sqrt(2*pi)*sigma_D_39_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%公式1.15，K39第1条线峰值处的吸收截面
sigma_0_39_2 = (1 / (sqrt(2*pi)*sigma_D_39_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%公式1.15，K39第2条线峰值处的吸收截面
sigma_0_39_3 = (1 / (sqrt(2*pi)*sigma_D_39_3)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%公式1.15，K39第3条线峰值处的吸收截面
sigma_0_39_4 = (1 / (sqrt(2*pi)*sigma_D_39_4)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%公式1.15，K39第4条线峰值处的吸收截面

%K41的4条线的峰值处的吸收截面
sigma_0_41_1 = (1 / (sqrt(2*pi)*sigma_D_41_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%公式1.15，K41第1条线峰值处的吸收截面
sigma_0_41_2 = (1 / (sqrt(2*pi)*sigma_D_41_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%公式1.15，K41第2条线峰值处的吸收截面
sigma_0_41_3 = (1 / (sqrt(2*pi)*sigma_D_41_3)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%公式1.15，K41第3条线峰值处的吸收截面
sigma_0_41_4 = (1 / (sqrt(2*pi)*sigma_D_41_4)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%公式1.15，K41第4条线峰值处的吸收截面

%% 绘图精度定义区，不绘图时记得注释掉

%nu = nu_L-2e9:1e6:nu_L+2e9;%绘图时，将nu变量换成等距向量，下公式要点乘

%% 单个吸收截面，1.13
%K39的4条线的吸收截面
sigma_abs_39_1 = sigma_0_39_1 * exp(-(nu-nu_0_39_1).^2 / (2*sigma_D_39_1^2));%公式1.13，K39第1条线吸收截面
sigma_abs_39_2 = sigma_0_39_2 * exp(-(nu-nu_0_39_2).^2 / (2*sigma_D_39_2^2));%公式1.13，K39第2条线吸收截面
sigma_abs_39_3 = sigma_0_39_3 * exp(-(nu-nu_0_39_3).^2 / (2*sigma_D_39_3^2));%公式1.13，K39第3条线吸收截面
sigma_abs_39_4 = sigma_0_39_4 * exp(-(nu-nu_0_39_4).^2 / (2*sigma_D_39_4^2));%公式1.13，K39第4条线吸收截面

%K41的4条线的吸收截面
sigma_abs_41_1 = sigma_0_41_1 * exp(-(nu-nu_0_41_1).^2 / (2*sigma_D_41_1^2));%公式1.13，K41第1条线吸收截面
sigma_abs_41_2 = sigma_0_41_2 * exp(-(nu-nu_0_41_2).^2 / (2*sigma_D_41_2^2));%公式1.13，K41第2条线吸收截面
sigma_abs_41_3 = sigma_0_41_3 * exp(-(nu-nu_0_41_3).^2 / (2*sigma_D_41_3^2));%公式1.13，K41第3条线吸收截面
sigma_abs_41_4 = sigma_0_41_4 * exp(-(nu-nu_0_41_4).^2 / (2*sigma_D_41_4^2));%公式1.13，K41第4条线吸收截面

%% K39的4个吸收截面叠加绘图区，不绘图时注释掉收起来
%{
plot(nu-nu_L,sigma_abs_39_1,'.-b', 'Linewidth', 2);%绘制K39第1条线的吸收截面
hold on
plot(nu-nu_L,sigma_abs_39_2,'.-r', 'Linewidth', 2);%绘制K39第2条线的吸收截面
hold on
plot(nu-nu_L,sigma_abs_39_3,'.-g', 'Linewidth', 2);%绘制K39第3条线的吸收截面
hold on
plot(nu-nu_L,sigma_abs_39_4,'.-m', 'Linewidth', 2);%绘制K39第4条线的吸收截面

max_sigma_abs_39_1 = max(sigma_abs_39_1);
max_sigma_abs_39_2 = max(sigma_abs_39_2);
max_sigma_abs_39_3 = max(sigma_abs_39_3);
max_sigma_abs_39_4 = max(sigma_abs_39_4);

hold on
%记住画直线的起终点是[x1 x2], [y1 y2],也就是先定义两个x，再定义两个y，而不是[x1 y1], [x2 y2]
plot([0.310e9,0.310e9],[0,max_sigma_abs_39_1],'.-b')
hold on
plot([0.254e9,0.254e9],[0,max_sigma_abs_39_2],'.-r')
hold on
plot([-0.152e9,-0.152e9],[0,max_sigma_abs_39_3],'.-g')
hold on
plot([-0.208e9,-0.208e9],[0,max_sigma_abs_39_4],'.-m')

set(gca,'FontSize',17);
xlabel('Frequency offset(Hz)');
ylabel('Absorption cross-section(m^2)');
legend('K39 line 1','K39 line 2','K39 line 3','K39 line 4');
%}

%% K41的4个吸收截面叠加绘图区，不绘图时注释掉收起来
%{
plot(nu-nu_L,sigma_abs_41_1,'.-b', 'Linewidth', 2);%绘制K41第1条线的吸收截面
hold on
plot(nu-nu_L,sigma_abs_41_2,'.-r', 'Linewidth', 2);%绘制K41第2条线的吸收截面
hold on
plot(nu-nu_L,sigma_abs_41_3,'.-g', 'Linewidth', 2);%绘制K41第3条线的吸收截面
hold on
plot(nu-nu_L,sigma_abs_41_4,'.-m', 'Linewidth', 2);%绘制K41第4条线的吸收截面

max_sigma_abs_41_1 = max(sigma_abs_41_1);
max_sigma_abs_41_2 = max(sigma_abs_41_2);
max_sigma_abs_41_3 = max(sigma_abs_41_3);
max_sigma_abs_41_4 = max(sigma_abs_41_4);

hold on
%记住画直线的起终点是[x1 x2], [y1 y2],也就是先定义两个x，再定义两个y，而不是[x1 y1], [x2 y2]
plot([0.405e9,0.405e9],[0,max_sigma_abs_41_1],'-b', 'Linewidth', 1)
hold on
plot([0.375e9,0.375e9],[0,max_sigma_abs_41_2],'-r', 'Linewidth', 1)
hold on
plot([0.151e9,0.151e9],[0,max_sigma_abs_41_3],'-g', 'Linewidth', 1)
hold on
plot([0.121e9,0.121e9],[0,max_sigma_abs_41_4],'-m', 'Linewidth', 1)

set(gca,'FontSize',17);
xlabel('Frequency offset(Hz)');
ylabel('Absorption cross-section(m^2)');
legend('K41 line 1','K41 line 2','K41 line 3','K41 line 4');
%}

%% K39的4条吸收截面乘以超精细跃迁强度因子绘图区
sigma_abs_39_1_fi = sigma_abs_39_1*f1;
sigma_abs_39_2_fi = sigma_abs_39_2*f2;
sigma_abs_39_3_fi = sigma_abs_39_3*f3;
sigma_abs_39_4_fi = sigma_abs_39_4*f4;

%{
plot(nu-nu_L,sigma_abs_39_1_fi,'.-b', 'Linewidth', 2);%绘制K39第1条线的吸收截面乘以超精细跃迁强度因子绘图区
hold on
plot(nu-nu_L,sigma_abs_39_2_fi,'.-r', 'Linewidth', 2);%绘制K39第2条线的吸收截面乘以超精细跃迁强度因子绘图区
hold on
plot(nu-nu_L,sigma_abs_39_3_fi,'.-g', 'Linewidth', 2);%绘制K39第3条线的吸收截面乘以超精细跃迁强度因子绘图区
hold on
plot(nu-nu_L,sigma_abs_39_4_fi,'.-m', 'Linewidth', 2);%绘制K39第4条线的吸收截面乘以超精细跃迁强度因子绘图区

max_sigma_abs_39_1_fi = max(sigma_abs_39_1_fi);
max_sigma_abs_39_2_fi = max(sigma_abs_39_2_fi);
max_sigma_abs_39_3_fi = max(sigma_abs_39_3_fi);
max_sigma_abs_39_4_fi = max(sigma_abs_39_4_fi);

hold on
%记住画直线的起终点是[x1 x2], [y1 y2],也就是先定义两个x，再定义两个y，而不是[x1 y1], [x2 y2]
plot([0.310e9,0.310e9],[0,max_sigma_abs_39_1_fi],'.-b')
hold on
plot([0.254e9,0.254e9],[0,max_sigma_abs_39_2_fi],'.-r')
hold on
plot([-0.152e9,-0.152e9],[0,max_sigma_abs_39_3_fi],'.-g')
hold on
plot([-0.208e9,-0.208e9],[0,max_sigma_abs_39_4_fi],'.-m')

set(gca,'FontSize',17);
xlabel('Frequency offset(Hz)');
ylabel('Absorption cross-section(m^2)');
legend('K39 line 1','K39 line 2','K39 line 3','K39 line 4');
%}

%% K41的4条吸收截面乘以超精细跃迁强度因子绘图区
sigma_abs_41_1_fi = sigma_abs_41_1*f1;
sigma_abs_41_2_fi = sigma_abs_41_2*f2;
sigma_abs_41_3_fi = sigma_abs_41_3*f3;
sigma_abs_41_4_fi = sigma_abs_41_4*f4;

%{
plot(nu-nu_L,sigma_abs_41_1_fi,'.-b', 'Linewidth', 2);%绘制K41第1条线的吸收截面乘以超精细跃迁强度因子绘图区
hold on
plot(nu-nu_L,sigma_abs_41_2_fi,'.-r', 'Linewidth', 2);%绘制K41第2条线的吸收截面乘以超精细跃迁强度因子绘图区
hold on
plot(nu-nu_L,sigma_abs_41_3_fi,'.-g', 'Linewidth', 2);%绘制K41第3条线的吸收截面乘以超精细跃迁强度因子绘图区
hold on
plot(nu-nu_L,sigma_abs_41_4_fi,'.-m', 'Linewidth', 2);%绘制K41第4条线的吸收截面乘以超精细跃迁强度因子绘图区

max_sigma_abs_41_1_fi = max(sigma_abs_41_1_fi);
max_sigma_abs_41_2_fi = max(sigma_abs_41_2_fi);
max_sigma_abs_41_3_fi = max(sigma_abs_41_3_fi);
max_sigma_abs_41_4_fi = max(sigma_abs_41_4_fi);

hold on
%记住画直线的起终点是[x1 x2], [y1 y2],也就是先定义两个x，再定义两个y，而不是[x1 y1], [x2 y2]
plot([0.405e9,0.405e9],[0,max_sigma_abs_41_1_fi],'.-b')
hold on
plot([0.375e9,0.375e9],[0,max_sigma_abs_41_2_fi],'.-r')
hold on
plot([0.151e9,0.151e9],[0,max_sigma_abs_41_3_fi],'.-g')
hold on
plot([0.121e9,0.121e9],[0,max_sigma_abs_41_4_fi],'.-m')

set(gca,'FontSize',17);
xlabel('Frequency offset(Hz)');
ylabel('Absorption cross-section(m^2)');
legend('K41 line 1','K41 line 2','K41 line 3','K41 line 4');
%}

%% K39与K41的吸收截面叠加绘图区
sigma_abs_39 = IsoAbdn_39 * (sigma_abs_39_1_fi + sigma_abs_39_2_fi + sigma_abs_39_3_fi + sigma_abs_39_4_fi);
sigma_abs_41 = IsoAbdn_41 * (sigma_abs_41_1_fi + sigma_abs_41_2_fi + sigma_abs_41_3_fi + sigma_abs_41_4_fi);

%{
hold on
box on
plot(nu-nu_L,sigma_abs_39,'-b', 'Linewidth', 2);%绘制K39的吸收截面
hold on
plot(nu-nu_L,sigma_abs_41,'.-r', 'Linewidth', 2);%绘制K39的吸收截面

max_sigma_abs_39 = max(sigma_abs_39);
max_sigma_abs_41 = max(sigma_abs_41);

hold on
%记住画直线的起终点是[x1 x2], [y1 y2],也就是先定义两个x，再定义两个y，而不是[x1 y1], [x2 y2]
plot([-0.096e9,-0.096e9],[0,max_sigma_abs_39],'.-b')
hold on
plot([0.225e9,0.225e9],[0,max_sigma_abs_41],'.-r')

set(gca,'FontSize',17);
xlabel('Frequency offset(Hz)');
ylabel('Absorption cross-section(m^2)');
legend('K39','K41');
%}

%% 钾原子的吸收截面与散射截面叠加绘图区
sigma_abs = sigma_abs_39 + sigma_abs_41;

%fx = max(sigma_abs)/(4*pi)
%{
sigma_eff_fx = conv(sigma_abs,g_L);
sigma_eff_fx = sigma_eff_fx * 1e6;
X = nu_L-2e9:0.5e6:nu_L+2e9;
Y = sigma_eff_fx;

%
plot(nu-nu_L,sigma_abs,'-b','Linewidth',2);%绘制钾原子的吸收截面
hold on
box on
plot((X-nu_L)*2,Y,'-r','Linewidth',2);%绘制钾原子的散射截面
xlim([-2e9,2e9]);
set(gca,'FontSize',19);
xlabel('Frequency offset(Hz)');
ylabel('Cross-section(m^2)');
legend('Absorption cross-section','Effective scattering cross-section');
%}

%% 有效后向散射截面，画图时注释掉收起来,计算时记得关闭绘图精度

sigma_eff_fx = sigma_abs.*g_L;
sigma_eff = int(sigma_eff_fx,nu,-inf,inf);
sigma_K = sigma_eff / (4*pi);
sprintf('钾原子的D1线的有效后向散射截面为\n%e m^2/sr',sigma_K)
%{
plot(sigma_FWHM,sigma_K,'-r','Linewidth',2)
set(gca,'FontSize',18);
xlabel('Line width(Hz)');
ylabel('Cross-section(m^2/sr)');

%}