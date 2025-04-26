%{
��Ŀ���ƣ���ԭ��D2�ߵ���Ч����ɢ�����������
�������ƣ�������
���ߵ�λ�����տƼ���ѧ��ѧԺ
�������䣺wangkexin1998@qq.com
��ʦ���ƣ�������
�汾��ţ�MATLAB-3
����ʱ�䣺2022/09/15/14:31:51
�޸�ʱ�䣺2022/09/18/10:58:07
%}

%% ������ֵ��
%ͨ�ó���
c = 2.99792458e8;%����(m/s)
V = 150;
k_B = 1.3806505e-23;%������������(J/K)
e = 1.60217662e-19;%���ӵ���(C)
m_e = 9.10938215e-31;%��������(kg)
epsilon_0 = 8.854187817e-12;%��յ�����(F/m)

%���ⳣ��  589.1583e-9
lambda_L = 589.1576e-9;%���Ⲩ��(m),������������β���ļ����״�۲��о����Ʋ�ȫ��ʱ�۲⼼��,P12��
nu_L = (c) / lambda_L;%��������Ƶ��(Hz)
%ר�ó���
T = 200;%�¶�Ϊ200K�����¡�����������β���ļ����״�۲��о����Ʋ�ȫ��ʱ�۲⼼����P14��
M = 3.82e-26;%Naԭ�Ӿ�������(kg)���������ԭ����������ó���
fD2 = 0.6408;%Naԭ��D2�ߵ�����ǿ�ȡ��߿��Ʋ㡢�ز�ͬʱ̽��ļ����״
IsoAbdn = 1;%����Naû��ͬλ�أ����ȡ1
sigma_FWHM = 80e6;%��������߿�1.5GHz�����ȫ���߿��Ʋ㡢�ز�ͬʱ̽��ļ����״
sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));%��˹���͵�RMS���

% ���߲�����ֵ��
%Na��6��ԾǨ���ߵ�����Ƶ��ƫ��ֵ�����ţ���Ϊ�������Ͳ����������ǹ��ס� + 0.6302e9
nu01 = (nu_L - 0.7328e9) * ((c+V)/c);%Na��1��ԾǨ���ߵ�����Ƶ�ʡ�����������β���ļ����״�۲��о����Ʋ�ȫ��ʱ�۲⼼��,P12��
nu02 = (nu_L - 0.6962e9) * ((c+V)/c);%Na-2
nu03 = (nu_L - 0.6302e9) * ((c+V)/c);%Na-3
nu04 = (nu_L + 1.0333e9) * ((c+V)/c);%Na-4
nu05 = (nu_L + 1.0552e9) * ((c+V)/c);%Na-5
nu06 = (nu_L + 1.0919e9) * ((c+V)/c);%Na-6
%Na�ĳ���ϸԾǨǿ�����ӡ�����������β���ļ����״�۲��о����Ʋ�ȫ��ʱ�۲⼼��,P12��
f1 = 1/32;
f2 = 5/32;
f3 = 14/32;
f4 = 2/32;
f5 = 5/32;
f6 = 5/32;

% ��һ�����ͣ�1.11

syms nu;%�������Ƶ��
nu = nu_L-2.5e9:1e7:nu_L+2.5e9;
%��һ����������
g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));%��ʽ1.11

% �����������RMS��ȣ�1.14
%Na��6���ߵĶ����������RMS���
sigma_D_1 = nu01 * sqrt(k_B*T / (M*c^2));
sigma_D_2 = nu02 * sqrt(k_B*T / (M*c^2));
sigma_D_3 = nu03 * sqrt(k_B*T / (M*c^2));
sigma_D_4 = nu04 * sqrt(k_B*T / (M*c^2));
sigma_D_5 = nu05 * sqrt(k_B*T / (M*c^2));
sigma_D_6 = nu06 * sqrt(k_B*T / (M*c^2));

% ��ֵ�������ս��棬1.15
%Na��6���ߵķ�ֵ�������ս���
sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
sigma_0_2 = (1 / (sqrt(2*pi)*sigma_D_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
sigma_0_3 = (1 / (sqrt(2*pi)*sigma_D_3)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
sigma_0_4 = (1 / (sqrt(2*pi)*sigma_D_4)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
sigma_0_5 = (1 / (sqrt(2*pi)*sigma_D_5)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
sigma_0_6 = (1 / (sqrt(2*pi)*sigma_D_6)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;

% �������ս��棬1.13
%Na��6���ߵ����ս���
sigma_abs_1 = sigma_0_1 * exp(-(nu-nu01).^2 / (2*sigma_D_1^2));
sigma_abs_2 = sigma_0_2 * exp(-(nu-nu02).^2 / (2*sigma_D_2^2));
sigma_abs_3 = sigma_0_3 * exp(-(nu-nu03).^2 / (2*sigma_D_3^2));
sigma_abs_4 = sigma_0_4 * exp(-(nu-nu04).^2 / (2*sigma_D_4^2));
sigma_abs_5 = sigma_0_5 * exp(-(nu-nu05).^2 / (2*sigma_D_5^2));
sigma_abs_6 = sigma_0_6 * exp(-(nu-nu06).^2 / (2*sigma_D_6^2));

% Na��6�����ս�����Գ���ϸԾǨǿ�����ӻ�ͼ��
sigma_abs_1_fi = sigma_abs_1*f1;
sigma_abs_2_fi = sigma_abs_2*f2;
sigma_abs_3_fi = sigma_abs_3*f3;
sigma_abs_4_fi = sigma_abs_4*f4;
sigma_abs_5_fi = sigma_abs_5*f5;
sigma_abs_6_fi = sigma_abs_6*f6;

% Na�����ս�����ӻ�ͼ��
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

%% ��Ч����ɢ����棬��ͼʱע�͵�������,����ʱ�ǵùرջ�ͼ����

sigma_eff_fx = sigma_abs.*g_L;
sigma_eff = int(sigma_eff_fx,nu,-inf,inf);
sigma_K = sigma_eff / (4*pi);
sprintf('Naԭ�ӵ�D2�ߵ���Ч����ɢ�����Ϊ\n%e m^2/sr',sigma_K)
