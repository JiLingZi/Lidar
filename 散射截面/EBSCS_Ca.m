%{
��Ŀ���ƣ���ԭ�ӵ���Ч����ɢ�����������
�������ƣ�������
���ߵ�λ�����տƼ���ѧ��ѧԺ
�������䣺wangkexin1998@qq.com
��ʦ���ƣ�������
�汾��ţ�MATLAB-3
����ʱ�䣺2022/09/14/14:31:51
�޸�ʱ�䣺2022/09/15/10:58:07
%}

%% ������ֵ��
%ͨ�ó���
c = 2.99792458e8;%����(m/s)
k_B = 1.3806505e-23;%������������(J/K)
e = 1.60217662e-19;%���ӵ���(C)
m_e = 9.10938215e-31;%��������(kg)
epsilon_0 = 8.854187817e-12;%��յ�����(F/m)

%���ⳣ��
lambda_L = 422.673e-9;%���Ⲩ��(m)
nu_L = c / lambda_L;%��������Ƶ��(Hz)

%ר�ó���
T = 200;%�¶�Ϊ200K������
M = 6.665e-26;%ԭ�Ӿ�������(kg)
fD1 = 1.75;%����ǿ��
IsoAbdn = 1;%Fe��ͬλ�ط��
sigma_FWHM = 0.3e9;%��������߿����ȫ��
sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));%��˹���͵�RMS���

%% ���߲�����ֵ��
%ԾǨ���ߵ�����Ƶ��ƫ��ֵ
nu_0_1 =c/422.673e-9;
%����ϸԾǨǿ������
f1 = 1;

%% ��һ�����ͣ�1.11

syms nu;%�������Ƶ��

%��һ���������ͻ�ͼ���ȣ���ͼ�����ڼ���ʱ����ر�
%nu = nu_L-2e9:1e6:nu_L+2e9;%��ͼʱ����nu�������ɵȾ��������¹�ʽҪ���


g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));%��ʽ1.11
%{
int_g_L = int(g_L,nu,-inf,inf);%��g_L��-oo��+oo������֤��һ�����ó�ֵΪ1
sprintf('��һ����֤\n%e',int_g_L)
%}
%{
plot(nu-nu_L,g_L,'.-k');%���ƹ�һ���ļ�������
max_g_L = max(g_L);
hold on
plot([0 0], [0 max_g_L],'-k')
set(gca,'FontSize',20);
xlabel('Frequency offset(Hz)');
ylabel('Normalized intensity');
%}


%% �����������RMS��ȣ�1.14
%Fe��1���ߵĶ����������RMS���
sigma_D_1 = nu_0_1 * sqrt(k_B*T / (M*c^2));
%sigma_D_2 = nu_0_2 * sqrt(k_B*T / (M*c^2));

%% ��ֵ�������ս��棬1.15
%Fe��1���ߵķ�ֵ�������ս���
sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%sigma_0_2 = (1 / (sqrt(2*pi)*sigma_D_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;


%% ��ͼ���ȶ�����������ͼʱ�ǵ�ע�͵�

%nu = nu_L-2e9:1e5:nu_L+2e9;%��ͼʱ����nu�������ɵȾ��������¹�ʽҪ���

%% �������ս��棬1.13
%Fe��1���ߵ����ս���
sigma_abs_1 = sigma_0_1 * exp(-(nu-nu_0_1).^2 / (2*sigma_D_1^2));
%sigma_abs_2 = sigma_0_2 * exp(-(nu-nu_0_1).^2 / (2*sigma_D_2^2));

%% ���ս�����ӻ�ͼ��������ͼʱע�͵�������
%{
plot(nu-nu_L,sigma_abs_1,'.-b');%����K39��1���ߵ����ս���
hold on

max_sigma_abs_1 = max(sigma_abs_1);
hold on

%��ס��ֱ�ߵ����յ���[x1 x2], [y1 y2],Ҳ�����ȶ�������x���ٶ�������y��������[x1 y1], [x2 y2]
plot([c / 372.0993,c / 372.0993],[0,max_sigma_abs_1],'.-b')
hold on

set(gca,'FontSize',20);
xlabel('Frequency offset(Hz)');
ylabel('Absorption cross-section(m^2)');
%}

%% ���ս�����Գ���ϸԾǨǿ�����ӻ�ͼ��
sigma_abs_1_fi = sigma_abs_1*f1;
%sigma_abs_2_fi = sigma_abs_2*f2;

%{
plot(nu-nu_L,sigma_abs_1_fi,'.-b');%����K39��1���ߵ����ս���
hold on

max_sigma_abs_1_fi = max(sigma_abs_1_fi);
hold on

%��ס��ֱ�ߵ����յ���[x1 x2], [y1 y2],Ҳ�����ȶ�������x���ٶ�������y��������[x1 y1], [x2 y2]
plot([c / 372.0993,c / 372.0993],[0,max_sigma_abs_1_fi],'.-b')
hold on

set(gca,'FontSize',20);
xlabel('Frequency offset(Hz)');
ylabel('Absorption cross-section(m^2)');
%}


%% K39��K41�����ս�����ӻ�ͼ��
sigma_abs = IsoAbdn * sigma_abs_1_fi;
%sigma_abs_1 = IsoAbdn * (sigma_abs_1_fi);
%sigma_abs_2 = IsoAbdn * (sigma_abs_2_fi);

%fx = max(sigma_abs)/(4*pi)

%{
plot(nu-nu_L,sigma_abs,'.-b');%����K39�����ս���
hold on
plot(nu-nu_L,sigma_abs_41,'.-r');%����K39�����ս���

max_sigma_abs_39 = max(sigma_abs_39);
max_sigma_abs_41 = max(sigma_abs_41);

hold on
%��ס��ֱ�ߵ����յ���[x1 x2], [y1 y2],Ҳ�����ȶ�������x���ٶ�������y��������[x1 y1], [x2 y2]
plot([-0.096e9,-0.096e9],[0,max_sigma_abs_39],'.-b')
hold on
plot([0.225e9,0.225e9],[0,max_sigma_abs_41],'.-r')

set(gca,'FontSize',20);
xlabel('Ƶ��ƫ��(Hz)');
ylabel('k39��K41�ĵ����ս���(m^2)');
legend('K39�����ս���','K41�����ս���');
%}

%% ��ԭ�ӵ����ս�����ɢ�������ӻ�ͼ��
sigma_abs = sigma_abs;

%fx = max(sigma_abs)/(4*pi)
%{
sigma_eff_fx = conv(sigma_abs,g_L);
sigma_eff_fx = sigma_eff_fx * 1e6;
X = nu_L-2e9:0.5e6:nu_L+2e9;
Y = sigma_eff_fx;
%sigma_eff_fx = sigma_abs.*g_L;
%
plot(nu-nu_L,sigma_abs,'-b','Linewidth',2);%���Ƽ�ԭ�ӵ����ս���
hold on
box on
plot((X-nu_L)*2,Y,'.-r','Linewidth',2);%���Ƽ�ԭ�ӵ�ɢ�����
xlim([-2e9,2e9]);
set(gca,'FontSize',19);
xlabel('Frequency offset(Hz)');
ylabel('Cross-section(m^2)');
legend('Absorption','Effective scattering');
%}

%% ��Ч����ɢ����棬��ͼʱע�͵�������,����ʱ�ǵùرջ�ͼ����

sigma_eff_fx = sigma_abs.*g_L;
sigma_eff = int(sigma_eff_fx,nu,-inf,inf);
sigma_K = sigma_eff / (4*pi);
sprintf('Caԭ�ӵ���Ч����ɢ�����Ϊ\n%e m^2/sr',sigma_K)
%{
plot(sigma_FWHM,sigma_K,'-m','Linewidth',2)
set(gca,'FontSize',18);
xlabel('Line width(Hz)');
ylabel('Cross-section(m^2/sr)');

%}