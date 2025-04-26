%{
��Ŀ���ƣ���ԭ��D1�ߵ���Ч����ɢ�����������
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
lambda_L = 769.898e-9;%���Ⲩ��(m)
nu_L = c / lambda_L;%��������Ƶ��(Hz)
%ר�ó���
T = 200;%�¶�Ϊ200K������
M_39 = 6.4679e-26;%K39ԭ�Ӿ�������(kg)
M_41 = 6.7996e-26;%K41ԭ�Ӿ�������(kg)
fD1 = 0.3327;%��ԭ��D1�ߵ�����ǿ��
IsoAbdn_39 = 0.933;%K39��ͬλ�ط��
IsoAbdn_41 = 0.067;%K39��ͬλ�ط��
sigma_FWHM =1.5e9;%��������߿�1.44GHz�����ȫ��
sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));%��˹���͵�RMS���

%% ���߲�����ֵ��
%K39��4��ԾǨ���ߵ�����Ƶ��ƫ��ֵ
nu_0_39_1 = nu_L+0.310e9;%K39��1��ԾǨ���ߵ�����Ƶ��
nu_0_39_2 = nu_L+0.254e9;%K39��2��ԾǨ���ߵ�����Ƶ��
nu_0_39_3 = nu_L-0.152e9;%K39��3��ԾǨ���ߵ�����Ƶ��
nu_0_39_4 = nu_L-0.208e9;%K39��4��ԾǨ���ߵ�����Ƶ��
%K41��4��ԾǨ���ߵ�����Ƶ��ƫ��ֵ
nu_0_41_1 = nu_L+0.405e9;%K41��1��ԾǨ���ߵ�����Ƶ��
nu_0_41_2 = nu_L+0.375e9;%K41��2��ԾǨ���ߵ�����Ƶ��
nu_0_41_3 = nu_L+0.151e9;%K41��3��ԾǨ���ߵ�����Ƶ��
nu_0_41_4 = nu_L+0.121e9;%K41��4��ԾǨ���ߵ�����Ƶ��
%K�ĳ���ϸԾǨǿ������
f1 = 5/16;
f2 = 1/16;
f3 = 5/16;
f4 = 5/16;

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
plot([-0.75e9 0.75e9],[3.1e-10 3.1e-10],'.-k')
set(gca,'FontSize',22);
xlabel('Frequency offset(Hz)');
ylabel('Normalized intensity');
%}


%% �����������RMS��ȣ�1.14
%K39��4���ߵĶ����������RMS���
sigma_D_39_1 = nu_0_39_1 * sqrt(k_B*T / (M_39*c^2));
sigma_D_39_2 = nu_0_39_2 * sqrt(k_B*T / (M_39*c^2));
sigma_D_39_3 = nu_0_39_3 * sqrt(k_B*T / (M_39*c^2));
sigma_D_39_4 = nu_0_39_4 * sqrt(k_B*T / (M_39*c^2));

%K41��4���ߵĶ����������RMS���
sigma_D_41_1 = nu_0_41_1 * sqrt(k_B*T / (M_41*c^2));
sigma_D_41_2 = nu_0_41_2 * sqrt(k_B*T / (M_41*c^2));
sigma_D_41_3 = nu_0_41_3 * sqrt(k_B*T / (M_41*c^2));
sigma_D_41_4 = nu_0_41_4 * sqrt(k_B*T / (M_41*c^2));

%% ��ֵ�������ս��棬1.15
%K39��4���ߵķ�ֵ�������ս���
sigma_0_39_1 = (1 / (sqrt(2*pi)*sigma_D_39_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%��ʽ1.15��K39��1���߷�ֵ�������ս���
sigma_0_39_2 = (1 / (sqrt(2*pi)*sigma_D_39_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%��ʽ1.15��K39��2���߷�ֵ�������ս���
sigma_0_39_3 = (1 / (sqrt(2*pi)*sigma_D_39_3)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%��ʽ1.15��K39��3���߷�ֵ�������ս���
sigma_0_39_4 = (1 / (sqrt(2*pi)*sigma_D_39_4)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%��ʽ1.15��K39��4���߷�ֵ�������ս���

%K41��4���ߵķ�ֵ�������ս���
sigma_0_41_1 = (1 / (sqrt(2*pi)*sigma_D_41_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%��ʽ1.15��K41��1���߷�ֵ�������ս���
sigma_0_41_2 = (1 / (sqrt(2*pi)*sigma_D_41_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%��ʽ1.15��K41��2���߷�ֵ�������ս���
sigma_0_41_3 = (1 / (sqrt(2*pi)*sigma_D_41_3)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%��ʽ1.15��K41��3���߷�ֵ�������ս���
sigma_0_41_4 = (1 / (sqrt(2*pi)*sigma_D_41_4)) * (e^2 / (4*epsilon_0*m_e*c)) * fD1;
%��ʽ1.15��K41��4���߷�ֵ�������ս���

%% ��ͼ���ȶ�����������ͼʱ�ǵ�ע�͵�

%nu = nu_L-2e9:1e6:nu_L+2e9;%��ͼʱ����nu�������ɵȾ��������¹�ʽҪ���

%% �������ս��棬1.13
%K39��4���ߵ����ս���
sigma_abs_39_1 = sigma_0_39_1 * exp(-(nu-nu_0_39_1).^2 / (2*sigma_D_39_1^2));%��ʽ1.13��K39��1�������ս���
sigma_abs_39_2 = sigma_0_39_2 * exp(-(nu-nu_0_39_2).^2 / (2*sigma_D_39_2^2));%��ʽ1.13��K39��2�������ս���
sigma_abs_39_3 = sigma_0_39_3 * exp(-(nu-nu_0_39_3).^2 / (2*sigma_D_39_3^2));%��ʽ1.13��K39��3�������ս���
sigma_abs_39_4 = sigma_0_39_4 * exp(-(nu-nu_0_39_4).^2 / (2*sigma_D_39_4^2));%��ʽ1.13��K39��4�������ս���

%K41��4���ߵ����ս���
sigma_abs_41_1 = sigma_0_41_1 * exp(-(nu-nu_0_41_1).^2 / (2*sigma_D_41_1^2));%��ʽ1.13��K41��1�������ս���
sigma_abs_41_2 = sigma_0_41_2 * exp(-(nu-nu_0_41_2).^2 / (2*sigma_D_41_2^2));%��ʽ1.13��K41��2�������ս���
sigma_abs_41_3 = sigma_0_41_3 * exp(-(nu-nu_0_41_3).^2 / (2*sigma_D_41_3^2));%��ʽ1.13��K41��3�������ս���
sigma_abs_41_4 = sigma_0_41_4 * exp(-(nu-nu_0_41_4).^2 / (2*sigma_D_41_4^2));%��ʽ1.13��K41��4�������ս���

%% K39��4�����ս�����ӻ�ͼ��������ͼʱע�͵�������
%{
plot(nu-nu_L,sigma_abs_39_1,'.-b', 'Linewidth', 2);%����K39��1���ߵ����ս���
hold on
plot(nu-nu_L,sigma_abs_39_2,'.-r', 'Linewidth', 2);%����K39��2���ߵ����ս���
hold on
plot(nu-nu_L,sigma_abs_39_3,'.-g', 'Linewidth', 2);%����K39��3���ߵ����ս���
hold on
plot(nu-nu_L,sigma_abs_39_4,'.-m', 'Linewidth', 2);%����K39��4���ߵ����ս���

max_sigma_abs_39_1 = max(sigma_abs_39_1);
max_sigma_abs_39_2 = max(sigma_abs_39_2);
max_sigma_abs_39_3 = max(sigma_abs_39_3);
max_sigma_abs_39_4 = max(sigma_abs_39_4);

hold on
%��ס��ֱ�ߵ����յ���[x1 x2], [y1 y2],Ҳ�����ȶ�������x���ٶ�������y��������[x1 y1], [x2 y2]
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

%% K41��4�����ս�����ӻ�ͼ��������ͼʱע�͵�������
%{
plot(nu-nu_L,sigma_abs_41_1,'.-b', 'Linewidth', 2);%����K41��1���ߵ����ս���
hold on
plot(nu-nu_L,sigma_abs_41_2,'.-r', 'Linewidth', 2);%����K41��2���ߵ����ս���
hold on
plot(nu-nu_L,sigma_abs_41_3,'.-g', 'Linewidth', 2);%����K41��3���ߵ����ս���
hold on
plot(nu-nu_L,sigma_abs_41_4,'.-m', 'Linewidth', 2);%����K41��4���ߵ����ս���

max_sigma_abs_41_1 = max(sigma_abs_41_1);
max_sigma_abs_41_2 = max(sigma_abs_41_2);
max_sigma_abs_41_3 = max(sigma_abs_41_3);
max_sigma_abs_41_4 = max(sigma_abs_41_4);

hold on
%��ס��ֱ�ߵ����յ���[x1 x2], [y1 y2],Ҳ�����ȶ�������x���ٶ�������y��������[x1 y1], [x2 y2]
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

%% K39��4�����ս�����Գ���ϸԾǨǿ�����ӻ�ͼ��
sigma_abs_39_1_fi = sigma_abs_39_1*f1;
sigma_abs_39_2_fi = sigma_abs_39_2*f2;
sigma_abs_39_3_fi = sigma_abs_39_3*f3;
sigma_abs_39_4_fi = sigma_abs_39_4*f4;

%{
plot(nu-nu_L,sigma_abs_39_1_fi,'.-b', 'Linewidth', 2);%����K39��1���ߵ����ս�����Գ���ϸԾǨǿ�����ӻ�ͼ��
hold on
plot(nu-nu_L,sigma_abs_39_2_fi,'.-r', 'Linewidth', 2);%����K39��2���ߵ����ս�����Գ���ϸԾǨǿ�����ӻ�ͼ��
hold on
plot(nu-nu_L,sigma_abs_39_3_fi,'.-g', 'Linewidth', 2);%����K39��3���ߵ����ս�����Գ���ϸԾǨǿ�����ӻ�ͼ��
hold on
plot(nu-nu_L,sigma_abs_39_4_fi,'.-m', 'Linewidth', 2);%����K39��4���ߵ����ս�����Գ���ϸԾǨǿ�����ӻ�ͼ��

max_sigma_abs_39_1_fi = max(sigma_abs_39_1_fi);
max_sigma_abs_39_2_fi = max(sigma_abs_39_2_fi);
max_sigma_abs_39_3_fi = max(sigma_abs_39_3_fi);
max_sigma_abs_39_4_fi = max(sigma_abs_39_4_fi);

hold on
%��ס��ֱ�ߵ����յ���[x1 x2], [y1 y2],Ҳ�����ȶ�������x���ٶ�������y��������[x1 y1], [x2 y2]
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

%% K41��4�����ս�����Գ���ϸԾǨǿ�����ӻ�ͼ��
sigma_abs_41_1_fi = sigma_abs_41_1*f1;
sigma_abs_41_2_fi = sigma_abs_41_2*f2;
sigma_abs_41_3_fi = sigma_abs_41_3*f3;
sigma_abs_41_4_fi = sigma_abs_41_4*f4;

%{
plot(nu-nu_L,sigma_abs_41_1_fi,'.-b', 'Linewidth', 2);%����K41��1���ߵ����ս�����Գ���ϸԾǨǿ�����ӻ�ͼ��
hold on
plot(nu-nu_L,sigma_abs_41_2_fi,'.-r', 'Linewidth', 2);%����K41��2���ߵ����ս�����Գ���ϸԾǨǿ�����ӻ�ͼ��
hold on
plot(nu-nu_L,sigma_abs_41_3_fi,'.-g', 'Linewidth', 2);%����K41��3���ߵ����ս�����Գ���ϸԾǨǿ�����ӻ�ͼ��
hold on
plot(nu-nu_L,sigma_abs_41_4_fi,'.-m', 'Linewidth', 2);%����K41��4���ߵ����ս�����Գ���ϸԾǨǿ�����ӻ�ͼ��

max_sigma_abs_41_1_fi = max(sigma_abs_41_1_fi);
max_sigma_abs_41_2_fi = max(sigma_abs_41_2_fi);
max_sigma_abs_41_3_fi = max(sigma_abs_41_3_fi);
max_sigma_abs_41_4_fi = max(sigma_abs_41_4_fi);

hold on
%��ס��ֱ�ߵ����յ���[x1 x2], [y1 y2],Ҳ�����ȶ�������x���ٶ�������y��������[x1 y1], [x2 y2]
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

%% K39��K41�����ս�����ӻ�ͼ��
sigma_abs_39 = IsoAbdn_39 * (sigma_abs_39_1_fi + sigma_abs_39_2_fi + sigma_abs_39_3_fi + sigma_abs_39_4_fi);
sigma_abs_41 = IsoAbdn_41 * (sigma_abs_41_1_fi + sigma_abs_41_2_fi + sigma_abs_41_3_fi + sigma_abs_41_4_fi);

%{
hold on
box on
plot(nu-nu_L,sigma_abs_39,'-b', 'Linewidth', 2);%����K39�����ս���
hold on
plot(nu-nu_L,sigma_abs_41,'.-r', 'Linewidth', 2);%����K39�����ս���

max_sigma_abs_39 = max(sigma_abs_39);
max_sigma_abs_41 = max(sigma_abs_41);

hold on
%��ס��ֱ�ߵ����յ���[x1 x2], [y1 y2],Ҳ�����ȶ�������x���ٶ�������y��������[x1 y1], [x2 y2]
plot([-0.096e9,-0.096e9],[0,max_sigma_abs_39],'.-b')
hold on
plot([0.225e9,0.225e9],[0,max_sigma_abs_41],'.-r')

set(gca,'FontSize',17);
xlabel('Frequency offset(Hz)');
ylabel('Absorption cross-section(m^2)');
legend('K39','K41');
%}

%% ��ԭ�ӵ����ս�����ɢ�������ӻ�ͼ��
sigma_abs = sigma_abs_39 + sigma_abs_41;

%fx = max(sigma_abs)/(4*pi)
%{
sigma_eff_fx = conv(sigma_abs,g_L);
sigma_eff_fx = sigma_eff_fx * 1e6;
X = nu_L-2e9:0.5e6:nu_L+2e9;
Y = sigma_eff_fx;

%
plot(nu-nu_L,sigma_abs,'-b','Linewidth',2);%���Ƽ�ԭ�ӵ����ս���
hold on
box on
plot((X-nu_L)*2,Y,'-r','Linewidth',2);%���Ƽ�ԭ�ӵ�ɢ�����
xlim([-2e9,2e9]);
set(gca,'FontSize',19);
xlabel('Frequency offset(Hz)');
ylabel('Cross-section(m^2)');
legend('Absorption cross-section','Effective scattering cross-section');
%}

%% ��Ч����ɢ����棬��ͼʱע�͵�������,����ʱ�ǵùرջ�ͼ����

sigma_eff_fx = sigma_abs.*g_L;
sigma_eff = int(sigma_eff_fx,nu,-inf,inf);
sigma_K = sigma_eff / (4*pi);
sprintf('��ԭ�ӵ�D1�ߵ���Ч����ɢ�����Ϊ\n%e m^2/sr',sigma_K)
%{
plot(sigma_FWHM,sigma_K,'-r','Linewidth',2)
set(gca,'FontSize',18);
xlabel('Line width(Hz)');
ylabel('Cross-section(m^2/sr)');

%}