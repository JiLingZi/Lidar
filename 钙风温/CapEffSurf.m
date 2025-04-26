%% �����ӷ���������

% E���¶� 200 K ; F���¶� 950 K

%% E��ɢ�����

% ͨ�ó���
c = 2.99792458e8;                           %����(m/s)
k_B = 1.3806505e-23;                        %������������(J/K)
e = 1.60217662e-19;                         %���ӵ���(C)
m_e = 9.10938215e-31;                       %��������(kg)
epsilon_0 = 8.854187817e-12;                %��յ�����(F/m)

% ���ⳣ��
lambda_L = 393.3663e-9;                     %���Ⲩ��(m)
nu_L = c / lambda_L;                        %��������Ƶ��(Hz)
sigma_FWHM = 180e6;                         %�����߿�
sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));  %��˹���͵�RMS���

% ƫ����
DnuR = 700;
DnuL = -700;

% �����ӳ���
M = 6.665e-26;                              %�����Ӿ�������(kg)
fD2 = 0.69;                                 %����������ǿ��

% �Ϳշ�������
T = 120:1:320;
V = -120:1:120;

% ��Ч����ɢ�����ѭ������
Data_eff_0 = zeros(size(T,2),size(V,2));
Data_eff_R = zeros(size(T,2),size(V,2));
Data_eff_L = zeros(size(T,2),size(V,2));

% �������ѭ��
for j = 1:length(V)
    
    nu_0_1 = nu_L * ((c+V(j))/c);
    f1 = 1;
    nu = nu_L-3e9:1e6:nu_L+3e9;
    g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));
    
    % �����¶�ѭ��
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
    strj = string(j);
    disp(['E�����ǵ�',strj]);
end

% ������Ϸ��ٺ��¶�����
T = T';
T = repmat(T,1,size(V,2));
V = repmat(V,size(T,1),1);

% ��Ϸ��±��ʾ���
R_T = (Data_eff_R + Data_eff_L) ./ (2 * Data_eff_0);
R_V = (Data_eff_R - Data_eff_L) ./ (Data_eff_0);

% ���¶Ƚ��ж�ά�������
ft = fittype( 'poly55' );
[fit_T, gof_T] = fit( [R_V(:), R_T(:)], T(:), ft );
coe_fit_t=coeffvalues(fit_T);
f_TL = fit_T;

% �����¶��������
figure('Name','�¶��������');
plot(f_TL);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Fitting Surface of Temperature');
xlabel('R_{V}');
ylabel('R_{T}');
zlabel('Temperature (K)');

% �Է��ٽ��ж�ά�������
ft = fittype( 'poly55' );
[fit_V, gof_V] = fit( [R_V(:), R_T(:)], V(:), ft );
coe_fit_v=coeffvalues(fit_V);
f_VL = fit_V;

% ���������������
figure('Name','�����������');
plot(f_VL);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Fitting Surface of Wind Velocity');
xlabel('R_{V}');
ylabel('R_{T}');
zlabel('Wind Velocity (m/s)');


%% F��ɢ�����-�����������

% ͨ�ó���
c = 2.99792458e8;                           %����(m/s)
k_B = 1.3806505e-23;                        %������������(J/K)
e = 1.60217662e-19;                         %���ӵ���(C)
m_e = 9.10938215e-31;                       %��������(kg)
epsilon_0 = 8.854187817e-12;                %��յ�����(F/m)

% ���ⳣ��
lambda_L = 393.477469e-9;                     %���Ⲩ��(m)
nu_L = c / lambda_L;                        %��������Ƶ��(Hz)
sigma_FWHM = 180e6;                         %�����߿�
sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));  %��˹���͵�RMS���

% ƫ����
DnuR = 1500;
DnuL = -1500;

% �����ӳ���
M = 6.665e-26;                              %�����Ӿ�������(kg)
fD2 = 0.69;                                 %����������ǿ��

% �Ϳշ�������
T = 850:1:1050;
V = -120:1:120;

% ��Ч����ɢ�����ѭ������
Data_eff_0 = zeros(size(T,2),size(V,2));
Data_eff_R = zeros(size(T,2),size(V,2));
Data_eff_L = zeros(size(T,2),size(V,2));

% �������ѭ��
for j = 1:length(V)
    
    nu_0_1 = nu_L * ((c+V(j))/c);
    f1 = 1;
    nu = nu_L-3e9:1e6:nu_L+3e9;
    g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));
    
    % �����¶�ѭ��
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
        disp(['F������',j,':',i]);
    end
    strj = string(j);
    disp(['F�����ǵ�',strj]);
end

% ������Ϸ��ٺ��¶�����
T = T';
T = repmat(T,1,size(V,2));
V = repmat(V,size(T,1),1);

% ��Ϸ��±��ʾ���
R_T = (Data_eff_R + Data_eff_L) ./ (2 * Data_eff_0);
R_V = (Data_eff_R - Data_eff_L) ./ (Data_eff_0);

% ���¶Ƚ��ж�ά�������
ft = fittype( 'poly55' );
[fit_T, gof_T] = fit( [R_V(:), R_T(:)], T(:), ft );
coe_fit_t=coeffvalues(fit_T);
f_TH = fit_T;

% �����¶��������
figure('Name','�¶��������F');
plot(f_TH);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Fitting Surface of Temperature');
xlabel('R_{V}');
ylabel('R_{T}');
zlabel('Temperature (K)');

% �Է��ٽ��ж�ά�������
ft = fittype( 'poly55' );
[fit_V, gof_V] = fit( [R_V(:), R_T(:)], V(:), ft );
coe_fit_v=coeffvalues(fit_V);
f_VH = fit_V;

% ���������������
figure('Name','�����������F');
plot(f_VH);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Fitting Surface of Wind Velocity');
xlabel('R_{V}');
ylabel('R_{T}');
zlabel('Wind Velocity (m/s)');


