%-------------------------------------------------------------------------
%=================�и߲�����Ʒ��¼����״���·��ݳ���======================
%-------------------------------------------------------------------------
% ��Ŀ���ƣ�Į����ԭ�ӷ��·��ݳ���
% �������ƣ�������
% ���ߵ�λ�����տƼ���ѧ
% ���ߵ�ʦ��������
% �������䣺wang.kexin@stu.just.edu.cn
% ����ʱ�䣺2023-07-09 17:30(UTC+8)
% �޸�ʱ�䣺2023-07-11 00:30(UTC+8)
% ������Դ�����繤�̶���Į�Ӽ����״�̨վ 
%===============================ע������==================================
% ��������ɷ��ݵ����Ʒ�����Ƶ�����ļ������跴�ݷ�����ʱ��仯���ѭ��
% ��������ΪӦ��ʹ�ã�ƽʱ������ʹ�ã�
% ����������һ����Ҫ5748�룬������ŵ��Թ���ʱ�䣬Ԥ�����ڴ�
% ƽʱ����ʹ�á�Na_Wind_Temperature�����ԡ�Na_Wind_Temperature_EFF��Ϊ����
% �мǲ�Ҫʹ�� clear all���˷�ʱ��
% �߾����ţ�bin_num = 8912
% ʱ��ֱ��ʣ�2 min
% �ռ�ֱ��ʣ�30 m
% �����߿�80 MHz
% ��Ƶ�л�Ƶ����:585 MHz
% V����ֱ����
% N���������춥���30��
% W�������򣨱�Ǵ��󣩣��춥���30��
% ����ģ�����ڿ��ټ��㣺
%                     ���꣺03��+060��04��+091��05��+121��06��+152��
% % % % 01��+000��02��+031��03��+059��04��+090��05��+120��06��+151��
% ���꣺07��+182��08��+213��09��+244��10��+274��11��+305��12��+335��
% % % % 07��+181��08��+212��09��+243��10��+273��11��+304��12��+334��
%-------------------------------------------------------------------------

%% ������ֵ��

% ͨ�ó���
c = 2.99792458e8;                           %����(m/s)
k_B = 1.3806505e-23;                        %������������(J/K)
e = 1.60217662e-19;                         %���ӵ���(C)
m_e = 9.10938215e-31;                       %��������(kg)
epsilon_0 = 8.854187817e-12;                %��յ�����(F/m)

% ���ⳣ��
lambda_L = 589.1583e-9;                     %���Ⲩ��(m)
nu_L = c / lambda_L;                        %��������Ƶ��(Hz)
sigma_FWHM = 80e6;                          %�����߿�
sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));  %��˹���͵�RMS���

% ��ԭ�ӳ���
M = 3.82e-26;                               %Naԭ�Ӿ�������(kg)
fD2 = 0.6408;                               %Naԭ��D2�ߵ�����ǿ��
%% ����������

% �¶�����
T = [80 100 120 140 160 180 200 220 240 260 280 300 320];

% ��������
V = [-120 -100 -80 -60 -40 -20 0 20 40 60 80 100 120];

%% F_0��Ч����ɢ��������ѭ��

Data_eff_0 = [];

% �������ѭ��
for j = 1:length(V)    
    %% ���߲�����ֵ��
    
    % Na��6��ԾǨ���ߵ�����Ƶ��ƫ��ֵ(Hz)
    nu01 = (nu_L - 0.7328e9 + 0.6302e9) * ((c+V(j))/c);
    nu02 = (nu_L - 0.6962e9 + 0.6302e9) * ((c+V(j))/c);
    nu03 = (nu_L - 0.6302e9 + 0.6302e9) * ((c+V(j))/c);
    nu04 = (nu_L + 1.0333e9 + 0.6302e9) * ((c+V(j))/c);
    nu05 = (nu_L + 1.0552e9 + 0.6302e9) * ((c+V(j))/c);
    nu06 = (nu_L + 1.0919e9 + 0.6302e9) * ((c+V(j))/c);
    
    % Na�ĳ���ϸԾǨǿ������
    f1 = 1/32;
    f2 = 5/32;
    f3 = 14/32;
    f4 = 2/32;
    f5 = 5/32;
    f6 = 5/32;
    
    %% ��һ������
    syms nu;
    g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));
    
    % �¶�ѭ��
    for i = 1:length(T)
        %% �����������RMS���
        sigma_D_1 = nu01 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_2 = nu02 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_3 = nu03 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_4 = nu04 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_5 = nu05 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_6 = nu06 * sqrt(k_B*T(i) / (M*c^2));
       
        %% ��ֵ�������ս���
        sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_2 = (1 / (sqrt(2*pi)*sigma_D_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_3 = (1 / (sqrt(2*pi)*sigma_D_3)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_4 = (1 / (sqrt(2*pi)*sigma_D_4)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_5 = (1 / (sqrt(2*pi)*sigma_D_5)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_6 = (1 / (sqrt(2*pi)*sigma_D_6)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        
        %% ����ԾǨ�ߵĵ������ս���
        sigma_abs_1 = sigma_0_1 * exp(-(nu-nu01).^2 / (2*sigma_D_1^2));
        sigma_abs_2 = sigma_0_2 * exp(-(nu-nu02).^2 / (2*sigma_D_2^2));
        sigma_abs_3 = sigma_0_3 * exp(-(nu-nu03).^2 / (2*sigma_D_3^2));
        sigma_abs_4 = sigma_0_4 * exp(-(nu-nu04).^2 / (2*sigma_D_4^2));
        sigma_abs_5 = sigma_0_5 * exp(-(nu-nu05).^2 / (2*sigma_D_5^2));
        sigma_abs_6 = sigma_0_6 * exp(-(nu-nu06).^2 / (2*sigma_D_6^2));
       
        %% �������ս�����Գ���ϸԾǨǿ������
        sigma_abs_1_fi = sigma_abs_1*f1;
        sigma_abs_2_fi = sigma_abs_2*f2;
        sigma_abs_3_fi = sigma_abs_3*f3;
        sigma_abs_4_fi = sigma_abs_4*f4;
        sigma_abs_5_fi = sigma_abs_5*f5;
        sigma_abs_6_fi = sigma_abs_6*f6;
       
        %% ���ս������
        sigma_abs = sigma_abs_1_fi + sigma_abs_2_fi + sigma_abs_3_fi + sigma_abs_4_fi + sigma_abs_5_fi + sigma_abs_6_fi;
        
        %% ��Ч����ɢ�����
        sigma_eff = int(sigma_abs*g_L,nu,-Inf,Inf);
        sigma_eff = sigma_eff / (4*pi);
        sigma_eff = double(sigma_eff);
        Data_eff_0(i,j) = sigma_eff;
    end
end

%% F_+��Ч����ɢ��������ѭ��

Data_eff_R = [];

% �������ѭ��
for j = 1:length(V)    
    %% ���߲�����ֵ��
    
    % Na��6��ԾǨ���ߵ�����Ƶ��ƫ��ֵ(Hz)
    nu01 = (nu_L - 0.7328e9 + 0.6302e9) * ((c+V(j))/c);
    nu02 = (nu_L - 0.6962e9 + 0.6302e9) * ((c+V(j))/c);
    nu03 = (nu_L - 0.6302e9 + 0.6302e9) * ((c+V(j))/c);
    nu04 = (nu_L + 1.0333e9 + 0.6302e9) * ((c+V(j))/c);
    nu05 = (nu_L + 1.0552e9 + 0.6302e9) * ((c+V(j))/c);
    nu06 = (nu_L + 1.0919e9 + 0.6302e9) * ((c+V(j))/c);
    
    % Na�ĳ���ϸԾǨǿ������
    f1 = 1/32;
    f2 = 5/32;
    f3 = 14/32;
    f4 = 2/32;
    f5 = 5/32;
    f6 = 5/32;
    
    %% ��һ������
    syms nu;
    g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-(nu_L+585e6)).^2 ./ (2.*sigma_L.^2));
    
    % �¶�ѭ��
    for i = 1:length(T)
        %% �����������RMS���
        sigma_D_1 = nu01 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_2 = nu02 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_3 = nu03 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_4 = nu04 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_5 = nu05 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_6 = nu06 * sqrt(k_B*T(i) / (M*c^2));
        
        %% ��ֵ�������ս���
        sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_2 = (1 / (sqrt(2*pi)*sigma_D_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_3 = (1 / (sqrt(2*pi)*sigma_D_3)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_4 = (1 / (sqrt(2*pi)*sigma_D_4)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_5 = (1 / (sqrt(2*pi)*sigma_D_5)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_6 = (1 / (sqrt(2*pi)*sigma_D_6)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        
        %% ����ԾǨ�ߵĵ������ս���
        sigma_abs_1 = sigma_0_1 * exp(-(nu-nu01).^2 / (2*sigma_D_1^2));
        sigma_abs_2 = sigma_0_2 * exp(-(nu-nu02).^2 / (2*sigma_D_2^2));
        sigma_abs_3 = sigma_0_3 * exp(-(nu-nu03).^2 / (2*sigma_D_3^2));
        sigma_abs_4 = sigma_0_4 * exp(-(nu-nu04).^2 / (2*sigma_D_4^2));
        sigma_abs_5 = sigma_0_5 * exp(-(nu-nu05).^2 / (2*sigma_D_5^2));
        sigma_abs_6 = sigma_0_6 * exp(-(nu-nu06).^2 / (2*sigma_D_6^2));
        
        %% �������ս�����Գ���ϸԾǨǿ������
        sigma_abs_1_fi = sigma_abs_1*f1;
        sigma_abs_2_fi = sigma_abs_2*f2;
        sigma_abs_3_fi = sigma_abs_3*f3;
        sigma_abs_4_fi = sigma_abs_4*f4;
        sigma_abs_5_fi = sigma_abs_5*f5;
        sigma_abs_6_fi = sigma_abs_6*f6;
        
        %% ���ս������
        sigma_abs = sigma_abs_1_fi + sigma_abs_2_fi + sigma_abs_3_fi + sigma_abs_4_fi + sigma_abs_5_fi + sigma_abs_6_fi;
        
        %% ��Ч����ɢ�����
        sigma_eff = int(sigma_abs*g_L,nu,-Inf,Inf);
        sigma_eff = sigma_eff / (4*pi);
        sigma_eff = double(sigma_eff);
        Data_eff_R(i,j) = sigma_eff;
    end
end

%% F_-��Ч����ɢ��������ѭ��

Data_eff_L = [];

% �������ѭ��
for j = 1:length(V)    
    %% ���߲�����ֵ��
    
    % Na��6��ԾǨ���ߵ�����Ƶ��ƫ��ֵ(Hz)
    nu01 = (nu_L - 0.7328e9 + 0.6302e9) * ((c+V(j))/c);
    nu02 = (nu_L - 0.6962e9 + 0.6302e9) * ((c+V(j))/c);
    nu03 = (nu_L - 0.6302e9 + 0.6302e9) * ((c+V(j))/c);
    nu04 = (nu_L + 1.0333e9 + 0.6302e9) * ((c+V(j))/c);
    nu05 = (nu_L + 1.0552e9 + 0.6302e9) * ((c+V(j))/c);
    nu06 = (nu_L + 1.0919e9 + 0.6302e9) * ((c+V(j))/c);
    
    % Na�ĳ���ϸԾǨǿ������
    f1 = 1/32;
    f2 = 5/32;
    f3 = 14/32;
    f4 = 2/32;
    f5 = 5/32;
    f6 = 5/32;
    
    %% ��һ������
    syms nu;
    g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-(nu_L-585e6)).^2 ./ (2.*sigma_L.^2));
    
    % �¶�ѭ��
    for i = 1:length(T)
        %% �����������RMS���
        sigma_D_1 = nu01 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_2 = nu02 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_3 = nu03 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_4 = nu04 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_5 = nu05 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_6 = nu06 * sqrt(k_B*T(i) / (M*c^2));
        
        %% ��ֵ�������ս���
        sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_2 = (1 / (sqrt(2*pi)*sigma_D_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_3 = (1 / (sqrt(2*pi)*sigma_D_3)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_4 = (1 / (sqrt(2*pi)*sigma_D_4)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_5 = (1 / (sqrt(2*pi)*sigma_D_5)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_6 = (1 / (sqrt(2*pi)*sigma_D_6)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        
        %% ����ԾǨ�ߵĵ������ս���
        sigma_abs_1 = sigma_0_1 * exp(-(nu-nu01).^2 / (2*sigma_D_1^2));
        sigma_abs_2 = sigma_0_2 * exp(-(nu-nu02).^2 / (2*sigma_D_2^2));
        sigma_abs_3 = sigma_0_3 * exp(-(nu-nu03).^2 / (2*sigma_D_3^2));
        sigma_abs_4 = sigma_0_4 * exp(-(nu-nu04).^2 / (2*sigma_D_4^2));
        sigma_abs_5 = sigma_0_5 * exp(-(nu-nu05).^2 / (2*sigma_D_5^2));
        sigma_abs_6 = sigma_0_6 * exp(-(nu-nu06).^2 / (2*sigma_D_6^2));
        
        %% �������ս�����Գ���ϸԾǨǿ������
        sigma_abs_1_fi = sigma_abs_1*f1;
        sigma_abs_2_fi = sigma_abs_2*f2;
        sigma_abs_3_fi = sigma_abs_3*f3;
        sigma_abs_4_fi = sigma_abs_4*f4;
        sigma_abs_5_fi = sigma_abs_5*f5;
        sigma_abs_6_fi = sigma_abs_6*f6;
        
        %% ���ս������
        sigma_abs = sigma_abs_1_fi + sigma_abs_2_fi + sigma_abs_3_fi + sigma_abs_4_fi + sigma_abs_5_fi + sigma_abs_6_fi;
        
        %% ��Ч����ɢ�����
        sigma_eff = int(sigma_abs*g_L,nu,-Inf,Inf);
        sigma_eff = sigma_eff / (4*pi);
        sigma_eff = double(sigma_eff);
        Data_eff_L(i,j) = sigma_eff;
    end
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
f_T = fit_T;

% �����¶��������
figure('Name','�¶��������');
plot(f_T);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Fitting Surface of Temperature');
xlabel('R_{V}');
ylabel('R_{T}');
zlabel('Temperature (K)');

% �Է��ٽ��ж�ά�������
ft = fittype( 'poly55' );
[fit_V, gof_V] = fit( [R_V(:), R_T(:)], V(:), ft );
coe_fit_v=coeffvalues(fit_V);
f_V = fit_V;

% ���������������
figure('Name','�����������');
plot(f_V);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Fitting Surface of Wind Velocity');
xlabel('R_{V}');
ylabel('R_{T}');
zlabel('Wind Velocity (m/s)');

% ��ȡ�Ʒ����״�ɼ����ݾ���
Na_table = readtable('G:\Mohe_Data\589Na\Na20230709150515310.dat');
Na_data = table2array(Na_table);
height_num = Na_data(:,1);
Na_photon = Na_data(:,2:10);
Na_photon = Na_photon - mean(Na_photon(4232:4558,:),1);
Na_photon = Na_photon ./ max(Na_photon);

% ��1kmΪ�ռ�ֱ��ʽ��й������кϲ�
Na_i_read = [];
for i = 1:256
    Na_i_read(i,:) = sum(Na_photon(1+32*(i-1):32*i,:),1);
end
Na_photon = Na_i_read;

% �������ʵ�ʸ߶�����
height = 1:256;
height = height * 0.9830390625;
height = height(1:221);

% �Բ�ͬ����ͬƵ�ʽ��в��
V_F_0 = Na_photon(:,1);
V_F_R = Na_photon(:,2);
V_F_L = Na_photon(:,3);
N_F_0 = (Na_photon(:,4))*(sqrt(3)/2);
N_F_R = (Na_photon(:,5))*(sqrt(3)/2);
N_F_L = (Na_photon(:,6))*(sqrt(3)/2);
W_F_0 = (Na_photon(:,7))*(sqrt(3)/2);
W_F_R = (Na_photon(:,8))*(sqrt(3)/2);
W_F_L = (Na_photon(:,9))*(sqrt(3)/2);

% ����RV��ֱ���򡢱��򡢶���
R_V_real_V = (V_F_R - V_F_L) ./ (V_F_0);  % Ǧ��RV
R_V_real_N = (N_F_R - N_F_L) ./ (N_F_0);  % ����RV
R_V_real_W = (W_F_R - W_F_L) ./ (W_F_0);  % ����RV

% ����RT����
R_T_real_V = (V_F_R + V_F_L) ./ (2 * V_F_0);
R_T_real_N = (N_F_R + N_F_L) ./ (2 * N_F_0);
R_T_real_W = (W_F_R + W_F_L) ./ (2 * W_F_0);

% ���㴹ֱ����
V_real_V = f_V(R_V_real_V,R_T_real_V);
V_real_V = V_real_V(1:221,:);

% ���㱱�����,������sin��30�㣩
V_real_N = zeros(221,1);
for i = 1:221
    V_real_N(i,:) = f_V(R_V_real_N(round(i*(2/sqrt(3))),:),R_T_real_N(round(i*(2/sqrt(3))),:));
end
V_real_N = V_real_N * 0.5;

% ���㶫�����,������sin��30�㣩
V_real_W = zeros(221,1);
for i = 1:221
    V_real_W(i,:) = f_V(R_V_real_W(round(i*(2/sqrt(3))),:),R_T_real_W(round(i*(2/sqrt(3))),:));
end
V_real_W = V_real_W * 0.5;

% �������
V_real_Array = [V_real_N,V_real_W,V_real_V];
V_real_length = zeros(221,1);
V_real_direction = zeros(221,3);
for i = 1:221
    V_real_length(i,:) = norm(V_real_Array(i,:));
    V_real_direction(i,:) = V_real_Array(i,:) / norm(V_real_Array(i,:));
end

% �����¶�
T_real = f_T(R_V_real_V,R_T_real_V);
T_real = T_real(1:221,:);

% ���¶�ͼ
yyyy = 2023;
dd = 181+9;
hh = 15;
[T_model,Den_model] = atmosnrlmsise00([height*1e3],53.5,122.3,yyyy,dd,hh);
T_model = T_model(:,2);
figure('Name','��ԭ���¶ȷ���ͼ');
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

% ������ͼ
figure('Name','��ԭ�ӷ��ٷ���ͼ');
plot(V_real_length,height,'-b','linewidth',2);
xlim([0 120]);
ylim([80 105]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Wind Velocity 20230709 23:05');
xlabel('Wind Velocity (m/s)');
ylabel('Altitude (km)');

% ���糡ͼ
figure('Name','��ԭ�ӷ糡����ͼ');
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


