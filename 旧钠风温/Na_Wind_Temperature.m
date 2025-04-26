%-------------------------------------------------------------------------
%=================�и߲�����Ʒ��¼����״���·��ݳ���======================
%-------------------------------------------------------------------------
% ��Ŀ���ƣ�Į����ԭ�ӷ��·��ݳ���
% �������ƣ�������
% ���ߵ�λ�����տƼ���ѧ
% ���ߵ�ʦ��������
% �������䣺wang.kexin@stu.just.edu.cn
% ����ʱ�䣺2023-07-09 17:30(UTC+8)
% �޸�ʱ�䣺2023-07-10 17:30(UTC+8)
% ������Դ�����繤�̶���Į�Ӽ����״�̨վ 
%===============================ע������==================================
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

% ������Ƶ��Ч����ɢ��������
load Data_eff_0;
load Data_eff_R;
load Data_eff_L;

% �¶�����
T = [80 100 120 140 160 180 200 220 240 260 280 300 320]';
T = repmat(T,1,size(T,1));

% ��������
V = [-120 -100 -80 -60 -40 -20 0 20 40 60 80 100 120];
V = repmat(V,size(V,2),1);

% ���±��ʾ���
R_T = (Data_eff_R + Data_eff_L) ./ (2 * Data_eff_0);
R_V = (Data_eff_R - Data_eff_L) ./ (Data_eff_0);

% ���¶Ƚ��ж�ά�������
ft = fittype( 'poly55' );
[fit_T, gof_T] = fit( [R_V(:), R_T(:)], T(:), ft );
coe_fit_t=coeffvalues(fit_T);
f_T = fit_T;

% �����¶��������
figure(1);
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
figure(2);
plot(f_V);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Fitting Surface of Wind Velocity');
xlabel('R_{V}');
ylabel('R_{T}');
zlabel('Wind Velocity (m/s)');

%%
% ��ȡ�Ʒ����״�ɼ����ݾ���
% Na_table = readtable('D:\ZWDATA\Na20230707\20230707\Na\Na20230707171542258.dat');
% Na_table = readtable('D:\ZWDATA\20230709Na\Na\Na20230709060144901.dat');
Na_table = readtable('G:\Mohe_Data\589Na\Na20230709150515310.dat');%�״η��ݳ�������
Na_data = table2array(Na_table);
height_num = Na_data(:,1);
Na_photon = Na_data(:,2:10);
Na_photon = Na_photon - mean(Na_photon(4232:4558,:),1);
Na_photon = Na_photon ./ max(Na_photon);

% ÿ1km�ϲ���������
Na_i_read = [];
for i = 1:256
    Na_i_read(i,:) = sum(Na_photon(1+32*(i-1):32*i,:),1);
end
Na_photon = Na_i_read;
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

% ����RT
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

% % �����뱱��߶��������ӣ�Ǧ���������ر�
% height = height * (sqrt(3)/2);

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

% ������ͼ,V_real_length������Ϸ��٣�ֻ����ֵ
Meridional_wind = V_real_N;     % �����
Zonal_wind = V_real_W;          % γ���
Vertical_wind = V_real_V;       %��ֱ��
figure('Name','��ԭ�ӷ��ٷ���ͼ');
plot(V_real_length,height,'-b','linewidth',2);
xlim([0 120]);
ylim([80 105]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Wind Velocity 20230709 23:05');
xlabel('Wind Velocity (m/s)');
ylabel('Altitude (km)');

% ���糡ͼ,53.49N,122.34E
figure('Name','��ԭ�ӷ糡����ͼ');
X = [53.49];
X = repmat(X,221,1);
Y = [122.34];
Y = repmat(Y,221,1);
quiver3(X, Y, height', V_real_direction(:,1), V_real_direction(:,2), V_real_direction(:,3), 'AutoScale','off','MaxHeadSize', 1 , 'LineWidth', 2);
zlim([85 100]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Wind Direction 20230709 23:05');
xlabel('Latitude (N)');
ylabel('Longitude (E)');
zlabel('Altitude (km)');








