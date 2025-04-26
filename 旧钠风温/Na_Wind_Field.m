%-------------------------------------------------------------------------
%======================�и߲�����Ʒ糡���ݳ���============================
%-------------------------------------------------------------------------
% ��Ŀ���ƣ�Į����ԭ�ӷ糡���ݳ���
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

%% ��Ч����ɢ��������ѭ��

Data_eff_0 = zeros(size(T,2),size(V,2));
Data_eff_R = zeros(size(T,2),size(V,2));
Data_eff_L = zeros(size(T,2),size(V,2));

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
    nu = nu_L-3e9:1e6:nu_L+3e9;
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
        sigma_eff_fx = conv(sigma_abs,g_L);
        sigma_eff_fx = sigma_eff_fx * 1e6;
        Data_eff_0(i,j) = sigma_eff_fx(6001)/(4*pi);
        Data_eff_R(i,j) = sigma_eff_fx(6001+585)/(4*pi);
        Data_eff_L(i,j) = sigma_eff_fx(6001-585)/(4*pi);
    end
end

% ������Ϸ��ٺ��¶�����
T = T';
T = repmat(T,1,size(V,2));
V = repmat(V,size(T,1),1);

% ���±��ʾ���
R_T = (Data_eff_R + Data_eff_L) ./ (2 * Data_eff_0);
R_V = (Data_eff_R - Data_eff_L) ./ (Data_eff_0);

% ���¶Ƚ��ж�ά�������
ft = fittype( 'poly55' );
[fit_T, gof_T] = fit( [R_V(:), R_T(:)], T(:), ft );
coe_fit_t=coeffvalues(fit_T);
f_T = fit_T;

% �Է��ٽ��ж�ά�������
ft = fittype( 'poly55' );
[fit_V, gof_V] = fit( [R_V(:), R_T(:)], V(:), ft );
coe_fit_v=coeffvalues(fit_V);
f_V = fit_V;

% ��ȡ�Ʒ����״�ɼ����ݾ���
folder = 'D:\ZWDATA\Mohe_20230625_0712\Field\';
% folder = 'D:\ZWDATA\Mohe_20230625_0712\Na07112300_0030\';
files = dir(fullfile(folder, '*.dat'));
num_files = length(files);  
V_F_0 = zeros(8192, num_files);
V_F_R = zeros(8192, num_files);
V_F_L = zeros(8192, num_files);
N_F_0 = zeros(8192, num_files);
N_F_R = zeros(8192, num_files);
N_F_L = zeros(8192, num_files);
W_F_0 = zeros(8192, num_files);
W_F_R = zeros(8192, num_files);
W_F_L = zeros(8192, num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    Na_table = readtable(filename);
    Na_data = table2array(Na_table);
    V_F_0(:,j) = Na_data(:,2);
    V_F_R(:,j) = Na_data(:,3);
    V_F_L(:,j) = Na_data(:,4);
    N_F_0(:,j) = Na_data(:,5);
    N_F_R(:,j) = Na_data(:,6);
    N_F_L(:,j) = Na_data(:,7);
    W_F_0(:,j) = Na_data(:,8);
    W_F_R(:,j) = Na_data(:,9);
    W_F_L(:,j) = Na_data(:,10);
end
height_num = Na_data(:,1);
H_FUCK = height_num*(sqrt(3)/2);

% �۳�130 to 150 km ����,������γ�����������cos(30)�����36km�����źŹ�һ��*(sqrt(3)/2))
V_F_0 = (V_F_0 - mean(V_F_0(4232:4558,:),1))./mean(V_F_0(977:1140,:),1);
V_F_R = (V_F_R - mean(V_F_R(4232:4558,:),1))./mean(V_F_R(977:1140,:),1);
V_F_L = (V_F_L - mean(V_F_L(4232:4558,:),1))./mean(V_F_L(977:1140,:),1);
N_F_0 = (N_F_0 - mean(N_F_0(4232:4558,:),1))./mean(N_F_0(977:1140,:),1);
N_F_R = (N_F_R - mean(N_F_R(4232:4558,:),1))./mean(N_F_R(977:1140,:),1);
N_F_L = (N_F_L - mean(N_F_L(4232:4558,:),1))./mean(N_F_L(977:1140,:),1);
W_F_0 = (W_F_0 - mean(W_F_0(4232:4558,:),1))./mean(W_F_0(977:1140,:),1);
W_F_R = (W_F_R - mean(W_F_R(4232:4558,:),1))./mean(W_F_R(977:1140,:),1);
W_F_L = (W_F_L - mean(W_F_L(4232:4558,:),1))./mean(W_F_L(977:1140,:),1);

% ÿ1km�ϲ���������
V_F_0_READ = zeros(256, size(V_F_0,2));
V_F_R_READ = zeros(256, size(V_F_0,2));
V_F_L_READ = zeros(256, size(V_F_0,2));
N_F_0_READ = zeros(256, size(V_F_0,2));
N_F_R_READ = zeros(256, size(V_F_0,2));
N_F_L_READ = zeros(256, size(V_F_0,2));
W_F_0_READ = zeros(256, size(V_F_0,2));
W_F_R_READ = zeros(256, size(V_F_0,2));
W_F_L_READ = zeros(256, size(V_F_0,2));
for i = 1:256
    V_F_0_READ(i,:) = sum(V_F_0(1+32*(i-1):32*i,:),1);
    V_F_R_READ(i,:) = sum(V_F_R(1+32*(i-1):32*i,:),1);
    V_F_L_READ(i,:) = sum(V_F_L(1+32*(i-1):32*i,:),1);
    N_F_0_READ(i,:) = sum(N_F_0(1+32*(i-1):32*i,:),1);
    N_F_R_READ(i,:) = sum(N_F_R(1+32*(i-1):32*i,:),1);
    N_F_L_READ(i,:) = sum(N_F_L(1+32*(i-1):32*i,:),1);
    W_F_0_READ(i,:) = sum(W_F_0(1+32*(i-1):32*i,:),1);
    W_F_R_READ(i,:) = sum(W_F_R(1+32*(i-1):32*i,:),1);
    W_F_L_READ(i,:) = sum(W_F_L(1+32*(i-1):32*i,:),1);    
end
V_F_0 = V_F_0_READ;
V_F_R = V_F_R_READ;
V_F_L = V_F_L_READ;
N_F_0 = N_F_0_READ;
N_F_R = N_F_R_READ;
N_F_L = N_F_L_READ;
W_F_0 = W_F_0_READ;
W_F_R = W_F_R_READ;
W_F_L = W_F_L_READ;
height = 1:256;
height = height * 0.9830390625;
height = height(1:221);

% ����RV��ֱ���򡢱��򡢶���
R_V_real_V = (V_F_R - V_F_L) ./ (V_F_0);        % Ǧ��RV
R_V_real_N = (N_F_R - N_F_L) ./ (N_F_0);        % ����(����)RV
R_V_real_W = (W_F_R - W_F_L) ./ (W_F_0);        % ����(γ��)RV

% ����RT
R_T_real_V = (V_F_R + V_F_L) ./ (2 * V_F_0);    % Ǧ��RT
R_T_real_N = (N_F_R + N_F_L) ./ (2 * N_F_0);    % ����(����)RT
R_T_real_W = (W_F_R + W_F_L) ./ (2 * W_F_0);    % ����(γ��)RT

% ���㴹ֱ����
V_real_V = f_V(R_V_real_V,R_T_real_V);
V_real_V = V_real_V(1:221,:);

% ����ֱ��ʱ��ͼ
figure('Name','��ԭ�Ӵ�ֱ��');
X = 1:num_files;
Y = height';
window = fspecial('average', [2,5]);
V_V_Window = imfilter(V_real_V, window, 'symmetric', 'same');
[Map, Line]=contourf(X,Y(87:102,:),V_V_Window(87:102,:),7);
set(Line,'LineColor','k');
clabel(Map, Line);
% xlim([120 320]);
ylim([85 100]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Vertical Wind 20230709/10');
xlabel('Local Time');
ylabel('Altitude (km)');
colorbar;
ylabel(colorbar,'Wind Velocity (m/s)');
colormap(jet);

% ���㱱�򣨾���磩����,������sin��30�㣩
V_real_N = zeros(221,size(V_F_0,2));
for i = 1:221
    V_real_N(i,:) = f_V(R_V_real_N(round(i*(2/sqrt(3))),:),R_T_real_N(round(i*(2/sqrt(3))),:));
end
V_real_N = V_real_N * 0.5;

% �������ʱ��ͼ
figure('Name','��ԭ�Ӿ����');
X = 1:num_files;
Y = height';
window = fspecial('average', [1,5]);
V_N_Window = imfilter(V_real_N, window, 'symmetric', 'same');
[Map, Line]=contourf(X,Y(87:102,:),V_N_Window(87:102,:),7);
set(Line,'LineColor','k');
% clabel(Map, Line);
% xlim([120 320]);
ylim([85 100]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Meridional Wind 20230709/10');
xlabel('Local Time');
ylabel('Altitude (km)');
colorbar;
ylabel(colorbar,'Wind Velocity (m/s)');
colormap(jet);

% ���㶫��γ��磩����,������sin��30�㣩
V_real_W = zeros(221,size(V_F_0,2));
for i = 1:221
    V_real_W(i,:) = f_V(R_V_real_W(round(i*(2/sqrt(3))),:),R_T_real_W(round(i*(2/sqrt(3))),:));
end
V_real_W = V_real_W * 0.5;

% ��γ���ʱ��ͼ
figure('Name','��ԭ��γ���');
X = 1:num_files;
Y = height';
window = fspecial('average', [1,5]);
V_W_Window = imfilter(V_real_W, window, 'symmetric', 'same');
[Map, Line]=contourf(X,Y(87:102,:),V_W_Window(87:102,:),7);
set(Line,'LineColor','k');
% clabel(Map, Line);
% xlim([120 320]);
ylim([85 100]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Zonal Wind 20230709/10');
xlabel('Local Time');
ylabel('Altitude (km)');
colorbar;
ylabel(colorbar,'Wind Velocity (m/s)');
colormap(jet);

% �������ٴ�С�ͷ���0����
% V_real_length = zeros(221,1);
V_real_direction = zeros(221,3);

% ���糡ͼ,53.49N,122.34E[1:size(V_F_0,2)(:,j)]��һ��X��quiver3�ڵ�X���޸�
fig = figure('Name','��ԭ�ӷ糡����ͼ','position',[500 50 400 580]);
X = zeros(221,1);
Y = zeros(221,1);
V_Vertical = zeros(221,1);
Z = height';
for j = 1:size(V_F_0,2)
    V_real_Array = [V_real_N(:,j),V_real_W(:,j),V_real_V(:,j)];
    for i = 1:221
        V_real_direction(i,:) = V_real_Array(i,:) / norm(V_real_Array(i,:));
    end
    clf;
    % ���������Ҫ��ɸ��ĲŶ�
    quiver3(X, Y, Z, V_real_direction(:,2), V_real_direction(:,1), V_Vertical, 'AutoScale','off','MaxHeadSize', 0.3, 'LineWidth', 2);
    xlim([-3 3]);
    ylim([-3 2]);
    zlim([85 100]);
    view([-20 15]);
    set(gca,'FontSize',12,'FontName','Times New Roman');
    title('Na Wind Field 20230709/10');
    xlabel('Zonal');
    ylabel('Meridional');
    zlabel('Altitude (km)');
    hold on;
    pause(0.3);
    frames(j) = getframe(fig);% ¼�ƶ���
end

figure('Name','���Ŷ���','Position',[500 0 500 680]);
axis off;
% implay(frames,10); 
movie(frames);

% ����Ϊ��Ƶ
writerObj = VideoWriter('Na_Wind_Field.mp4');
writerObj.FrameRate = 5; % ����֡��Ϊÿ��10֡
open(writerObj);
writeVideo(writerObj, frames);
close(writerObj);

% ���糡��ͼ,53.49N,122.34E
figure('Name','��ԭ�ӷ糡����ͼ');
for j = 1:size(V_F_0,2)
    plot(height,0.1.*V_real_N(:,j)+j-1)
    hold on;
end
% 26(43) 28(44) 30(45) 32(46) 34(47) 36(48) 38(49) 40(50) 42(51) 44(52)
xlim([86 100]);
ylim([40 55]);
set(gca,'FontSize',12,'FontName','Times New Roman');
title('Na Wind Field 20230709/10');
xlabel('Local Time');
ylabel('53.49N,122.34E');
zlabel('Altitude (km)');

% �������糡
plot(height,V_real_N(:,45))
height_pis = height';

