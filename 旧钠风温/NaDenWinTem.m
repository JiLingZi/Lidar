% Na�ܶ�
% ��ȡ����������
ReadDate = '20231105';
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Na\',ReadDate,'\'];
folder = [fldstart,'Na\'];
files = dir(fullfile(folder, '*.dat'));
num_files = length(files)-1;% Ϊ��ͳһ�ļ�������һ

V_F_0 = zeros(8192, num_files);
V_F_R = zeros(8192, num_files);
V_F_L = zeros(8192, num_files);
N_F_0 = zeros(8192, num_files);
N_F_R = zeros(8192, num_files);
N_F_L = zeros(8192, num_files);
E_F_0 = zeros(8192, num_files);
E_F_R = zeros(8192, num_files);
E_F_L = zeros(8192, num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    Na_table = readtable(filename);
    Na_data = table2array(Na_table);
    V_F_0(:,j) = Na_data(1:8192,2);
    V_F_R(:,j) = Na_data(1:8192,3);
    V_F_L(:,j) = Na_data(1:8192,4);
    N_F_0(:,j) = Na_data(1:8192,5);
    N_F_R(:,j) = Na_data(1:8192,6);
    N_F_L(:,j) = Na_data(1:8192,7);
    E_F_0(:,j) = Na_data(1:8192,8);
    E_F_R(:,j) = Na_data(1:8192,9);
    E_F_L(:,j) = Na_data(1:8192,10);
end

DenPh = V_F_0 ;

% ��ȡʱ������
TimeValues = cell(1,num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    fid = fopen(filename);
    TXT = textscan(fid, '%s', 'Delimiter', '\n');
    RowTime = TXT{1}{4};
    TimeData = strsplit(RowTime);
    fclose(fid);
    TimeValues(:,j) = cellstr(TimeData{4});
end

RunRead = '���ݶ�ȡOK��'

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

%% �ϲ�ʱ���ļ���
jNum = 1;
iNum = 1;

TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
	TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

TimeX = datetime(TimeList, 'InputFormat', ':yyyy-MM-dd''T''HH:mm:ss.SSS');

% ��ȡԭʼ�߶Ⱦ���
height_num_origin = Na_data(1:8192,1);

% �ϲ��У�3�кϲ���96��
PhiSum = zeros(floor(8192/iNum),num_files);
VF_0i = zeros(floor(8192/iNum),num_files);
VF_Ri = zeros(floor(8192/iNum),num_files);
VF_Li = zeros(floor(8192/iNum),num_files);
NF_0i = zeros(floor(8192/iNum),num_files);
NF_Ri = zeros(floor(8192/iNum),num_files);
NF_Li = zeros(floor(8192/iNum),num_files);
EF_0i = zeros(floor(8192/iNum),num_files);
EF_Ri = zeros(floor(8192/iNum),num_files);
EF_Li = zeros(floor(8192/iNum),num_files);
for i = 1:floor(8192/iNum)
    PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
    VF_0i(i,:) = sum(V_F_0(iNum*(i-1)+1:iNum*i,:),1);
    VF_Ri(i,:) = sum(V_F_R(iNum*(i-1)+1:iNum*i,:),1);
    VF_Li(i,:) = sum(V_F_L(iNum*(i-1)+1:iNum*i,:),1);
    NF_0i(i,:) = sum(N_F_0(iNum*(i-1)+1:iNum*i,:),1);
    NF_Ri(i,:) = sum(N_F_R(iNum*(i-1)+1:iNum*i,:),1);
    NF_Li(i,:) = sum(N_F_L(iNum*(i-1)+1:iNum*i,:),1);
    EF_0i(i,:) = sum(E_F_0(iNum*(i-1)+1:iNum*i,:),1);
    EF_Ri(i,:) = sum(E_F_R(iNum*(i-1)+1:iNum*i,:),1);
    EF_Li(i,:) = sum(E_F_L(iNum*(i-1)+1:iNum*i,:),1);
end

%% �ϲ��У�15�кϲ���15����
PhjSum = zeros(floor(8192/iNum),floor(num_files/jNum));
VF_0 = zeros(floor(8192/iNum),floor(num_files/jNum));
VF_R = zeros(floor(8192/iNum),floor(num_files/jNum));
VF_L = zeros(floor(8192/iNum),floor(num_files/jNum));
NF_0 = zeros(floor(8192/iNum),floor(num_files/jNum));
NF_R = zeros(floor(8192/iNum),floor(num_files/jNum));
NF_L = zeros(floor(8192/iNum),floor(num_files/jNum));
EF_0 = zeros(floor(8192/iNum),floor(num_files/jNum));
EF_R = zeros(floor(8192/iNum),floor(num_files/jNum));
EF_L = zeros(floor(8192/iNum),floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    PhjSum(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    VF_0(:,j) = sum(VF_0i(:,jNum*(j-1)+1:jNum*j),2);
    VF_R(:,j) = sum(VF_Ri(:,jNum*(j-1)+1:jNum*j),2);
    VF_L(:,j) = sum(VF_Li(:,jNum*(j-1)+1:jNum*j),2);
    NF_0(:,j) = sum(NF_0i(:,jNum*(j-1)+1:jNum*j),2);
    NF_R(:,j) = sum(NF_Ri(:,jNum*(j-1)+1:jNum*j),2);
    NF_L(:,j) = sum(NF_Li(:,jNum*(j-1)+1:jNum*j),2);
    EF_0(:,j) = sum(EF_0i(:,jNum*(j-1)+1:jNum*j),2);
    EF_R(:,j) = sum(EF_Ri(:,jNum*(j-1)+1:jNum*j),2);
    EF_L(:,j) = sum(EF_Li(:,jNum*(j-1)+1:jNum*j),2);
end

%% ����120-150��������������߶�
Altitude = height_num_origin(1:floor(8192/iNum))*iNum;
KM30 = size(Altitude(Altitude<30),1)+1;
KM35 = size(Altitude(Altitude<35),1)+1;
KM75 = size(Altitude(Altitude<75),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM85 = size(Altitude(Altitude<85),1)+1;
KM98 = size(Altitude(Altitude<98),1)+1;
KM100 = size(Altitude(Altitude<100),1)+1;
KM105 = size(Altitude(Altitude<105),1)+1;
KM115 = size(Altitude(Altitude<115),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM125 = size(Altitude(Altitude<125),1)+1;

% ����
Noise = mean(PhjSum(KM120:KM125,:),1);
SumPh = PhjSum - Noise;
SNR = SumPh./Noise;

% ��ȡ��������ɢ��������Ч����ɢ�����
ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));
RayScatter = Scatter(5,1);
EffScatter = Scatter(5,2);

% ��ȡ30 km ������ģ���¶����ܶ�
[TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,2024,103,17);

%% �����ܶ�
Z = Altitude;                       % �߶Ⱦ���
ZR = 30;                            % �ο��߶�
SigmaRay = RayScatter;              % ��������ɢ�����
SigmaEff = EffScatter;              % ��Ч����ɢ�����
NZ = PhjSum;                        % ����������
NB = Noise;                         % ��������
NZR = PhjSum(KM30,:);                % �ο��߶ȴ�������
NumberRay = sum(DenRay)*1e-6;       % �ο��߶ȴ�����ģ�����ܶ�
NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

%% ���ܶ�ʱ��ͼ
% NumberZ(:,1:2) = NaN;

DenZ = NumberZ(KM75:KM115,:);
TimeN = datenum(TimeX);
AltitudeY = Altitude(KM75:KM115,:);
figure('Name','Na Density')
window = fspecial('average', [1,5]);
DenWin = imfilter(DenZ, window, 'symmetric', 'same');
[Map,Line] = contourf(TimeN,AltitudeY,DenWin,30);
datetick('x', 'HH:MM');
% tiknum = 1/24;
% xticks(TimeN(1):tiknum:TimeN(end));
set(Line,'LineColor','none')
set(gca,'FontSize',12,'FontName','Times New Roman');
TTSTR = ['Na Density (cm^{-3}) ' ReadDate ', Mohe'];
title(TTSTR);
%set(gca, 'ColorScale', 'log');
colormap(jet); 
colorbar;
xlabel('Universal Time');
ylabel('Altitude (km)');
ylim([80 105]);

%% ���ݷ糡

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
%����������
% �¶�����
T = [80 100 120 140 160 180 200 220 240 260 280 300 320];
% ��������
V = [-120 -100 -80 -60 -40 -20 0 20 40 60 80 100 120];
%��Ч����ɢ��������ѭ��
Data_eff_0 = zeros(size(T,2),size(V,2));
Data_eff_R = zeros(size(T,2),size(V,2));
Data_eff_L = zeros(size(T,2),size(V,2));
% �������ѭ��
for j = 1:length(V)    
    %���߲�����ֵ��
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
    %��һ������
    nu = nu_L-3e9:1e6:nu_L+3e9;
    g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));
    % �¶�ѭ��
    for i = 1:length(T)
        %�����������RMS���
        sigma_D_1 = nu01 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_2 = nu02 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_3 = nu03 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_4 = nu04 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_5 = nu05 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_6 = nu06 * sqrt(k_B*T(i) / (M*c^2));
        %��ֵ�������ս���
        sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_2 = (1 / (sqrt(2*pi)*sigma_D_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_3 = (1 / (sqrt(2*pi)*sigma_D_3)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_4 = (1 / (sqrt(2*pi)*sigma_D_4)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_5 = (1 / (sqrt(2*pi)*sigma_D_5)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_6 = (1 / (sqrt(2*pi)*sigma_D_6)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        %����ԾǨ�ߵĵ������ս���
        sigma_abs_1 = sigma_0_1 * exp(-(nu-nu01).^2 / (2*sigma_D_1^2));
        sigma_abs_2 = sigma_0_2 * exp(-(nu-nu02).^2 / (2*sigma_D_2^2));
        sigma_abs_3 = sigma_0_3 * exp(-(nu-nu03).^2 / (2*sigma_D_3^2));
        sigma_abs_4 = sigma_0_4 * exp(-(nu-nu04).^2 / (2*sigma_D_4^2));
        sigma_abs_5 = sigma_0_5 * exp(-(nu-nu05).^2 / (2*sigma_D_5^2));
        sigma_abs_6 = sigma_0_6 * exp(-(nu-nu06).^2 / (2*sigma_D_6^2));
        %�������ս�����Գ���ϸԾǨǿ������
        sigma_abs_1_fi = sigma_abs_1*f1;
        sigma_abs_2_fi = sigma_abs_2*f2;
        sigma_abs_3_fi = sigma_abs_3*f3;
        sigma_abs_4_fi = sigma_abs_4*f4;
        sigma_abs_5_fi = sigma_abs_5*f5;
        sigma_abs_6_fi = sigma_abs_6*f6;
        %���ս������
        sigma_abs = sigma_abs_1_fi + sigma_abs_2_fi + sigma_abs_3_fi + sigma_abs_4_fi + sigma_abs_5_fi + sigma_abs_6_fi;
        %��Ч����ɢ�����
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

%% ��ʼ�����¶�

VF_0 = VF_0./VF_0(KM35,:);
VF_R = VF_R./VF_R(KM35,:);
VF_L = VF_L./VF_L(KM35,:);
NF_0 = NF_0./NF_0(KM35,:);
NF_R = NF_R./NF_R(KM35,:);
NF_L = NF_L./NF_L(KM35,:);
EF_0 = EF_0./EF_0(KM35,:);
EF_R = EF_R./EF_R(KM35,:);
EF_L = EF_L./EF_L(KM35,:);


%% 

ChPf = 280;
tFIND = TimeList';

figure('name','V��Ƶ���','position',[250 50 350 800])
pf = ChPf;
plot(VF_0(:,pf), Altitude, '-k', 'linewidth',1.5)
hold on
plot(VF_R(:,pf), Altitude, '-b', 'linewidth',1.5)
hold on
plot(VF_L(:,pf), Altitude, '-r', 'linewidth',1.5)
grid on
% set(gca, 'XScale', 'log');
ylim([30 125])
% ttstr = string(pf) + ' | ' + TimeList(pf);
Tits = 'Vertical ';
ttstr = "Mohe Na " + TimeList(pf);
% title(ttstr)
title('(a) Vertical')
set(gca,'FontSize',12.5);
% ylim([3000 7000])
xlabel('Photon Counts')
ylabel('Altitude (km)')
legend('\nu_{0}','\nu_{+}','\nu_{-}')

% 
figure('name','N��Ƶ���','position',[600 50 350 800])
pf = ChPf;
plot(NF_0(:,pf), Altitude*(sqrt(3)/2), '-k', 'linewidth',1.5)
hold on
plot(NF_R(:,pf), Altitude*(sqrt(3)/2), '-b', 'linewidth',1.5)
hold on
plot(NF_L(:,pf), Altitude*(sqrt(3)/2), '-r', 'linewidth',1.5)
grid on
% set(gca, 'XScale', 'log');
ylim([30 125])
% ttstr = string(pf) + ' | ' + TimeList(pf);
Tits = 'Vertical ';
ttstr = "Mohe Na " + TimeList(pf);
title('(b) North')
set(gca,'FontSize',12.5);
% ylim([3000 7000])
xlabel('Photon Counts')
ylabel('Altitude (km)')
legend('\nu_{0}','\nu_{+}','\nu_{-}')

% 
figure('name','E��Ƶ���','position',[950 50 350 800])
pf = ChPf;
plot(EF_0(:,pf), Altitude*(sqrt(3)/2), '-k', 'linewidth',1.5)
hold on
plot(EF_R(:,pf), Altitude*(sqrt(3)/2), '-b', 'linewidth',1.5)
hold on
plot(EF_L(:,pf), Altitude*(sqrt(3)/2), '-r', 'linewidth',1.5)
grid on
% set(gca, 'XScale', 'log');
ylim([30 125])
% ttstr = string(pf) + ' | ' + TimeList(pf);
Tits = 'Vertical ';
ttstr = "Mohe Na " + TimeList(pf)
title('(c) East')
set(gca,'FontSize',12.5);
% ylim([3000 7000])
xlabel('Photon Counts')
ylabel('Altitude (km)')
legend('\nu_{0}','\nu_{+}','\nu_{-}')


%% 
figure('name','��Ƶ���','position',[300 300 800 250])
pf = 280;
plot(Altitude*(sqrt(3)/2),EF_0(:,pf), '-k', 'linewidth',1.5)
hold on
plot(Altitude*(sqrt(3)/2),EF_R(:,pf), '-b', 'linewidth',1.5)
hold on
plot(Altitude*(sqrt(3)/2),EF_L(:,pf), '-r', 'linewidth',1.5)
grid on
set(gca, 'YScale', 'log');
xlim([30 125])
% ttstr = string(pf) + ' | ' + TimeList(pf);
Tits = 'Na Photon Counts East';
ttstr = string(Tits) + TimeList(pf);
title(ttstr)
set(gca,'FontSize',12.5);
% ylim([3000 7000])
xlabel('Altitude (km)')
ylabel('Photon Counts')
legend('f_{0}','f_{+}','f_{-}')



%%

% ����RV
R_V_real_V = (VF_R - VF_L) ./ (VF_0);
R_V_real_N = (NF_R - NF_L) ./ (NF_0);
R_V_real_E = (EF_R - EF_L) ./ (EF_0);

% ����RT
R_T_real_V = (VF_R + VF_L) ./ (2 * VF_0);
R_T_real_N = (NF_R + NF_L) ./ (2 * NF_0);
R_T_real_E = (EF_R + EF_L) ./ (2 * EF_0);

% �����¶�
T_real = f_T(R_V_real_V,R_T_real_V);
T_real = T_real(KM85:KM100,:);

% ���¶�ͼ
% T_real(:,1:2) = NaN;

AltitudeT = Altitude(KM85:KM100,:);
figure('Name','��ԭ���¶ȷ���ͼ');
window = fspecial('average', [3,3]);
T_real = imfilter(T_real, window, 'symmetric', 'same');
[Map, Line]=contourf(TimeN,AltitudeT,T_real,10);
set(Line,'LineColor','k')
datetick('x', 'HH:MM');
% tiknum = 1/24;
% xticks(TimeN(1):tiknum:TimeN(end));
set(Line,'LineColor','none')
set(gca,'FontSize',12,'FontName','Times New Roman');
TTSTRt = ['Na Temperature (K) ' ReadDate ', Mohe'];
title(TTSTRt);
%set(gca, 'ColorScale', 'log');
colormap(jet); 
colorbar;
xlabel('Universal Time');
ylabel('Altitude (km)');
ylim([85 100]);

% ���㴹ֱ����
V_real_V = f_V(R_V_real_V,R_T_real_V);
V_real_V = V_real_V(KM85:KM100,:);

% ����ֱ��ʱ��ͼ
% V_real_V(:,1:2) = NaN;

figure('Name','��ԭ�Ӵ�ֱ�練��ͼ');
window = fspecial('average', [3,3]);
V_real_V = imfilter(V_real_V, window, 'symmetric', 'same');
[Map, Line]=contourf(TimeN,AltitudeT,V_real_V,10);
set(Line,'LineColor','k')
datetick('x', 'HH:MM');
% tiknum = 1/24;
% xticks(TimeN(1):tiknum:TimeN(end));
set(Line,'LineColor','none')
set(gca,'FontSize',12,'FontName','Times New Roman');
TTSTRVV = ['Na Wind Vertical (m/s) ' ReadDate ', Mohe'];
title(TTSTRVV);
%set(gca, 'ColorScale', 'log');
colormap(jet); 
colorbar;
xlabel('Universal Time');
ylabel('Altitude (km)');
ylim([85 100]);

% ����N����
V_real_N = f_V(R_V_real_N,R_T_real_N);
V_real_N = V_real_N(KM98:KM115,:);

% ����ֱ��ʱ��ͼ
% V_real_N(:,1:2) = NaN;

AltitudeNE = Altitude(KM98:KM115,:).*(sqrt(3)/2);
figure('Name','��ԭ��N�練��ͼ');
window = fspecial('average', [3,3]);
V_real_N = imfilter(V_real_N, window, 'symmetric', 'same');
[Map, Line]=contourf(TimeN,AltitudeNE,V_real_N,10);
set(Line,'LineColor','k')
datetick('x', 'HH:MM');
% tiknum = 1/24;
% xticks(TimeN(1):tiknum:TimeN(end));
set(Line,'LineColor','none')
set(gca,'FontSize',12,'FontName','Times New Roman');
TTSTRVV = ['Na Wind North (m/s) ' ReadDate ', Mohe'];
title(TTSTRVV);
%set(gca, 'ColorScale', 'log');
colormap(jet); 
colorbar;
xlabel('Universal Time');
ylabel('Altitude (km)');
ylim([85 100]);

% ����E����
V_real_E = f_V(R_V_real_E,R_T_real_E);
V_real_E = V_real_E(KM98:KM115,:);

% ����ֱ��ʱ��ͼ
% V_real_E(:,1:2) = NaN;

figure('Name','��ԭ��E�練��ͼ');
window = fspecial('average', [3,3]);
V_real_E = imfilter(V_real_E, window, 'symmetric', 'same');
[Map, Line]=contourf(TimeN,AltitudeNE,V_real_E,10);
set(Line,'LineColor','k')
datetick('x', 'HH:MM');
% tiknum = 1/24;
% xticks(TimeN(1):tiknum:TimeN(end));
set(Line,'LineColor','none')
set(gca,'FontSize',12,'FontName','Times New Roman');
TTSTRVV = ['Na Wind East (m/s) ' ReadDate ', Mohe'];
title(TTSTRVV);
%set(gca, 'ColorScale', 'log');
colormap(jet); 
colorbar;
xlabel('Universal Time');
ylabel('Altitude (km)');
ylim([85 100]);


