% Na�ܶ�

%% ��ȡ����������
ReadDate = '20231106';
fldstart = ['H:\ZWDATA\MoheNEW\Na\',ReadDate,'\'];
folder = [fldstart,'Na\'];
files = dir(fullfile(folder, '*.dat'));
num_files = length(files)-1;% Ϊ��ͳһ�ļ�������һ
% ��һ����Ϊ�Ǹ��ƻ�����������һ���ļ������²������ĸ�ͨ���ļ�����ͬ
% ��ʽ��Ҫ�ٸ��ˣ�̽��߶�����Ҳ��Ҫ�ٸ��ˣ�����ʵ���Զ���������������

VF0 = zeros(8192, num_files);
VFR = zeros(8192, num_files);
VFL = zeros(8192, num_files);
NF0 = zeros(8192, num_files);
NFR = zeros(8192, num_files);
NFL = zeros(8192, num_files);
EF0 = zeros(8192, num_files);
EFR = zeros(8192, num_files);
EFL = zeros(8192, num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    Na_table = readtable(filename);
    Na_data = table2array(Na_table);
    VF0(:,j) = Na_data(:,2);
    VFR(:,j) = Na_data(:,3);
    VFL(:,j) = Na_data(:,4);
    NF0(:,j) = Na_data(:,5);
    NFR(:,j) = Na_data(:,6);
    NFL(:,j) = Na_data(:,7);
    EF0(:,j) = Na_data(:,8);
    EFR(:,j) = Na_data(:,9);
    EFL(:,j) = Na_data(:,10);
end
height_num = Na_data(:,1);

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

% �ϲ�ʱ���ļ���
jNum = 30;
TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

% ��ȡԭʼ�߶Ⱦ���
height_num_origin = height_num(1:8192);

% �ϲ��У�32�кϲ���1000��
PhiV0 = zeros(floor(8192/32),num_files);
PhiVR = zeros(floor(8192/32),num_files);
PhiVL = zeros(floor(8192/32),num_files);
PhiN0 = zeros(floor(8192/32),num_files);
PhiNR = zeros(floor(8192/32),num_files);
PhiNL = zeros(floor(8192/32),num_files);
PhiE0 = zeros(floor(8192/32),num_files);
PhiER = zeros(floor(8192/32),num_files);
PhiEL = zeros(floor(8192/32),num_files);
for i = 1:floor(8192/32)
    PhiV0(i,:) = sum(VF0(32*(i-1)+1:32*i,:),1);
    PhiVR(i,:) = sum(VFR(32*(i-1)+1:32*i,:),1);
    PhiVL(i,:) = sum(VFL(32*(i-1)+1:32*i,:),1);
    PhiN0(i,:) = sum(NF0(32*(i-1)+1:32*i,:),1);
    PhiNR(i,:) = sum(NFR(32*(i-1)+1:32*i,:),1);
    PhiNL(i,:) = sum(NFL(32*(i-1)+1:32*i,:),1);
    PhiE0(i,:) = sum(EF0(32*(i-1)+1:32*i,:),1);
    PhiER(i,:) = sum(EFR(32*(i-1)+1:32*i,:),1);
    PhiEL(i,:) = sum(EFL(32*(i-1)+1:32*i,:),1);
end

% �ϲ��У�30�кϲ���60����
PhjV0 = zeros(floor(8192/32),floor(num_files/jNum));
PhjVR = zeros(floor(8192/32),floor(num_files/jNum));
PhjVL = zeros(floor(8192/32),floor(num_files/jNum));
PhjN0 = zeros(floor(8192/32),floor(num_files/jNum));
PhjNR = zeros(floor(8192/32),floor(num_files/jNum));
PhjNL = zeros(floor(8192/32),floor(num_files/jNum));
PhjE0 = zeros(floor(8192/32),floor(num_files/jNum));
PhjER = zeros(floor(8192/32),floor(num_files/jNum));
PhjEL = zeros(floor(8192/32),floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    PhjV0(:,j) = sum(PhiV0(:,jNum*(j-1)+1:jNum*j),2);
    PhjVR(:,j) = sum(PhiVR(:,jNum*(j-1)+1:jNum*j),2);
    PhjVL(:,j) = sum(PhiVL(:,jNum*(j-1)+1:jNum*j),2);
    PhjN0(:,j) = sum(PhiN0(:,jNum*(j-1)+1:jNum*j),2);
    PhjNR(:,j) = sum(PhiNR(:,jNum*(j-1)+1:jNum*j),2);
    PhjNL(:,j) = sum(PhiNL(:,jNum*(j-1)+1:jNum*j),2);
    PhjE0(:,j) = sum(PhiE0(:,jNum*(j-1)+1:jNum*j),2);
    PhjER(:,j) = sum(PhiER(:,jNum*(j-1)+1:jNum*j),2);
    PhjEL(:,j) = sum(PhiEL(:,jNum*(j-1)+1:jNum*j),2);
end

% ����120-150��������������߶�
Altitude = height_num_origin(1:floor(8192/32))*32;
KM30 = size(Altitude(Altitude<30),1)+1;
KM75 = size(Altitude(Altitude<75),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
KM115 = size(Altitude(Altitude<115),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM125 = size(Altitude(Altitude<125),1)+1;

Altitude30 = Altitude .* (sqrt(3)/2);
YKM30 = size(Altitude30(Altitude30<30),1)+1;
YKM75 = size(Altitude30(Altitude30<75),1)+1;
YKM80 = size(Altitude30(Altitude30<80),1)+1;
YKM110 = size(Altitude30(Altitude30<110),1)+1;
YKM115 = size(Altitude30(Altitude30<115),1)+1;
YKM120 = size(Altitude30(Altitude30<120),1)+1;
YKM125 = size(Altitude30(Altitude30<125),1)+1;



%% �����������
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

% �¶�����
T = [80 100 120 140 160 180 200 220 240 260 280 300 320];
% T = 100:1:300;
% ��������
V = [-120 -100 -80 -60 -40 -20 0 20 40 60 80 100 120];
% V = -100:1:100;

Data_eff_0 = zeros(size(T,2),size(V,2));
Data_eff_R = zeros(size(T,2),size(V,2));
Data_eff_L = zeros(size(T,2),size(V,2));

% �������ѭ��
for j = 1:length(V)    
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
    nu = nu_L-3e9:1e6:nu_L+3e9;
    g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));
    % �¶�ѭ��
    for i = 1:length(T)
        sigma_D_1 = nu01 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_2 = nu02 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_3 = nu03 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_4 = nu04 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_5 = nu05 * sqrt(k_B*T(i) / (M*c^2));
        sigma_D_6 = nu06 * sqrt(k_B*T(i) / (M*c^2));
        sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_2 = (1 / (sqrt(2*pi)*sigma_D_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_3 = (1 / (sqrt(2*pi)*sigma_D_3)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_4 = (1 / (sqrt(2*pi)*sigma_D_4)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_5 = (1 / (sqrt(2*pi)*sigma_D_5)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_0_6 = (1 / (sqrt(2*pi)*sigma_D_6)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
        sigma_abs_1 = sigma_0_1 * exp(-(nu-nu01).^2 / (2*sigma_D_1^2));
        sigma_abs_2 = sigma_0_2 * exp(-(nu-nu02).^2 / (2*sigma_D_2^2));
        sigma_abs_3 = sigma_0_3 * exp(-(nu-nu03).^2 / (2*sigma_D_3^2));
        sigma_abs_4 = sigma_0_4 * exp(-(nu-nu04).^2 / (2*sigma_D_4^2));
        sigma_abs_5 = sigma_0_5 * exp(-(nu-nu05).^2 / (2*sigma_D_5^2));
        sigma_abs_6 = sigma_0_6 * exp(-(nu-nu06).^2 / (2*sigma_D_6^2));
        sigma_abs_1_fi = sigma_abs_1*f1;
        sigma_abs_2_fi = sigma_abs_2*f2;
        sigma_abs_3_fi = sigma_abs_3*f3;
        sigma_abs_4_fi = sigma_abs_4*f4;
        sigma_abs_5_fi = sigma_abs_5*f5;
        sigma_abs_6_fi = sigma_abs_6*f6;
        sigma_abs = sigma_abs_1_fi + sigma_abs_2_fi + sigma_abs_3_fi + sigma_abs_4_fi + sigma_abs_5_fi + sigma_abs_6_fi;
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
% ��Ϸ��±��ʾ���
R_T = (Data_eff_R + Data_eff_L) ./ (2 * Data_eff_0);
R_V = (Data_eff_R - Data_eff_L) ./ (Data_eff_0);

% ���¶Ƚ��ж�ά�������
ft = fittype( 'poly55' );
[fit_T, gof_T] = fit( [R_V(:), R_T(:)], T(:), ft );
coe_fit_t=coeffvalues(fit_T);
f_T = fit_T;

% % �����¶��������
% figure('Name','�¶��������');
% plot(f_T);
% set(gca,'FontSize',12,'FontName','Times New Roman');
% title('Fitting Surface of Temperature');
% xlabel('R_{V}');
% ylabel('R_{T}');
% zlabel('Temperature (K)');

% �Է��ٽ��ж�ά�������
ft = fittype( 'poly55' );
[fit_V, gof_V] = fit( [R_V(:), R_T(:)], V(:), ft );
coe_fit_v=coeffvalues(fit_V);
f_V = fit_V;

% % ���������������
% figure('Name','�����������');
% plot(f_V);
% set(gca,'FontSize',12,'FontName','Times New Roman');
% title('Fitting Surface of Wind Velocity');
% xlabel('R_{V}');
% ylabel('R_{T}');
% zlabel('Wind Velocity (m/s)');

V_F_0 = (PhjV0 - mean(PhjV0(KM120:KM125,:),1))./PhjV0(KM30,:);
V_F_R = (PhjVR - mean(PhjVR(KM120:KM125,:),1))./PhjVR(KM30,:);
V_F_L = (PhjVL - mean(PhjVL(KM120:KM125,:),1))./PhjVL(KM30,:);
N_F_0 = (PhjN0 - mean(PhjN0(YKM120:YKM125,:),1))./PhjN0(YKM30,:);
N_F_R = (PhjNR - mean(PhjNR(YKM120:YKM125,:),1))./PhjNR(YKM30,:);
N_F_L = (PhjNL - mean(PhjNL(YKM120:YKM125,:),1))./PhjNL(YKM30,:);
E_F_0 = (PhjE0 - mean(PhjE0(YKM120:YKM125,:),1))./PhjE0(YKM30,:);
E_F_R = (PhjER - mean(PhjER(YKM120:YKM125,:),1))./PhjER(YKM30,:);
E_F_L = (PhjEL - mean(PhjEL(YKM120:YKM125,:),1))./PhjEL(YKM30,:);


% ����RV��ֱ���򡢱��򡢶���
R_V_real_V = (V_F_R - V_F_L) ./ (V_F_0);  % Ǧ��RV
R_V_real_N = (N_F_R - N_F_L) ./ (N_F_0);  % ����RV
R_V_real_E = (E_F_R - E_F_L) ./ (E_F_0);  % ����RV

% ����RT����
R_T_real_V = (V_F_R + V_F_L) ./ (2 * V_F_0);
R_T_real_N = (N_F_R + N_F_L) ./ (2 * N_F_0);
R_T_real_E = (E_F_R + E_F_L) ./ (2 * E_F_0);

% �������
V_real_V = f_V(R_V_real_V,R_T_real_V);
V_real_N = f_V(R_V_real_N,R_T_real_N) .* 2;
V_real_E = f_V(R_V_real_E,R_T_real_E) .* 2;

% �����¶�
T_real = f_T(R_V_real_V,R_T_real_V);


%% ���ȼ���
TemError = abs(real(T_real./sqrt(V_F_0+V_F_R+V_F_L))./100);
VNError = abs(real(V_real_N./sqrt(V_F_0+V_F_R+V_F_L))./10);
VEError = abs(real(V_real_E./sqrt(V_F_0+V_F_R+V_F_L))./50);


%% �¶ȴ����ļ���
folderName = ['Proccessed\NaTem\',ReadDate,'\'];  % ָ���ļ�������
mkdir(folderName);       % �����ļ���

% �����¶����
DenError = TemError;
DenError(DenError>20)=NaN;
DataAltitude = Altitude(KM80:KM110,:);
for j = 1:size(T_real,2)
    % �ϲ��߶ȣ��ܶȣ����
    DataDensity = T_real(KM80:KM110,j);
    for i = 1:size(T_real,1)
        if T_real(i,j)>300
            DenError(i,j)=NaN;
        end
        if T_real(i,j)<180
            DenError(i,j)=NaN;
        end
    end
    DataDensity(DataDensity>300)=NaN;
    DataDensity(DataDensity<180)=NaN;
    DataError = real(DenError(KM80:KM110,j));
    data = [DataAltitude,DataDensity,DataError];

    % �½�TXT�ļ�
    datetimeStr = TimeList{j};  % �ӵ�Ԫ���л�ȡ����ʱ���ַ���
    
    datetimeObj = datetime(datetimeStr, 'InputFormat', ':yyyy-MM-dd''T''HH:mm:ss.SSS');  % ���ַ���ת��Ϊdatetime����
    WriteStr = datestr(datetimeObj, 'yyyymmddHHMMSS');  % ��datetime����

    PsdName = [folderName,'OMOHE_MUCL01_TMSL_L2_STP_',WriteStr,'_V01.00.TXT'];
    fileID = fopen(PsdName, 'w');
    RecordNumber = size(data,1);
    % ��ȡ��ǰʱ��
    dt_now = datetime('now');
    dt_str_now = datestr(dt_now, 'yyyy-mm-ddTHH:MM:SS');
    % д���ͷ
    fprintf(fileID,'#DataName: Temperature of Sodium Layer\n');
    fprintf(fileID,'#CopyRight: Chinese Meridian Project\n');
    fprintf(fileID,'#Station: OMOHE(122.3E,53.6N,298m)\n');
    fprintf(fileID,'#Instrument: Middle-upper Atmosphere Wind-Temperature-Metal-Constituents Lidar\n');
    fprintf(fileID,'#Producer: National Space Science Center, CAS\n');
    fprintf(fileID,'#FileCreateTime: %s\n', dt_str_now);
    fprintf(fileID,'#DataLevel: L2\n');
    fprintf(fileID,'#DataVersion: V01.00\n');
    fprintf(fileID,'#DataTime: %s\n', datetimeStr);
    fprintf(fileID,'#RecordNumber: %d\n',RecordNumber);
    fprintf(fileID,'#QualityFlag: TBD\n');
    fprintf(fileID,'#DeviceState: BeamDirect=0.0 degree(Zenith), 30.0 degree(North), 30.0 degree(East)\n');
    fprintf(fileID,'#DeviceSpec: WaveLen=589nm, PRF=15Hz, PlsEnergy=10mJ\n');
    fprintf(fileID,'#ObsParameters: PlsAccu=27000\n');
    fprintf(fileID,'#Quantities: Temperature of Sodium Layer (K)\n');
    fprintf(fileID,'#Elev(km): Height, F7.3, missingdata=NaN\n');
    fprintf(fileID,'#TempNa(K): Sodium Layer Temperature, F6.1, missingdata=NaN\n');
    fprintf(fileID,'#TemEr(%%): Temperature Error, F5.1, missingdata=NaN\n');
    fprintf(fileID,'#---------------------------------\n');
    % д���������
    fprintf(fileID, '%+7s %+6s %+5s\n', 'Elev', 'TempNa', 'DenEr');
    % д�����ݱ��
    fprintf(fileID, '%7.3f %6.1f %5.1f\n', data.');
    % �ر��ļ�
    fclose(fileID);
end


%% �糡�����ļ���
folderName = ['Proccessed\NaWin\',ReadDate,'\'];  % ָ���ļ�������
mkdir(folderName);       % �����ļ���

% �����¶����
DenNError = VNError;
% DenNError(DenNError>20)=NaN;
DenEError = VEError;
% DenEError(DenEError>20)=NaN;
DataAltitude = Altitude(KM80:KM110,:);
for j = 1:size(T_real,2)
    % �ϲ��߶ȣ��ܶȣ����
    DataDensityN = V_real_N(KM80:KM110,j);
    for i = 1:size(V_real_N,1)
        if V_real_N(i,j)>150
            DenNError(i,j)=NaN;
        end
        if V_real_N(i,j)<-150
            DenNError(i,j)=NaN;
        end
    end
    DataDensityN(DataDensityN>150)=NaN;
    DataDensityN(DataDensityN<-150)=NaN;
    DataNError = real(DenNError(KM80:KM110,j));
    
    DataDensityE = V_real_E(KM80:KM110,j);
    for i = 1:size(V_real_E,1)
        if V_real_E(i,j)>150
            DenEError(i,j)=NaN;
        end
        if V_real_E(i,j)<-150
            DenEError(i,j)=NaN;
        end
    end
    DataDensityE(DataDensityE>150)=NaN;
    DataDensityE(DataDensityE<-150)=NaN;
    DataEError = real(DenEError(KM80:KM110,j));
    
    data = [DataAltitude,DataDensityE,DataDensityN,DataEError,DataNError];

    % �½�TXT�ļ�
    datetimeStr = TimeList{j};  % �ӵ�Ԫ���л�ȡ����ʱ���ַ���
    
    datetimeObj = datetime(datetimeStr, 'InputFormat', ':yyyy-MM-dd''T''HH:mm:ss.SSS');  % ���ַ���ת��Ϊdatetime����
    WriteStr = datestr(datetimeObj, 'yyyymmddHHMMSS');  % ��datetime����

    PsdName = [folderName,'OMOHE_MUCL01_HWSL_L2_STP_',WriteStr,'_V01.00.TXT'];
    fileID = fopen(PsdName, 'w');
    RecordNumber = size(data,1);
    % ��ȡ��ǰʱ��
    dt_now = datetime('now');
    dt_str_now = datestr(dt_now, 'yyyy-mm-ddTHH:MM:SS');
    % д���ͷ
    fprintf(fileID,'#DataName: Horizontal Wind of Sodium Layer\n');
    fprintf(fileID,'#CopyRight: Chinese Meridian Project\n');
    fprintf(fileID,'#Station: OMOHE(122.3E,53.6N,298m)\n');
    fprintf(fileID,'#Instrument: Middle-upper Atmosphere Wind-Temperature-Metal-Constituents Lidar\n');
    fprintf(fileID,'#Producer: National Space Science Center, CAS\n');
    fprintf(fileID,'#FileCreateTime: %s\n', dt_str_now);
    fprintf(fileID,'#DataLevel: L2\n');
    fprintf(fileID,'#DataVersion: V01.00\n');
    fprintf(fileID,'#DataTime: %s\n', datetimeStr);
    fprintf(fileID,'#RecordNumber: %d\n',RecordNumber);
    fprintf(fileID,'#QualityFlag: TBD\n');
    fprintf(fileID,'#DeviceState: BeamDirect=0.0 degree(Zenith), 30.0 degree(North), 30.0 degree(East)\n');
    fprintf(fileID,'#DeviceSpec: WaveLen=589nm, PRF=15Hz, PlsEnergy=10mJ\n');
    fprintf(fileID,'#ObsParameters: PlsAccu=27000\n');
    fprintf(fileID,'#Quantities: Horizontal wind of sodium layer (m/s)\n');
    fprintf(fileID,'#Elev(km): Height, F7.3, missingdata=NaN\n');
    fprintf(fileID,'#ZonalWind(m/s): Zonal Wind Speed in Sodium Layer, F9.1, missingdata=NaN\n');
    fprintf(fileID,'#MeridiWind(m/s): Meridional Wind Speed in Sodium Layer, F10.1, missingdata=NaN\n');
    fprintf(fileID,'#ZonalWEr(m/s): Zonal Wind Speed Error of Sodium Layer, F8.1, missingdata=NaN\n');
    fprintf(fileID,'#MeridiWEr(m/s): Meridional Wind Speed Error of Sodium Layer, F9.1, missingdata=NaN\n');
    fprintf(fileID,'#---------------------------------\n');
    % д���������
    fprintf(fileID, '%+7s %+9s %+10s %+8s %+9s\n', 'Elev', 'ZonalWind', 'MeridiWind' ,'ZonalWEr', 'MeridiWEr');
    % д�����ݱ��
    fprintf(fileID, '%7.3f %9.1f %10.1f %8.1f %9.1f\n', data.');
    % �ر��ļ�
    fclose(fileID);
end






% %% ���ȼ���
% N0 = V_F_0;
% varN0 = N0;
% NR = V_F_R;
% varNR = NR;
% NL = V_F_L;
% varNL= NL;
% 
% N0M = N_F_0;
% varN0M = N0M;
% NRM = N_F_R;
% varNRM = NRM;
% NLM = N_F_L;
% varNLM= NLM;
% 
% N0Z = E_F_0;
% varN0Z = N0Z;
% NRZ = E_F_R;
% varNRZ = NRZ;
% NLZ = E_F_L;
% varNLZ= NLZ;
% 
% % ����RT_real RV_real
% RT_real = (NR+NL)./(2*N0);
% RV_real = (NR-NL)./N0;
% 
% RTM_real = (NRM+NLM)./(2*N0M);
% RVM_real = (NRM-NL)./N0M;
% 
% RTZ_real = (NRZ+NLZ)./(2*N0Z);
% RVZ_real = (NRZ-NLZ)./N0Z;
% 
% % ����fT �� fV ���������ƫ��
% 
% % ���� RT RV ����
% syms RT RV;
% 
% % ���� fT ����
% t00 =       37.82;
% t10 =         330;
% t01 =       5.273;
% t20 =      -208.7;
% t11 =      -58.13;
% t02 =      -569.2;
% t30 =         570;
% t21 =         295;
% t12 =        2631;
% t03 =      -73.41;
% t40 =      -42.78;
% t31 =      -655.1;
% t22 =       -5104;
% t13 =       453.6;
% t04 =      -249.6;
% t50 =       66.43;
% t41 =        1090;
% t32 =        2825;
% t23 =      -799.6;
% t14 =       588.6;
% t05 =        6.18;
% fT(RT,RV) = t00+t10*RT+t01*RV+t20*RT^2+t11*RT*RV+t02*RV^2+t30*RT^3 ...
%          +t21*RT^2*RV+t12*RT*RV^2+t03*RV^3+t40*RT^4+t31*RT^3*RV ...
%          +t22*RT^2*RV^2+t13*RT*RV^3+t04*RV^4+t50*RT^5+t41*RT^4*RV ... 
%          +t32*RT^3*RV^2+t23*RT^2*RV^3+t14*RT*RV^4+t05*RV^5;
% 
% grad_fT = gradient(fT, [RT,RV]);
% dvalueTRTV = subs(grad_fT(RT_real,RV_real));
% dvalueTRT = double(dvalueTRTV(1:size(dvalueTRTV,1)/2,:));
% dvalueTRV = double(dvalueTRTV((size(dvalueTRTV,1)/2)+1:size(dvalueTRTV,1),:));
% 
% dvalueTRTVM = subs(grad_fT(RTM_real,RVM_real));
% dvalueTRTM = double(dvalueTRTVM(1:size(dvalueTRTVM,1)/2,:));
% dvalueTRVM = double(dvalueTRTVM((size(dvalueTRTVM,1)/2)+1:size(dvalueTRTVM,1),:));
% 
% dvalueTRTVZ = subs(grad_fT(RTZ_real,RVZ_real));
% dvalueTRTZ = double(dvalueTRTVZ(1:size(dvalueTRTVZ,1)/2,:));
% dvalueTRVZ = double(dvalueTRTVZ((size(dvalueTRTVZ,1)/2)+1:size(dvalueTRTVZ,1),:));
% 
% % ���� fV ����
% v00 =       11.66;
% v10 =       22.17;
% v01 =       538.2;
% v20 =      -129.4;
% v11 =       -2327;
% v02 =      -52.02;
% v30 =       266.4;
% v21 =        6336;
% v12 =       390.2;
% v03 =       185.4;
% v40 =      -311.4;
% v31 =       -7249;
% v22 =       -1126;
% v13 =      -970.4;
% v04 =        87.8;
% v50 =      -49.18;
% v41 =        3358;
% v32 =        1315;
% v23 =       679.7;
% v14 =      -235.3;
% v05 =        32.6;
% fV(RT,RV) = v00+v10*RT+v01*RV+v20*RT^2+v11*RT*RV+v02*RV^2+v30*RT^3 ...
%          +v21*RT^2*RV+v12*RT*RV^2+v03*RV^3+v40*RT^4+v31*RT^3*RV ...
%          +v22*RT^2*RV^2+v13*RT*RV^3+v04*RV^4+v50*RT^5+v41*RT^4*RV ... 
%          +v32*RT^3*RV^2+v23*RT^2*RV^3+v14*RT*RV^4+v05*RV^5;
% 
% grad_fV = gradient(fV, [RT,RV]);
% dvalueVRTV = subs(grad_fV(RT_real,RV_real));
% dvalueVRT = double(dvalueVRTV(1:size(dvalueVRTV,1)/2,:));
% dvalueVRV = double(dvalueVRTV((size(dvalueVRTV,1)/2)+1:size(dvalueVRTV,1),:));
% 
% dvalueVRTVM = subs(grad_fV(RTM_real,RVM_real));
% dvalueVRTM = double(dvalueVRTVM(1:size(dvalueVRTVM,1)/2,:));
% dvalueVRVM = double(dvalueVRTVM((size(dvalueVRTVM,1)/2)+1:size(dvalueVRTVM,1),:));
% 
% dvalueVRTVZ = subs(grad_fV(RTZ_real,RVZ_real));
% dvalueVRTZ = double(dvalueVRTVZ(1:size(dvalueVRTVZ,1)/2,:));
% dvalueVRVZ = double(dvalueVRTVZ((size(dvalueVRTVZ,1)/2)+1:size(dvalueVRTVZ,1),:));
% 
% % ���岢����V��ֱ��Pre_V
% Pre_V = (varNR.*(0.5.*dvalueVRT+dvalueVRV).^2+varNL.*(0.5.*dvalueVRT-dvalueVRV).^2+varN0.*(RT_real.*dvalueVRT+RV_real.*dvalueVRV).^2)./(N0.^2);
% Pre_V = double(Pre_V);
% 
% % ���岢����M�����Pre_V
% Pre_VM = (varNRM.*(0.5.*dvalueVRTM+dvalueVRVM).^2+varNLM.*(0.5.*dvalueVRTM-dvalueVRVM).^2+varN0M.*(RTM_real.*dvalueVRTM+RVM_real.*dvalueVRVM).^2)./(N0M.^2);
% Pre_VM = double(Pre_VM);
% 
% % ���岢����Zγ���Pre_V
% Pre_VZ = (varNRZ.*(0.5.*dvalueVRTZ+dvalueVRVZ).^2+varNLZ.*(0.5.*dvalueVRTZ-dvalueVRVZ).^2+varN0Z.*(RTZ_real.*dvalueVRTZ+RVZ_real.*dvalueVRVZ).^2)./(N0Z.^2);
% Pre_VZ = double(Pre_VZ);
% 
% % ���岢����Pre_T
% Pre_T = (varNR.*(0.5.*dvalueTRT+dvalueTRV).^2+varNL.*(0.5.*dvalueTRT-dvalueTRV).^2+varN0.*(RT_real.*dvalueTRT+RV_real.*dvalueTRV).^2)./(N0.^2);
% Pre_T = double(Pre_T);


