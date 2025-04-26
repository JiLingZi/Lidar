% Ray�ܶ�
% ��ȡ����������
ReadDate = '20231106';
fldstart = ['G:\ZWDATA\MOHEnew\K532\',ReadDate,'\'];
folder = [fldstart,'CH2\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;% Ϊ��ͳһ�ļ�������һ
% ��һ����Ϊ�Ǹ��ƻ�����������һ���ļ������²������ĸ�ͨ���ļ�����ͬ
% ��ʽ��Ҫ�ٸ��ˣ�̽��߶�����Ҳ��Ҫ�ٸ��ˣ�����ʵ���Զ���������������

DenPh = zeros(4096, num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    fid = fopen(filename);
    PhData = textscan(fid, '%f %f','HeaderLines',20);   % 7��22�ռ�֮ǰ��Ϊ23����Ӧ��Ϊ20
    fclose(fid);
    DenPh(:,j) = PhData{1,2}(1:4096);
end

% ��ȡʱ������
TimeValues = cell(1,num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    fid = fopen(filename);
    TXT = textscan(fid, '%s', 'Delimiter', '\n');
    RowTime = TXT{1}{9};
    TimeData = strsplit(RowTime);
    fclose(fid);
    TimeValues(:,j) = cellstr(TimeData{2:end});
end

% �ϲ�ʱ���ļ���
jNum = 30;
TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

% ��ȡԭʼ�߶Ⱦ���
height_num_origin = PhData{1,1}(1:4096);

% �ϲ��У�3�кϲ���96��
PhiSum = zeros(floor(4096/3),num_files);
for i = 1:floor(4096/3)
    PhiSum(i,:) = sum(DenPh(3*(i-1)+2:3*i+1,:),1);
end

% �ϲ��У�15�кϲ���15����
PhjSum = zeros(floor(4096/3),floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    PhjSum(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
end

% ����120-150��������������߶�
Altitude = height_num_origin(1:floor(4096/3))*3;
KM30 = size(Altitude(Altitude<30),1)+1;
KM40 = size(Altitude(Altitude<40),1)+1;
KM70 = size(Altitude(Altitude<70),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM125 = size(Altitude(Altitude<125),1)+1;

% ��һ�����������
CheckH = 50;% ����������һ���߶ȣ�KM
KMCheck = size(Altitude(Altitude<CheckH),1);
CheckPh = PhjSum./PhjSum(KMCheck,:);
figure('Name','���������','position',[10 100 1350 500])
for j = 1:size(CheckPh,2)
    plot(CheckPh(:,j)+j,Altitude,'-','linewidth',1.5);
    hold on;
end
xlim([0 31]);
ylim([75 115]);
xInterval = 1;  % x��̶ȼ��



% ��������
R = 8.314;                  % ����ѧ���� J/(K*mol)
M = 28.959e-3;              % ����Ħ������ kg/mol
NA = 6.023e23;              % ����٤���޳���

% �������ٶ�
G_all = 6.67259e-11;        % ������������
M_Earth = 5.965e24;         % �������� kg
R_Earth = 6371004;          % ����뾶 m
G = (G_all*M_Earth)./((R_Earth+(Altitude.*1e3)).^2);

% ��ȡNRLMSIS����ģ���¶����ܶ�,Rho_Z_0ȡ40km��T_Z_0ȡ70km
Z_0_NRL = Altitude(KM40)*1e3;
Z_top = Altitude(KM70)*1e3;
[TemRay,DenRay] = atmosnrlmsise00([Z_0_NRL,Z_top],40.5,116.0,2023,306,17);
Rho_Z_0 = (sum(DenRay(1,:))./NA) * M;
T_Z_top = TemRay(2,2);

% �ܶȷ���
Ph = PhjSum;              % ȡ20:00 ������
Z = Altitude;        % �߶Ⱦ���
Z_0 = Z_0_NRL*1e-3;               % �ܶȲο��߶�
N_Z = Ph;          % 30-60 km������
N_B = mean(Ph(KM120:KM125,:));  % ������ƽ������ˮƽ
% N_STD = std(Ph(126:147,:));
% Ratio_STD = 1-(N_STD./N_B);.*Ratio_STD
N_Z_0 = Ph(KM40,:);           % �ܶȲο��߶ȴ�������
Rho_Z = (Rho_Z_0 .* ((Z.^2)./(Z_0^2)) .* ((N_Z-N_B)./(N_Z_0-N_B)));

% �ܶ�������
DenErr = Rho_Z.*(1./sqrt(N_Z));

TimeX = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
TimeY = 1:size(TimeX,2);
hours = hour(TimeX);
minutes = minute(TimeX);
Time = (((hours(1,:)+8)*3600+minutes(1,:)*60)+((((hours(end,:)+8)*3600+minutes(end,:)*60)-((hours(1,:)+8)*3600+minutes(1,:)*60)).*(TimeY./TimeY(end))))./3600;
AltitudeY = Altitude(KM40:KM70,:);
DenZ = Rho_Z(KM40:KM70,:);


figure(1);
[Map,Line] = contourf(Time,AltitudeY,DenZ,30);
set(Line,'LineColor','none')
set(gca,'FontSize',12,'FontName','Times New Roman');
title(fldstart);
% set(gca, 'ColorScale', 'log');
colormap(jet); 
colorbar;
ylabel(colorbar,'Density (cm^{-3})');
xlabel('Local Time');
ylabel('Altitude (km)');
ylim([35 75]);


% �ܶȴ����ļ���
folderName = ['Proccessed\RayDen\',ReadDate,'\'];  % ָ���ļ�������
mkdir(folderName);       % �����ļ���

% �����ܶ����
DenError = 1./sqrt(PhjSum);
DenError(DenError>1)=NaN;
DataAltitude = Altitude(KM30:KM80,:);
for j = 1:size(TimeList,2)
    % �ϲ��߶ȣ��ܶȣ����
    DataDensity = Rho_Z(KM30:KM80,j);
    for i = 1:size(Altitude,1)
        if Rho_Z(i,:)<0
            DenError(i,:)=NaN;
        end
    end
    DataDensity(DataDensity<0)=NaN;
    DataError = real(DenError(KM30:KM80,j)).*100;
    data = [DataAltitude,DataDensity,DataError];

    % �½�TXT�ļ�
    datetimeStr = TimeList{j};  % �ӵ�Ԫ���л�ȡ����ʱ���ַ���
    
    datetimeObj = datetime(datetimeStr, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');  % ���ַ���ת��Ϊdatetime����
    WriteStr = datestr(datetimeObj, 'yyyymmddHHMMSS');  % ��datetime����

    PsdName = [folderName,'OMOHE_MUCL01_DUSM_L2_STP_',WriteStr,'_V01.00.TXT'];
    fileID = fopen(PsdName, 'w');
    RecordNumber = size(data,1);
    % ��ȡ��ǰʱ��
    dt_now = datetime('now');
    dt_str_now = datestr(dt_now, 'yyyy-mm-ddTHH:MM:SS');
    % д���ͷ
    fprintf(fileID,'#DataName: Density of Upper Stratosphere and Mesosphere\n');
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
    fprintf(fileID,'#DeviceState: BeamDirect=0.0 degree\n');
    fprintf(fileID,'#DeviceSpec: WaveLen=532nm, PRF=30Hz, PlsEnergy=300mJ\n');
    fprintf(fileID,'#ObsParameters: PlsAccu=27000\n');
    fprintf(fileID,'#Quantities: Density of upper stratosphere and mesosphere (kg/m^3)\n');
    fprintf(fileID,'#Elev(km): Height, F7.3, missingdata=NaN\n');
    fprintf(fileID,'#DenRay(kg/m^3): Upper Stratospheric and Mesospheric Density, F8.5, missingdata=NaN\n');
    fprintf(fileID,'#DenEr(%%): Density Error, F5.1, missingdata=NaN\n');
    fprintf(fileID,'#---------------------------------\n');
    % д���������
    fprintf(fileID, '%+7s %+8s %+5s\n', 'Elev', 'DenRay', 'DenEr');
    % д�����ݱ��
    fprintf(fileID, '%7.3f %8.5f %5.1f\n', data.');
    % �ر��ļ�
    fclose(fileID);
end



% ����
Curve = Rho_Z.*G;
% plot(Altitude,Curve);


% �¶ȷ���
% Rho_Z_top = (sum(DenRay(2,:))./NA) * M;
Rho_Z_top = Rho_Z(KM70,:);
Tem = zeros((KM70-1),size(Rho_Z_top,2));
for j = 1:size(Rho_Z_top,2)
    T_Z = (1:(KM70-1))';
    for i = 1:(KM70-1)
        T_Z(i,:) = (Rho_Z_top(:,j)./Rho_Z(i,j)).*T_Z_top + ((M/R)./Rho_Z(i,j)).*((sum(Curve(i:KM70,j)))*96);
    end
    Tem(:,j) = T_Z;
end
TemZ = Tem(KM30:KM70-1,:);
AltitudeY2 = Altitude(KM30:KM70-1,:);
TemZ(TemZ>350) = NaN;

% ���¶�����ͼ
figure(2);
[Map,Line] = contourf(Time,AltitudeY2,TemZ,30);
set(Line,'LineColor','none')
set(gca,'FontSize',12,'FontName','Times New Roman');
title(fldstart);
% set(gca, 'ColorScale', 'log');
colormap(jet); 
colorbar;
ylabel(colorbar,'Temperature (K)');
xlabel('Local Time');
ylabel('Altitude (km)');
ylim([25 85]);
xlim([21 28]);


% �¶ȴ����ļ���
folderName = ['Proccessed\RayTem\',ReadDate,'\'];  % ָ���ļ�������
mkdir(folderName);       % �����ļ���

% �����¶����
DenError = Tem./sqrt(PhjSum(1:KM70-1,:));
DenError(DenError>20)=NaN;
DataAltitude = Altitude(KM30:KM70-1,:);
for j = 1:size(TimeList,2)
    % �ϲ��߶ȣ��ܶȣ����
    DataDensity = Tem(KM30:KM70-1,j);
    for i = 1:size(Tem,1)
        if Tem(i,j)>350
            DenError(i,j)=NaN;
        end
    end
    DataDensity(DataDensity>350)=NaN;
    DataError = real(DenError(KM30:KM70-1,j));
    data = [DataAltitude,DataDensity,DataError];

    % �½�TXT�ļ�
    datetimeStr = TimeList{j};  % �ӵ�Ԫ���л�ȡ����ʱ���ַ���
    
    datetimeObj = datetime(datetimeStr, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');  % ���ַ���ת��Ϊdatetime����
    WriteStr = datestr(datetimeObj, 'yyyymmddHHMMSS');  % ��datetime����

    PsdName = [folderName,'OMOHE_MUCL01_TUSM_L2_STP_',WriteStr,'_V01.00.TXT'];
    fileID = fopen(PsdName, 'w');
    RecordNumber = size(data,1);
    % ��ȡ��ǰʱ��
    dt_now = datetime('now');
    dt_str_now = datestr(dt_now, 'yyyy-mm-ddTHH:MM:SS');
    % д���ͷ
    fprintf(fileID,'#DataName: Temperature of Upper Stratosphere and Mesosphere\n');
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
    fprintf(fileID,'#DeviceState: BeamDirect=0.0 degree\n');
    fprintf(fileID,'#DeviceSpec: WaveLen=532nm, PRF=30Hz, PlsEnergy=300mJ\n');
    fprintf(fileID,'#ObsParameters: PlsAccu=54000\n');
    fprintf(fileID,'#Quantities: Temperature of upper stratosphere and mesosphere (K)\n');
    fprintf(fileID,'#Elev(km): Height, F7.3, missingdata=NaN\n');
    fprintf(fileID,'#TempRay(K): Upper Stratospheric and Mesospheric Temperatures, F7.1, missingdata=NaN\n');
    fprintf(fileID,'#TemEr(%%): Density Error, F5.1, missingdata=NaN\n');
    fprintf(fileID,'#---------------------------------\n');
    % д���������
    fprintf(fileID, '%+7s %+7s %+5s\n', 'Elev', 'TempRay', 'DenEr');
    % д�����ݱ��
    fprintf(fileID, '%7.3f %7.1f %5.1f\n', data.');
    % �ر��ļ�
    fclose(fileID);
end






% Tem10 = zeros(size(T_Z,1),1);
% for t = 6:(size(T_Z,1)-5)
%     Tem10(t,:) = mean(T_Z(t-5:t+5,:),1);
% end
% tint = [1 2 3 4 5 (size(T_Z,1)-5) (size(T_Z,1)-4) (size(T_Z,1)-3) (size(T_Z,1)-2) (size(T_Z,1)-1) (size(T_Z,1)-0)];
% Tem10(tint,:) =NaN;
% Tem = zeros(floor(size(Tem10,1)/10),1);
% 
% for t = 1:floor(size(Tem10,1)/10)
%     Tem(t,:) = Tem10(10*t-5,:);
% end
% 
% 
% plot(TemRho(:,2),Z,'-k','Linewidth',1);
% hold on;
% plot(Tem,Altitude2(1:(i2_70-1),:),'-m','Linewidth',1);
% hold on;
% ylim([30,60]);
% TemTitleStr = strcat('HaiKou ',' ',YYYY,MM,DD,' 01:00 Temperature');
% title(TemTitleStr);
% xlabel('Temperature (K)');
% ylabel('Altitude (km)');
% set(gca,'FontName','Times New Roman','FontSize',15);
% legend('Location', 'northwest');
% legend('Yanqing','NRLMSIS');












