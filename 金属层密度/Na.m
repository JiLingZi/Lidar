% Na�ܶ�
% ��ȡ����������
ReadDate = '20231103';
fldstart = ['F:\RawData\ZWDATA\MoheNEW\Na\',ReadDate,'\'];
folder = [fldstart,'Na\'];
files = dir(fullfile(folder, '*.dat'));
num_files = length(files)-1;% Ϊ��ͳһ�ļ�������һ
% ��һ����Ϊ�Ǹ��ƻ�����������һ���ļ������²������ĸ�ͨ���ļ�����ͬ
% ��ʽ��Ҫ�ٸ��ˣ�̽��߶�����Ҳ��Ҫ�ٸ��ˣ�����ʵ���Զ���������������



DenPh = zeros(8192, num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    Na_table = readtable(filename);
    Na_data = table2array(Na_table);
    DenPh(:,j) = Na_data(:,2);
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
jNum = 1;
TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

% ��ȡԭʼ�߶Ⱦ���
height_num_origin = height_num(1:4096);

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
KM75 = size(Altitude(Altitude<75),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
KM115 = size(Altitude(Altitude<115),1)+1;
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

% ����
Noise = mean(PhjSum(KM120:KM125,:),1);
SumPh = PhjSum - Noise;

SNR = SumPh./Noise;

% ��ȡ��������ɢ��������Ч����ɢ����棬ɢ������ļ�Ӧ����ͬһ�ļ����£��������޸�·��
ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));
RayScatter = Scatter(5,1);
EffScatter = Scatter(5,2);

% ��ȡ30 km ������ģ���¶����ܶȣ����߶��ף�γ�ȣ����ȣ��꣬���еڼ��죬����ʱ��
[TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,2023,201,17);

% �����ܶ�
Z = Altitude;                       % �߶Ⱦ���
ZR = 30;                            % �ο��߶�
SigmaRay = RayScatter;              % ��������ɢ�����
SigmaEff = EffScatter;              % ��Ч����ɢ�����
NZ = PhjSum;                        % ����������
NB = Noise;                         % ��������
NZR = PhjSum(KM30,:);               % �ο��߶ȴ�������
NumberRay = sum(DenRay)*1e-6;       % �ο��߶ȴ�����ģ�����ܶ�
NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

% �ܶ����ߴּ�飬ȡ80-105km
% NaNJ = [1];
% NumberZ(:,NaNJ) = NaN; 
% 439

figure('Name','Check')
Density = NumberZ(KM75:KM125,:);
% Density(Density<0) = NaN;
Density(Density<0) = 0;
TimeX = datetime(TimeList, 'InputFormat', ':yyyy-MM-dd''T''HH:mm:ss.SSS');
TimeY = 1:size(TimeX,2);
hours = hour(TimeX);
minutes = minute(TimeX);
Time = (((hours(1,:)+8)*3600+minutes(1,:)*60)+((((hours(end,:)+8)*3600+minutes(end,:)*60)-((hours(1,:)+8)*3600+minutes(1,:)*60)).*(TimeY./TimeY(end))))./3600;
Time(1,1) = Time(1,1)-24;
AltitudeY = Altitude(KM75:KM125,:);
plot(Density,AltitudeY,'-k','linewidth',1.5);
title('(e) Na 20231103 23:00 Mohe')
xlabel('Density (cm^{-3})');
ylabel('Altitude (km)');
ylim([75 115]);



% ���ܶ�ʱ��ͼ
TimeXx = datenum(TimeX);

figure('Name','Varriation')
window = fspecial('average',[5,10]);
DenZ = imfilter(Density,window,'symmetric','same');
[Map,Line] = contourf(Time,AltitudeY,DenZ,200);
% datetick('x', 'HH:MM')
set(Line,'LineColor','none')
set(gca,'FontSize',12,'FontName','arial');
% title(fldstart);
title('(b) 20231103 Mohe');
% set(gca, 'ColorScale', 'log');
colormap(jet); 
colorbar;
ylabel(colorbar,'Na density (cm^{-3})');
xlabel('Time (LT)');
ylabel('Altitude (km)');
ylim([75 115]);
xlim([19 25])
hold on
plot([23 23],[75 115],'--w','linewidth',1)




AltPolor = AltitudeY';
TimeO = (TimeXx-TimeXx(1));
RangeO = 60;
TimePolor = (TimeO./TimeO(end))*RangeO -(0.5*RangeO);
DenP = DenZ;
window = fspecial('average',[1,1]);
DenPolor = imfilter(DenP,window,'symmetric','same');
polarPcolor( AltPolor, TimePolor, DenPolor);
colormap(jet); 




% �����ļ���
folderName = ['Proccessed\Na\',ReadDate,'\'];  % ָ���ļ�������
mkdir(folderName);       % �����ļ���

% �����ܶ����
DenError = 1./sqrt(PhjSum);
DenError(DenError>1)=NaN;
DataAltitude = Altitude(870:1195,:);
for j = 1:size(TimeList,2)
    % �ϲ��߶ȣ��ܶȣ����
    DataDensity = NumberZ(870:1195,j);
    for i = 1:size(Altitude,1)
        if NumberZ(i,:)<0
            DenError(i,:)=NaN;
        end
    end
    DataDensity(DataDensity<0)=NaN;
    DataError = real(DenError(870:1195,j)).*100;
    data = [DataAltitude,DataDensity,DataError];

    % �½�TXT�ļ�
    datetimeStr = TimeList{j};  % �ӵ�Ԫ���л�ȡ����ʱ���ַ���
    
    datetimeObj = datetime(datetimeStr, 'InputFormat', ':yyyy-MM-dd''T''HH:mm:ss.SSS');  % ���ַ���ת��Ϊdatetime����
    WriteStr = datestr(datetimeObj, 'yyyymmddHHMMSS');  % ��datetime����

    PsdName = [folderName,'OMOHE_MUCL01_NDNA_L2_STP_',WriteStr,'_V01.00.TXT'];
    fileID = fopen(PsdName, 'w');
    RecordNumber = size(data,1);
    % ��ȡ��ǰʱ��
    dt_now = datetime('now');
    dt_str_now = datestr(dt_now, 'yyyy-mm-ddTHH:MM:SS');
    % д���ͷ
    fprintf(fileID,'#DataName: Number Density of Na Atoms\n');
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
    fprintf(fileID,'#ObsParameters: PlsAccu=13500\n');
    fprintf(fileID,'#Quantities: Number density of Na atoms (cm^-3)\n');
    fprintf(fileID,'#Elev(km): Height, F7.3, missingdata=NaN\n');
    fprintf(fileID,'#NaDen(cm^-3): Atomic Number Density of Sodium, F9.1, missingdata=NaN\n');
    fprintf(fileID,'#DenEr(%%): Density Error, F5.1, missingdata=NaN\n');
    fprintf(fileID,'#---------------------------------\n');
    % д���������
    fprintf(fileID, '%+7s %+9s %+5s\n', 'Elev', 'NaDen', 'DenEr');
    % д�����ݱ��
    fprintf(fileID, '%7.3f %9.1f %5.1f\n', data.');
    % �ر��ļ�
    fclose(fileID);
end


% %polarPcolor(����,�Ƕ�,��Ҫ���Ƶı���,'Ncircles'��ԲȦ����,'Nspokes'�뾶����)
% polarPcolor(r,theta,ratial,'Ncircles',4,'Nspokes',13')
% %����ɫ����ɫ
% colormap jet;
% %����ɫ�귶Χ
% caxis([-10,10]);
% colorbar off %����ʾ�Դ�������
% c=colorbar;
% set(c,'YTick',-10:5:10); %ɫ��ֵ��Χ����ʾ���
% set(c,'YTickLabel',{'-10','-5','0','5','10'}) %����̶ȸ�ֵ
% ylabel(c,'radial wind speed (m/s)');
% hold on

