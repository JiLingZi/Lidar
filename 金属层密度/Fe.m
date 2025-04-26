% Fe�ܶ�
% ��ȡ����������
ReadDate = '20231019';
fldstart = ['E:\RawData\ZWDATA\MOHEnew\Fe\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;% Ϊ��ͳһ�ļ�������һ
% ��һ����Ϊ�Ǹ��ƻ�����������һ���ļ������²������ĸ�ͨ���ļ�����ͬ
% ��ʽ��Ҫ�ٸ��ˣ�̽��߶�����Ҳ��Ҫ�ٸ��ˣ�����ʵ���Զ���������������


DenPh = zeros(8192, num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    fid = fopen(filename);
    PhData = textscan(fid, '%f %f','HeaderLines',20);   % 7��22�ռ�֮ǰ��Ϊ23����Ӧ��Ϊ20
    fclose(fid);
    DenPh(:,j) = PhData{1,2}(1:8192);
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
jNum = 1;
iNum = 1;
TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

% ��ȡԭʼ�߶Ⱦ���
height_num_origin = PhData{1,1}(1:8192);

% �ϲ��У�3�кϲ���96��
PhiSum = zeros(floor(8192/iNum),num_files);
for i = 1:floor(8192/iNum)
    PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
end

% �ϲ��У�15�кϲ���15����
PhjSum = zeros(floor(8192/iNum),floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    PhjSum(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
end

% ����120-150��������������߶�
Altitude = height_num_origin(1:floor(8192/iNum))*iNum;
KM30 = size(Altitude(Altitude<30),1)+1;
KM75 = size(Altitude(Altitude<75),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
KM115 = size(Altitude(Altitude<115),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM125 = size(Altitude(Altitude<125),1)+1;
KM350 = size(Altitude(Altitude<350),1)+1;
KM390 = size(Altitude(Altitude<390),1)+1;

% % ��һ�����������
% CheckH = 55;% ����������һ���߶ȣ�KM
% KMCheck = size(Altitude(Altitude<CheckH),1);
% CheckPh = PhjSum./PhjSum(KMCheck,:);
% figure('Name','���������','position',[10 100 1350 500])
% for j = 1:size(CheckPh,2)
%     plot(CheckPh(:,j)+j,Altitude,'-k','linewidth',1.5);
%     hold on;
% end
% set(gca,'FontSize',12,'FontName','Times New Roman');
% title('Fe Density Profiles 20231102 @Mohe');
% xlabel('Local Time');
% ylabel('Altitude (km)');
% xlim([-1 41]);
% ylim([75 115]);
% xInterval = 1;  % x��̶ȼ��

% ����
Noise = mean(PhjSum(KM120:KM125,:),1);
SumPh = PhjSum - Noise;

SNR = SumPh./Noise;



% ��ȡ��������ɢ��������Ч����ɢ�����
ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));
RayScatter = Scatter(2,1);
EffScatter = Scatter(2,2);

% ��ȡ30 km ������ģ���¶����ܶ�
[TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,2023,201,17);

% �����ܶ�
Z = Altitude;                       % �߶Ⱦ���
ZR = 30;                            % �ο��߶�
SigmaRay = RayScatter;              % ��������ɢ�����
SigmaEff = EffScatter;              % ��Ч����ɢ�����
NZ = PhjSum;                        % ����������
NB = Noise;                         % ��������
NZR = PhjSum(KM30,:);                % �ο��߶ȴ�������
NumberRay = sum(DenRay)*1e-6;       % �ο��߶ȴ�����ģ�����ܶ�
NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;



% �ܶ����ߴּ��
% NaNJ = [1:2];
% NumberZ(:,NaNJ) = NaN;

TS ='OKfE'

%% ������
PSJ = 571
Density = NumberZ(KM70:KM125,PSJ);
Density(Density<0) = 0;
AltitudeY = Altitude(KM70:KM125,:);
plot(AltitudeY,Density,'-b', 'linewidth', 1.5);
xlim([70 125]); 
% ylim([0 20])
xlabel('Altitude (km)'); ylabel('Density (cm^{��3})')
legend('Fe'); grid on;
ttstr = ['(c) ',char(string(TimeList(PSJ)))];
title(ttstr)
set(gca,'fontsize',15);

%%
figure(1)

tFIND = TimeList';

Density = NumberZ(KM75:KM125,28);
Density(Density<0) = 0;
TimeX = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
TimeY = 1:size(TimeX,2);
hours = hour(TimeX);
minutes = minute(TimeX);
Time = (((hours(1,:)+8)*3600+minutes(1,:)*60)+((((hours(end,:)+8)*3600+minutes(end,:)*60)-((hours(1,:)+8)*3600+minutes(1,:)*60)).*(TimeY./TimeY(end))))./3600;
AltitudeY = Altitude(KM75:KM125,:);
AltitudeY = Altitude(KM75:KM125,:);
plot(Density, AltitudeY, '-k', 'linewidth', 1.5);
title('Fe Density (cm^{-3}) 20230715 22:09 Mohe', 'FontWeight', 'bold');
% xlim([75 125]);
ylim([75 115])
xlabel('Density (cm^{-3})', 'FontWeight', 'bold')
ylabel('Altitude (km)', 'FontWeight', 'bold')


%% ����������
figure(2)

tFIND = TimeList';

Density = SumPh(KM75:KM115,28);
Density(Density<0) = 0;
TimeX = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
TimeY = 1:size(TimeX,2);
hours = hour(TimeX);
minutes = minute(TimeX);
Time = (((hours(1,:)+8)*3600+minutes(1,:)*60)+((((hours(end,:)+8)*3600+minutes(end,:)*60)-((hours(1,:)+8)*3600+minutes(1,:)*60)).*(TimeY./TimeY(end))))./3600;
AltitudeY = Altitude(KM75:KM115,:);
AltitudeY = Altitude(KM75:KM115,:);
plot(Density, AltitudeY, '-', 'linewidth', 1.5);
grid on;
title('20240526 13:56:47', 'FontWeight', 'bold');
% xlim([75 125]);
ylim([75 115])
xlabel('Photon Counts', 'FontWeight', 'bold')
ylabel('Altitude (km)', 'FontWeight', 'bold')
%% ���ܶ�ʱ��ͼ
window = fspecial('average',[5,2]);
Density(:,1) = 0;
DenZ = imfilter(Density,window,'symmetric','same');
figure(2)
TimeA = datenum(TimeX);
[Map,Line] = contourf(Time,AltitudeY,DenZ,200);
% datetick('x', 'HH:MM')
set(Line,'LineColor','none')
set(gca,'FontSize',12,'FontName','arial');
% title(fldstart);
title('(f) 20231103 Mohe');
% set(gca, 'ColorScale', 'log');
colormap(jet); 
colorbar;
ylabel(colorbar,'Fe density (cm^{-3})');
xlabel('Time (LT)');
ylabel('Altitude (km)');
ylim([75 115]);
xlim([19 25])
hold on
plot([23 23],[75 115],'--w','linewidth',1)

%% % ����ƽ������ͼ
% DenMean = mean(DenZ,2);
% figure('Name','MeanDenDay');
% plot(DenMean,AltitudeY,'-b','LINEWIDTH',2);
% hold on;
% [maxV,maxDex] = max(DenMean);
% upDenMean = DenMean(DenMean>(0.5*maxV));
% upDenMeanLong = size(upDenMean,1);
% [maxV2,maxDex2] = max(upDenMean);
% RMSL = maxDex-maxDex2;
% RMSR = maxDex-maxDex2+upDenMeanLong;
% RMS = AltitudeY(RMSR)-AltitudeY(RMSL)
% % plot([0.5*maxV,0.5*maxV],[AltitudeY(RMSL),AltitudeY(RMSR)],'-r','linewidth',2)
% set(gca,'FontSize',12,'FontName','Times New Roman');
% titlestr2=['Fe Density Mean Profile ',ReadDate,' @Yanqing'];
% title(titlestr2);
% xlabel('Density (cm^{-3})');
% ylabel('Altitude (km)');


% �����ļ���
folderName = ['Proccessed\Fe\',ReadDate,'\'];  % ָ���ļ�������
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
    
    datetimeObj = datetime(datetimeStr, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');  % ���ַ���ת��Ϊdatetime����
    WriteStr = datestr(datetimeObj, 'yyyymmddHHMMSS');  % ��datetime����

    PsdName = [folderName,'OMOHE_MUCL01_NDFE_L2_STP_',WriteStr,'_V01.00.TXT'];
    fileID = fopen(PsdName, 'w');
    RecordNumber = size(data,1);
    % ��ȡ��ǰʱ��
    dt_now = datetime('now');
    dt_str_now = datestr(dt_now, 'yyyy-mm-ddTHH:MM:SS');
    % д���ͷ
    fprintf(fileID,'#DataName: Number Density of Fe Atoms\n');
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
    fprintf(fileID,'#DeviceSpec: WaveLen=372nm, PRF=30Hz, PlsEnergy=10mJ\n');
    fprintf(fileID,'#ObsParameters: PlsAccu=27000\n');
    fprintf(fileID,'#Quantities: Number density of Fe atoms (cm^-3)\n');
    fprintf(fileID,'#Elev(km): Height, F7.3, missingdata=NaN\n');
    fprintf(fileID,'#FeDen(cm^-3): Atomic Number Density of Iron, F9.1, missingdata=NaN\n');
    fprintf(fileID,'#DenEr(%%): Density Error, F5.1, missingdata=NaN\n');
    fprintf(fileID,'#---------------------------------\n');
    % д���������
    fprintf(fileID, '%+7s %+9s %+5s\n', 'Elev', 'FeDen', 'DenEr');
    % д�����ݱ��
    fprintf(fileID, '%7.3f %9.1f %5.1f\n', data.');
    % �ر��ļ�
    fclose(fileID);
end
