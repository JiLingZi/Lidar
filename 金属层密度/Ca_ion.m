% CaIon�ܶ�
% ��ȡ����������,��ԭ�ӡ���������4��ͨ��
ReadDate = '20231019';
fldstart = ['E:\RawData\ZWDATA\MOHEnew\Cap\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;% Ϊ��ͳһ�ļ�������һ
% ��һ����Ϊ�Ǹ��ƻ�����������һ���ļ������²������ĸ�ͨ���ļ�����ͬ
% ��ʽ��Ҫ�ٸ��ˣ�̽��߶�����Ҳ��Ҫ�ٸ��ˣ�����ʵ���Զ���������������


HL = 20;

DenPh1 = zeros(4096, num_files);
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
for j = 1:num_files
    filename = fullfile(folder,files(j).name);
    fid = fopen(filename);
    PhData1 = textscan(fid, '%f %f','HeaderLines',HL);
    fclose(fid);
    
    DenPh1(:,j) = PhData1{1,2}(1:4096,:);
end

DenPh2 = zeros(4096, num_files);
folder = [fldstart,'CH2\'];
files = dir(fullfile(folder, '*.txt'));
for j = 1:num_files
    filename = fullfile(folder,files(j).name);
    fid = fopen(filename);
    PhData2 = textscan(fid, '%f %f','HeaderLines',HL);
    fclose(fid);
    DenPh2(:,j) = PhData2{1,2}(1:4096,:);
end

DenPh3 = zeros(4096, num_files);
folder = [fldstart,'CH3\'];
files = dir(fullfile(folder, '*.txt'));
for j = 1:num_files
    filename = fullfile(folder,files(j).name);
    fid = fopen(filename);
    PhData3 = textscan(fid, '%f %f','HeaderLines',HL);
    fclose(fid);
    DenPh3(:,j) = PhData3{1,2}(1:4096,:);
end

DenPh4 = zeros(4096, num_files);
folder = [fldstart,'CH4\'];
files = dir(fullfile(folder, '*.txt'));
for j = 1:num_files
    filename = fullfile(folder,files(j).name);
    fid = fopen(filename);
    PhData4 = textscan(fid, '%f %f','HeaderLines',HL);
    fclose(fid);
    DenPh4(:,j) = PhData4{1,2}(1:4096,:);
end

DenPh = DenPh1+DenPh2+DenPh3+DenPh4;

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

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

disp(['���ݶ�ȡ���']);


%% �ϲ�ʱ���ļ���
jNum = 1;
iNum = 1;
TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

% ��ȡԭʼ�߶Ⱦ���
height_num_origin = PhData1{1,1};

% �ϲ��У�3�кϲ���96��
PhiSum = zeros(floor(4096/iNum),num_files);
for i = 1:floor(4096/iNum)
    PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
end

% �ϲ��У�15�кϲ���15����
PhjSum = zeros(floor(4096/iNum),floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    PhjSum(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
end

% ����120-150��������������߶�
Altitude = height_num_origin(1:floor(4096/iNum))*iNum;
KM30 = size(Altitude(Altitude<30),1)+1;
KM75 = size(Altitude(Altitude<75),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
KM115 = size(Altitude(Altitude<115),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM125 = size(Altitude(Altitude<125),1)+1;
KM180 = size(Altitude(Altitude<180),1)+1;
KM200 = size(Altitude(Altitude<200),1)+1;
KM350 = size(Altitude(Altitude<350),1)+1;
KM400 = size(Altitude(Altitude<400),1)+1;

% % ��һ�����������
% CheckH = 60;% ����������һ���߶ȣ�KM
% KMCheck = size(Altitude(Altitude<CheckH),1);
% CheckPh = PhjSum./PhjSum(KMCheck,:);
% figure('Name','���������','position',[10 100 1350 500])
% for j = 1:size(CheckPh,2)
%     plot(CheckPh(:,j)+j,Altitude,'-','linewidth',1.5);
%     hold on;
% end
% xlim([0 31]);
% ylim([75 115]);
% xInterval = 1;  % x��̶ȼ��

% ����
Noise = mean(PhjSum(KM120:KM125,:),1);
SumPh = PhjSum - Noise;

SNR = SumPh./Noise;



STD = std(PhjSum(KM120:KM125,:),1);
DetcRange = Noise + STD;



% ��ȡ��������ɢ��������Ч����ɢ�����
ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));
RayScatter = Scatter(3,1);
EffScatter = Scatter(3,2);

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
% NaNJ = [3:4];
% NumberZ(:,NaNJ) = NaN;

disp(['�ܶȷ������']);

%% ������
PSJ = 255
Density = NumberZ(KM70:KM125,PSJ);
Density(Density<0) = 0;
AltitudeY = Altitude(KM70:KM125,:);
plot(AltitudeY,Density,'-b', 'linewidth', 1.5);
xlim([70 125]); 
% ylim([0 20])
xlabel('Altitude (km)'); ylabel('Density (cm^{��3})')
legend('Ca^+'); grid on;
ttstr = ['(d) ',char(string(TimeList(PSJ)))];
title(ttstr)
set(gca,'fontsize',15);

%%
tFind = TimeList';

figure(1)
Density = NumberZ(KM75:KM350,:);
% Density(Density<2.5) = 0;
TimeX = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z');
TimeY = 1:size(TimeX,2);
hours = hour(TimeX);
minutes = minute(TimeX);
Time = (((hours(1,:)+8)*3600+minutes(1,:)*60)+((((hours(end,:)+8)*3600+minutes(end,:)*60)-((hours(1,:)+8)*3600+minutes(1,:)*60)).*(TimeY./TimeY(end))))./3600;
AltitudeY = Altitude(KM75:KM350,:);
plot(Density, AltitudeY, '-k', 'linewidth', 1.5);
ttstr1 = '(i) Ca ion ' + string(ReadDate) + ' 23:00 Mohe';
title(ttstr1, 'FontWeight', 'bold');
% xlim([75 125]);
ylim([75 115])
xlim([0 3000])
xlabel('Density (cm^{-3})', 'FontWeight', 'bold')
ylabel('Altitude (km)', 'FontWeight', 'bold')

%% ���ܶ�ʱ��ͼ,'position',[100 100 500 400]
FuHao = char(8722);

TimeC = datenum(TimeX)+8/24;

start_time = datenum('2023-07-20 23:00', 'yyyy-mm-dd HH:MM');
end_time = datenum('2023-7-21 02:00', 'yyyy-mm-dd HH:MM');
st_time = datenum('2023-07-20 23:00', 'yyyy-mm-dd HH:MM');
ed_time = datenum('2023-7-21 02:00', 'yyyy-mm-dd HH:MM');
ptTime = datenum('2023-11-03 23:00', 'yyyy-mm-dd HH:MM');
Ttk = 1/24;

figure('name','CaiTu')


window = fspecial('average',[1,1]);

% window = hann(5).*hann(5);

DenZ = imfilter(Density,window,'symmetric','same');
% DenZ(:,[1:15]) = 0;
DenZ(DenZ<=-0) = NaN;
DenZ(DenZ>=1000) = NaN;
[Map,Line] = contourf(TimeC(:,:),AltitudeY,DenZ(:,:),200);
xlim([start_time end_time])
x_ticks = st_time:Ttk:ed_time;
set(gca, 'XTick', x_ticks);
datetick('x','HH:MM','keeplimits','keepticks');
% xticks(18:1:23)
set(Line,'LineStyle', ':', 'LineWidth', 0.5, 'LineColor', 'none')
set(gca,'FontWeight', 'bold','FontSize',12,'fontname','arial');
ttstr1 = '(a) ' + string(ReadDate) + ' Mohe';
title(ttstr1, 'FontWeight', 'bold','FontSize',12,'fontname','arial');
set(gca, 'ColorScale', 'log');
colormap(jet); 
% colorbar;
% ylabel(colorbar,'Ca ion density (cm^{-3})','FontWeight', 'bold','FontSize',12,'fontname','arial');
xlabel('Time (LT)','FontWeight', 'bold','FontSize',12,'fontname','arial');
ylabel('Altitude (km)','FontWeight', 'bold','FontSize',12,'fontname','arial');
ylim([80 300]);
% set(gca,'YTick',[])

% hold on
% plot([ptTime ptTime],[75 115],'--w','linewidth',1)

%% �����ļ���
folderName = ['Proccessed\Cap\',ReadDate,'\'];  % ָ���ļ�������
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

    PsdName = [folderName,'OMOHE_MUCL01_NDCI_L2_STP_',WriteStr,'_V01.00.TXT'];
    fileID = fopen(PsdName, 'w');
    RecordNumber = size(data,1);
    % ��ȡ��ǰʱ��
    dt_now = datetime('now');
    dt_str_now = datestr(dt_now, 'yyyy-mm-ddTHH:MM:SS');
    % д���ͷ
    fprintf(fileID,'#DataName: Number Density of Ca Ions\n');
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
    fprintf(fileID,'#DeviceSpec: WaveLen=393nm, PRF=15Hz, PlsEnergy=90mJ\n');
    fprintf(fileID,'#ObsParameters: PlsAccu=13500\n');
    fprintf(fileID,'#Quantities: Number density of Ca ions (cm^-3)\n');
    fprintf(fileID,'#Elev(km): Height, F7.3, missingdata=NaN\n');
    fprintf(fileID,'#CaIonDen(cm^-3): Ionic Ordinal Density of Calcium, F9.1, missingdata=NaN\n');
    fprintf(fileID,'#DenEr(%%): Density Error, F5.1, missingdata=NaN\n');
    fprintf(fileID,'#---------------------------------\n');
    % д���������
    fprintf(fileID, '%+7s %+9s %+5s\n', 'Elev', 'CaIonDen', 'DenEr');
    % д�����ݱ��
    fprintf(fileID, '%7.3f %9.1f %5.1f\n', data.');
    % �ر��ļ�
    fclose(fileID);
end
