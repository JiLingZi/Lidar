% Ca密度
% 读取光子数矩阵,钙原子、钙离子有4个通道
ReadDate = '20231103';
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Ca\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;% 为了统一文件数，减一
% 减一是因为那个破机器动不动少一个文件，导致不晓得哪个通道文件数不同
% 格式不要再改了，探测高度上限也不要再改了，否则实现自动化反演困难重重


HL = 20;

DenPh1 = zeros(16384, num_files);
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
for j = 1:num_files
    filename = fullfile(folder,files(j).name);
    fid = fopen(filename);
    PhData1 = textscan(fid, '%f %f','HeaderLines',HL);
    fclose(fid);
    
    DenPh1(:,j) = PhData1{1,2}(1:16384,:);
end

DenPh2 = zeros(16384, num_files);
folder = [fldstart,'CH2\'];
files = dir(fullfile(folder, '*.txt'));
for j = 1:num_files
    filename = fullfile(folder,files(j).name);
    fid = fopen(filename);
    PhData2 = textscan(fid, '%f %f','HeaderLines',HL);
    fclose(fid);
    DenPh2(:,j) = PhData2{1,2}(1:16384,:);
end

DenPh3 = zeros(16384, num_files);
folder = [fldstart,'CH3\'];
files = dir(fullfile(folder, '*.txt'));
for j = 1:num_files
    filename = fullfile(folder,files(j).name);
    fid = fopen(filename);
    PhData3 = textscan(fid, '%f %f','HeaderLines',HL);
    fclose(fid);
    DenPh3(:,j) = PhData3{1,2}(1:16384,:);
end

DenPh4 = zeros(16384, num_files);
folder = [fldstart,'CH4\'];
files = dir(fullfile(folder, '*.txt'));
for j = 1:num_files
    filename = fullfile(folder,files(j).name);
    fid = fopen(filename);
    PhData4 = textscan(fid, '%f %f','HeaderLines',HL);
    fclose(fid);
    DenPh4(:,j) = PhData4{1,2}(1:16384,:);
end

DenPh = DenPh1+DenPh2+DenPh3+DenPh4;

%% 读取时间序列
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

%% 合并时间文件数
jNum = 10;
% jNum = 1;
TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

% 获取原始高度矩阵
height_num_origin = PhData1{1,1};

% 合并行，3行合并，96米
PhiSum = zeros(floor(16384/3),num_files);
for i = 1:floor(16384/3)
    PhiSum(i,:) = sum(DenPh(3*(i-1)+2:3*i+1,:),1);
end

% 合并列，15列合并，15分钟
PhjSum = zeros(floor(16384/3),floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    PhjSum(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
end

% 噪声120-150，计算光子数，高度
Altitude = height_num_origin(1:floor(16384/3))*3;
KM30 = size(Altitude(Altitude<30),1)+1;
KM75 = size(Altitude(Altitude<75),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
KM115 = size(Altitude(Altitude<115),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM125 = size(Altitude(Altitude<125),1)+1;
KM180 = size(Altitude(Altitude<180),1)+1;
KM300 = size(Altitude(Altitude<300),1)+1;
KM350 = size(Altitude(Altitude<350),1)+1;
KM400 = size(Altitude(Altitude<400),1)+1;

% 归一化光子数检查
CheckH = 50;% 光子数检查归一化高度，KM
KMCheck = size(Altitude(Altitude<CheckH),1);
CheckPh = PhjSum./PhjSum(KMCheck,:);
figure('Name','光子数检查','position',[10 100 1350 500])
for j = 1:size(CheckPh,2)
    plot(CheckPh(:,j)+j,Altitude,'-','linewidth',1.5);
    hold on;
end
xlim([0 31]);
ylim([75 115]);
xInterval = 1;  % x轴刻度间隔

% 噪声
Noise = mean(PhjSum(KM350:KM400,:),1);
SumPh = PhjSum - Noise;

SNR = SumPh./Noise;

% 获取瑞利后向散射截面和有效后向散射截面
ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));
RayScatter = Scatter(4,1);
EffScatter = Scatter(4,2);

% 获取30 km 处大气模型温度与密度
[TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,2023,201,17);

% 反演密度
Z = Altitude;                       % 高度矩阵
ZR = 30;                            % 参考高度
SigmaRay = RayScatter;              % 瑞利后向散射截面
SigmaEff = EffScatter;              % 有效后向散射截面
NZ = PhjSum;                        % 光子数矩阵
NB = Noise;                         % 噪声矩阵
NZR = PhjSum(KM30,:);                % 参考高度处光子数
NumberRay = sum(DenRay)*1e-6;       % 参考高度处大气模型数密度
NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

% 密度廓线粗检查
% NaNJ = [1];
% NumberZ(:,NaNJ) = NaN;

%%
figure(1)

tFIND = TimeList';

Density = NumberZ(KM75:KM350,:);
Density(Density<0) = 0;
TimeX = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
TimeY = 1:size(TimeX,2);
hours = hour(TimeX);
minutes = minute(TimeX);
Time = (((hours(1,:)+8)*3600+minutes(1,:)*60)+((((hours(end,:)+8)*3600+minutes(end,:)*60)-((hours(1,:)+8)*3600+minutes(1,:)*60)).*(TimeY./TimeY(end))))./3600;
AltitudeY = Altitude(KM75:KM350,:);
plot(Density, AltitudeY, '-k', 'linewidth', 1.5);
title('Ca Density (cm^{-3}) 20230722 01:13 Mohe', 'FontWeight', 'bold');
% xlim([75 125]);
ylim([75 150])
xlabel('Density (cm^{-3})', 'FontWeight', 'bold')
ylabel('Altitude (km)', 'FontWeight', 'bold')

%% 画密度时变图
figure('Name','Varriation')
% Density(:,24:end) = 0;
window = fspecial('average',[1,1]);
DenZ = imfilter(Density,window,'symmetric','same');
[Map,Line] = contourf(Time,AltitudeY,DenZ,200);
% datetick('x', 'HH:MM')
set(Line,'LineColor','none')
set(gca,'FontSize',12,'FontName','arial');
% title(fldstart);
title('(h) 20231103 Mohe');
% set(gca, 'ColorScale', 'log');
colormap(jet); 
colorbar;
ylabel(colorbar,'Ca density (cm^{-3})');
xlabel('Time (LT)');
ylabel('Altitude (km)');
ylim([75 115]);
xlim([19 25])
hold on
plot([23 23],[75 115],'--w','linewidth',1)

%% 创建文件夹
folderName = ['Proccessed\Ca\',ReadDate,'\'];  % 指定文件夹名称
mkdir(folderName);       % 创建文件夹

% 计算密度误差
DenError = 1./sqrt(PhjSum);
DataAltitude = Altitude(870:1195,:);

for j = 1:size(TimeList,2)
    % 合并高度，密度，误差
    DataDensity = NumberZ(870:1195,j);
    DataDensity(DataDensity<0)=0;
    DataError = real(DenError(870:1195,j)).*100;
    data = [DataAltitude,DataDensity,DataError];

    % 新建TXT文件
    datetimeStr = TimeList{j};  % 从单元格中获取日期时间字符串
    
    datetimeObj = datetime(datetimeStr, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');  % 将字符串转换为datetime对象
    WriteStr = datestr(datetimeObj, 'yyyymmddHHMMSS');  % 将datetime对象

    PsdName = [folderName,'OMOHE_MUCL01_NDCA_L2_STP_',WriteStr,'_V01.00.TXT'];
    fileID = fopen(PsdName, 'w');
    RecordNumber = size(data,1);
    % 获取当前时间
    dt_now = datetime('now');
    dt_str_now = datestr(dt_now, 'yyyy-mm-ddTHH:MM:SS');
    % 写入表头
    fprintf(fileID,'#DataName: Number Density of Ca Atoms\n');
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
    fprintf(fileID,'#DeviceSpec: WaveLen=423nm, PRF=15Hz, PlsEnergy=90mJ\n');
    fprintf(fileID,'#ObsParameters: PlsAccu=13500\n');
    fprintf(fileID,'#Quantities: Number density of Ca atoms (cm^-3)\n');
    fprintf(fileID,'#Elev(km): Height, F7.3, missingdata=NaN\n');
    fprintf(fileID,'#CaDen(cm^-3): Atomic Number Density of Calcium, F9.1, missingdata=NaN\n');
    fprintf(fileID,'#DenEr(%%): Density Error, F5.1, missingdata=NaN\n');
    fprintf(fileID,'#---------------------------------\n');
    % 写入参数名称
    fprintf(fileID, '%+7s %+9s %+5s\n', 'Elev', 'CaDen', 'DenEr');
    % 写入数据表格
    fprintf(fileID, '%7.3f %9.1f %5.1f\n', data.');
    % 关闭文件
    fclose(fileID);
end


