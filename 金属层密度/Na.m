% Na密度
% 读取光子数矩阵
ReadDate = '20231103';
fldstart = ['F:\RawData\ZWDATA\MoheNEW\Na\',ReadDate,'\'];
folder = [fldstart,'Na\'];
files = dir(fullfile(folder, '*.dat'));
num_files = length(files)-1;% 为了统一文件数，减一
% 减一是因为那个破机器动不动少一个文件，导致不晓得哪个通道文件数不同
% 格式不要再改了，探测高度上限也不要再改了，否则实现自动化反演困难重重



DenPh = zeros(8192, num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    Na_table = readtable(filename);
    Na_data = table2array(Na_table);
    DenPh(:,j) = Na_data(:,2);
end
height_num = Na_data(:,1);

% 读取时间序列
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

% 合并时间文件数
jNum = 1;
TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

% 获取原始高度矩阵
height_num_origin = height_num(1:4096);

% 合并行，3行合并，96米
PhiSum = zeros(floor(4096/3),num_files);
for i = 1:floor(4096/3)
    PhiSum(i,:) = sum(DenPh(3*(i-1)+2:3*i+1,:),1);
end

% 合并列，15列合并，15分钟
PhjSum = zeros(floor(4096/3),floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    PhjSum(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
end

% 噪声120-150，计算光子数，高度
Altitude = height_num_origin(1:floor(4096/3))*3;
KM30 = size(Altitude(Altitude<30),1)+1;
KM75 = size(Altitude(Altitude<75),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
KM115 = size(Altitude(Altitude<115),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM125 = size(Altitude(Altitude<125),1)+1;

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
Noise = mean(PhjSum(KM120:KM125,:),1);
SumPh = PhjSum - Noise;

SNR = SumPh./Noise;

% 获取瑞利后向散射截面和有效后向散射截面，散射截面文件应放在同一文件夹下，或自行修改路径
ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));
RayScatter = Scatter(5,1);
EffScatter = Scatter(5,2);

% 获取30 km 处大气模型温度与密度，（高度米，纬度，经度，年，年中第几天，世界时）
[TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,2023,201,17);

% 反演密度
Z = Altitude;                       % 高度矩阵
ZR = 30;                            % 参考高度
SigmaRay = RayScatter;              % 瑞利后向散射截面
SigmaEff = EffScatter;              % 有效后向散射截面
NZ = PhjSum;                        % 光子数矩阵
NB = Noise;                         % 噪声矩阵
NZR = PhjSum(KM30,:);               % 参考高度处光子数
NumberRay = sum(DenRay)*1e-6;       % 参考高度处大气模型数密度
NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

% 密度廓线粗检查，取80-105km
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



% 画密度时变图
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




% 创建文件夹
folderName = ['Proccessed\Na\',ReadDate,'\'];  % 指定文件夹名称
mkdir(folderName);       % 创建文件夹

% 计算密度误差
DenError = 1./sqrt(PhjSum);
DenError(DenError>1)=NaN;
DataAltitude = Altitude(870:1195,:);
for j = 1:size(TimeList,2)
    % 合并高度，密度，误差
    DataDensity = NumberZ(870:1195,j);
    for i = 1:size(Altitude,1)
        if NumberZ(i,:)<0
            DenError(i,:)=NaN;
        end
    end
    DataDensity(DataDensity<0)=NaN;
    DataError = real(DenError(870:1195,j)).*100;
    data = [DataAltitude,DataDensity,DataError];

    % 新建TXT文件
    datetimeStr = TimeList{j};  % 从单元格中获取日期时间字符串
    
    datetimeObj = datetime(datetimeStr, 'InputFormat', ':yyyy-MM-dd''T''HH:mm:ss.SSS');  % 将字符串转换为datetime对象
    WriteStr = datestr(datetimeObj, 'yyyymmddHHMMSS');  % 将datetime对象

    PsdName = [folderName,'OMOHE_MUCL01_NDNA_L2_STP_',WriteStr,'_V01.00.TXT'];
    fileID = fopen(PsdName, 'w');
    RecordNumber = size(data,1);
    % 获取当前时间
    dt_now = datetime('now');
    dt_str_now = datestr(dt_now, 'yyyy-mm-ddTHH:MM:SS');
    % 写入表头
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
    % 写入参数名称
    fprintf(fileID, '%+7s %+9s %+5s\n', 'Elev', 'NaDen', 'DenEr');
    % 写入数据表格
    fprintf(fileID, '%7.3f %9.1f %5.1f\n', data.');
    % 关闭文件
    fclose(fileID);
end


% %polarPcolor(距离,角度,需要绘制的变量,'Ncircles'，圆圈个数,'Nspokes'半径个数)
% polarPcolor(r,theta,ratial,'Ncircles',4,'Nspokes',13')
% %设置色标颜色
% colormap jet;
% %设置色标范围
% caxis([-10,10]);
% colorbar off %不显示自带的坐标
% c=colorbar;
% set(c,'YTick',-10:5:10); %色标值范围及显示间隔
% set(c,'YTickLabel',{'-10','-5','0','5','10'}) %具体刻度赋值
% ylabel(c,'radial wind speed (m/s)');
% hold on

