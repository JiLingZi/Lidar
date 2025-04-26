% Ray密度
% 读取光子数矩阵
ReadDate = '20231106';
fldstart = ['G:\ZWDATA\MOHEnew\K532\',ReadDate,'\'];
folder = [fldstart,'CH2\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;% 为了统一文件数，减一
% 减一是因为那个破机器动不动少一个文件，导致不晓得哪个通道文件数不同
% 格式不要再改了，探测高度上限也不要再改了，否则实现自动化反演困难重重

DenPh = zeros(4096, num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    fid = fopen(filename);
    PhData = textscan(fid, '%f %f','HeaderLines',20);   % 7月22日及之前均为23，后应改为20
    fclose(fid);
    DenPh(:,j) = PhData{1,2}(1:4096);
end

% 读取时间序列
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

% 合并时间文件数
jNum = 30;
TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

% 获取原始高度矩阵
height_num_origin = PhData{1,1}(1:4096);

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
KM40 = size(Altitude(Altitude<40),1)+1;
KM70 = size(Altitude(Altitude<70),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
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



% 基本常数
R = 8.314;                  % 热力学常数 J/(K*mol)
M = 28.959e-3;              % 大气摩尔质量 kg/mol
NA = 6.023e23;              % 阿伏伽德罗常数

% 重力加速度
G_all = 6.67259e-11;        % 万有引力常数
M_Earth = 5.965e24;         % 地球质量 kg
R_Earth = 6371004;          % 地球半径 m
G = (G_all*M_Earth)./((R_Earth+(Altitude.*1e3)).^2);

% 获取NRLMSIS大气模型温度与密度,Rho_Z_0取40km，T_Z_0取70km
Z_0_NRL = Altitude(KM40)*1e3;
Z_top = Altitude(KM70)*1e3;
[TemRay,DenRay] = atmosnrlmsise00([Z_0_NRL,Z_top],40.5,116.0,2023,306,17);
Rho_Z_0 = (sum(DenRay(1,:))./NA) * M;
T_Z_top = TemRay(2,2);

% 密度反演
Ph = PhjSum;              % 取20:00 光子数
Z = Altitude;        % 高度矩阵
Z_0 = Z_0_NRL*1e-3;               % 密度参考高度
N_Z = Ph;          % 30-60 km光子数
N_B = mean(Ph(KM120:KM125,:));  % 光子数平均噪声水平
% N_STD = std(Ph(126:147,:));
% Ratio_STD = 1-(N_STD./N_B);.*Ratio_STD
N_Z_0 = Ph(KM40,:);           % 密度参考高度处光子数
Rho_Z = (Rho_Z_0 .* ((Z.^2)./(Z_0^2)) .* ((N_Z-N_B)./(N_Z_0-N_B)));

% 密度相对误差
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


% 密度创建文件夹
folderName = ['Proccessed\RayDen\',ReadDate,'\'];  % 指定文件夹名称
mkdir(folderName);       % 创建文件夹

% 计算密度误差
DenError = 1./sqrt(PhjSum);
DenError(DenError>1)=NaN;
DataAltitude = Altitude(KM30:KM80,:);
for j = 1:size(TimeList,2)
    % 合并高度，密度，误差
    DataDensity = Rho_Z(KM30:KM80,j);
    for i = 1:size(Altitude,1)
        if Rho_Z(i,:)<0
            DenError(i,:)=NaN;
        end
    end
    DataDensity(DataDensity<0)=NaN;
    DataError = real(DenError(KM30:KM80,j)).*100;
    data = [DataAltitude,DataDensity,DataError];

    % 新建TXT文件
    datetimeStr = TimeList{j};  % 从单元格中获取日期时间字符串
    
    datetimeObj = datetime(datetimeStr, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');  % 将字符串转换为datetime对象
    WriteStr = datestr(datetimeObj, 'yyyymmddHHMMSS');  % 将datetime对象

    PsdName = [folderName,'OMOHE_MUCL01_DUSM_L2_STP_',WriteStr,'_V01.00.TXT'];
    fileID = fopen(PsdName, 'w');
    RecordNumber = size(data,1);
    % 获取当前时间
    dt_now = datetime('now');
    dt_str_now = datestr(dt_now, 'yyyy-mm-ddTHH:MM:SS');
    % 写入表头
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
    % 写入参数名称
    fprintf(fileID, '%+7s %+8s %+5s\n', 'Elev', 'DenRay', 'DenEr');
    % 写入数据表格
    fprintf(fileID, '%7.3f %8.5f %5.1f\n', data.');
    % 关闭文件
    fclose(fileID);
end



% 曲线
Curve = Rho_Z.*G;
% plot(Altitude,Curve);


% 温度反演
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

% 画温度廓线图
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


% 温度创建文件夹
folderName = ['Proccessed\RayTem\',ReadDate,'\'];  % 指定文件夹名称
mkdir(folderName);       % 创建文件夹

% 计算温度误差
DenError = Tem./sqrt(PhjSum(1:KM70-1,:));
DenError(DenError>20)=NaN;
DataAltitude = Altitude(KM30:KM70-1,:);
for j = 1:size(TimeList,2)
    % 合并高度，密度，误差
    DataDensity = Tem(KM30:KM70-1,j);
    for i = 1:size(Tem,1)
        if Tem(i,j)>350
            DenError(i,j)=NaN;
        end
    end
    DataDensity(DataDensity>350)=NaN;
    DataError = real(DenError(KM30:KM70-1,j));
    data = [DataAltitude,DataDensity,DataError];

    % 新建TXT文件
    datetimeStr = TimeList{j};  % 从单元格中获取日期时间字符串
    
    datetimeObj = datetime(datetimeStr, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');  % 将字符串转换为datetime对象
    WriteStr = datestr(datetimeObj, 'yyyymmddHHMMSS');  % 将datetime对象

    PsdName = [folderName,'OMOHE_MUCL01_TUSM_L2_STP_',WriteStr,'_V01.00.TXT'];
    fileID = fopen(PsdName, 'w');
    RecordNumber = size(data,1);
    % 获取当前时间
    dt_now = datetime('now');
    dt_str_now = datestr(dt_now, 'yyyy-mm-ddTHH:MM:SS');
    % 写入表头
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
    % 写入参数名称
    fprintf(fileID, '%+7s %+7s %+5s\n', 'Elev', 'TempRay', 'DenEr');
    % 写入数据表格
    fprintf(fileID, '%7.3f %7.1f %5.1f\n', data.');
    % 关闭文件
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












