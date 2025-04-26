%% 读取 比较粒子 光子数矩阵
% 日期数据
ReadDate = '20240512';
RunDate = ReadDate;

Year = str2double(char(RunDate(1:4)));
Month = str2double(char(RunDate(5:6)));
Day = str2double(char(RunDate(7:8)));
DayNum = Month*30-30+Day;

YYYYNum = Year;
MMNum = Month;
DDNum = Day;
DateNum = datenum(YYYYNum,MMNum,DDNum)-datenum(YYYYNum,1,1);

fldstart = ['F:\RawData\ZWDATA\MOHEnew\K532\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;


DenPh = zeros(8192, num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    fid = fopen(filename);
    PhData = textscan(fid, '%f %f','HeaderLines',20);   % 7月22日及之前均为23，后应改为20
    fclose(fid);
    DenPh(:,j) = PhData{1,2}(1:8192);
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
jNum = 54;
iNum = 32;
TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

% 获取原始高度矩阵
height_num_origin = PhData{1,1}(1:8192);

% 合并行，3行合并，96米
PhiSum = zeros(floor(8192/iNum),num_files);
for i = 1:floor(8192/iNum)
    PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
end

% 合并列，15列合并，15分钟
PhjSum = zeros(floor(8192/iNum),floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    PhjSum(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
end

% 噪声120-150，计算光子数，高度
Altitude = height_num_origin(1:floor(8192/iNum))*iNum;
KM30 = size(Altitude(Altitude<30),1)+1;
KM75 = size(Altitude(Altitude<75),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
KM115 = size(Altitude(Altitude<115),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM125 = size(Altitude(Altitude<125),1)+1;


% 噪声
Noise = mean(PhjSum(KM120:KM125,:),1);
SumPh = PhjSum - Noise;
SumPf = SumPh(:,1);

RunRead = '数据读取OK了'

figure('name','Other原始光子数','position',[50 50 1200 267])
plot(Altitude,SumPf,'-b','linewidth',2)
grid on;
xlim([20 115])
ylim([0 12e6])
xlabel('Altitude (km)');
ylabel('Counts')

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

%% Ray温密度段
Alt_Res = 983;
pf =1;
TlistRay = cellstr(TimeList)';

% 各个参考高度的行数
i_40_mix = find(Altitude<40);
i_70_mix = find(Altitude>70);
i_120_mix = find(Altitude>120);
i_140_mix = find(Altitude>140);
i_40 = i_40_mix(end);
i_70 = i_70_mix(1);
i_120 = i_120_mix(1);
i_140 = i_140_mix(1);

% i_70_mix2 = find(Altitude2>70);
% i2_70 = i_70_mix2(1);

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
Z_0_NRL = Altitude(i_40)*1e3;
Z_top = Altitude(i_70)*1e3;
[TemRay,DenRay] = atmosnrlmsise00([Z_0_NRL,Z_top],40.5,116.0,YYYYNum,DateNum,17);
Rho_Z_0 = (sum(DenRay(1,:))./NA) * M;
T_Z_top = TemRay(2,2);

% 密度反演
Ph = SumPh(:,pf);              % 取20:00 光子数
Ph = smoothdata(Ph, 'movmean', 5);
PhOt = Ph;
Z = Altitude;        % 高度矩阵
Z_0 = Z_0_NRL*1e-3;               % 密度参考高度
N_Z = Ph;          % 30-60 km光子数
N_B = mean(Ph(i_120:i_140,:))+std(Ph(i_120:i_140,:));  % 光子数平均噪声水平
% N_STD = std(Ph(i_120:i_140,:));
% Ratio_STD = 1-(N_STD./N_B);.*Ratio_STD
N_Z_0 = Ph(i_40,:);           % 密度参考高度处光子数
Rho_ZCheck = (Rho_Z_0 .* ((Z.^2)./(Z_0^2)) .* ((N_Z-N_B)./(N_Z_0-N_B)));

% 密度相对误差
DenErr = Rho_ZCheck.*(1./sqrt(N_Z));

% 曲线
Curve = Rho_ZCheck.*G;
% plot(Altitude,Curve);

% 温度反演
% Rho_Z_top = (sum(DenRay(2,:))./NA) * M;
Rho_Z_top = Rho_ZCheck(i_70,:);
T_ZCheck = (1:(i_70-1))';
for i = 1:(i_70-1)
    T_ZCheck(i,:) = (Rho_Z_top/Rho_ZCheck(i,:))*T_Z_top + ((M/R)/Rho_ZCheck(i,:))*((sum(Curve(i:i_70,:)))*Alt_Res);
end

AltTZ = Altitude(1:i_70-1,:);


plot(T_ZCheck,AltTZ,'-r','linewidth',2)
ylim([25 65])
xlabel('Temperature (K)');
ylabel('Altitude (km)');
ttstr = "Rayleigh Temperature of 532 " + string(TlistRay(pf));
title(ttstr);
box on;
grid on;




