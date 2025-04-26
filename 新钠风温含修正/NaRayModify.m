%% 读取 Na 光子数矩阵
% 日期数据
ReadDate = '20240128';
RunDate = ReadDate;

Year = str2double(char(RunDate(1:4)));
Month = str2double(char(RunDate(5:6)));
Day = str2double(char(RunDate(7:8)));
DayNum = Month*30-30+Day;

YYYYNum = Year;
MMNum = Month;
DDNum = Day;
DateNum = datenum(YYYYNum,MMNum,DDNum)-datenum(YYYYNum,1,1);


fldstart = ['F:\RawData\ZWDATA\MOHEnew\Na\',ReadDate,'\'];
folder = [fldstart,'Na\'];
files = dir(fullfile(folder, '*.dat'));
num_files = length(files)-1;% 为了统一文件数，减一

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

RunRead = '数据读取OK了'

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

%%

jNum = 1;
% jNum = 28;
iNum = 32;
iNumx = 37;

TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
	TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

TimeX = datetime(TimeList, 'InputFormat', ':yyyy-MM-dd''T''HH:mm:ss.SSS');

% 获取原始高度矩阵
height_num_origin = Na_data(1:8192,1);

% 合并行
VF_0i = zeros(floor(8192/iNum),num_files);
VF_Ri = zeros(floor(8192/iNum),num_files);
VF_Li = zeros(floor(8192/iNum),num_files);

NF_0i = zeros(floor(8192/iNumx),num_files);
NF_Ri = zeros(floor(8192/iNumx),num_files);
NF_Li = zeros(floor(8192/iNumx),num_files);
EF_0i = zeros(floor(8192/iNumx),num_files);
EF_Ri = zeros(floor(8192/iNumx),num_files);
EF_Li = zeros(floor(8192/iNumx),num_files);

for i = 1:floor(8192/iNum)
    VF_0i(i,:) = sum(V_F_0(iNum*(i-1)+1:iNum*i,:),1);
    VF_Ri(i,:) = sum(V_F_R(iNum*(i-1)+1:iNum*i,:),1);
    VF_Li(i,:) = sum(V_F_L(iNum*(i-1)+1:iNum*i,:),1);
end
    
for i = 1:floor(8192/iNumx)
    NF_0i(i,:) = sum(N_F_0(iNumx*(i-1)+1:iNumx*i,:),1);
    NF_Ri(i,:) = sum(N_F_R(iNumx*(i-1)+1:iNumx*i,:),1);
    NF_Li(i,:) = sum(N_F_L(iNumx*(i-1)+1:iNumx*i,:),1);
    EF_0i(i,:) = sum(E_F_0(iNumx*(i-1)+1:iNumx*i,:),1);
    EF_Ri(i,:) = sum(E_F_R(iNumx*(i-1)+1:iNumx*i,:),1);
    EF_Li(i,:) = sum(E_F_L(iNumx*(i-1)+1:iNumx*i,:),1);
end

% 合并列
VF_0 = zeros(floor(8192/iNum),floor(num_files/jNum));
VF_R = zeros(floor(8192/iNum),floor(num_files/jNum));
VF_L = zeros(floor(8192/iNum),floor(num_files/jNum));

NF_0 = zeros(floor(8192/iNumx),floor(num_files/jNum));
NF_R = zeros(floor(8192/iNumx),floor(num_files/jNum));
NF_L = zeros(floor(8192/iNumx),floor(num_files/jNum));
EF_0 = zeros(floor(8192/iNumx),floor(num_files/jNum));
EF_R = zeros(floor(8192/iNumx),floor(num_files/jNum));
EF_L = zeros(floor(8192/iNumx),floor(num_files/jNum));

for j = 1:floor(num_files/jNum)
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

% 噪声120-150，计算光子数，高度
Altitude = height_num_origin(1:floor(8192/iNum))*iNum;
KM30 = size(Altitude(Altitude<30),1)+1;
KM35 = size(Altitude(Altitude<35),1)+1;
KM40 = size(Altitude(Altitude<40),1)+1;
KM75 = size(Altitude(Altitude<75),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM85 = size(Altitude(Altitude<85),1)+1;
KM90 = size(Altitude(Altitude<90),1)+1;
KM95 = size(Altitude(Altitude<95),1)+1;
KM98 = size(Altitude(Altitude<98),1)+1;
KM100 = size(Altitude(Altitude<100),1)+1;
KM105 = size(Altitude(Altitude<105),1)+1;
KM115 = size(Altitude(Altitude<115),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM150 = size(Altitude(Altitude<150),1)+1;
Alt30 = (height_num_origin(1:floor(8192/iNumx))*iNumx)*(sqrt(3)/2);
KMx30 = size(Alt30(Alt30<30),1)+1;
KMx35 = size(Alt30(Alt30<35),1)+1;
KMx40 = size(Alt30(Alt30<40),1)+1;
KMx75 = size(Alt30(Alt30<75),1)+1;
KMx80 = size(Alt30(Alt30<80),1)+1;
KMx85 = size(Alt30(Alt30<85),1)+1;
KMx90 = size(Alt30(Alt30<90),1)+1;
KMx95 = size(Alt30(Alt30<95),1)+1;
KMx98 = size(Alt30(Alt30<98),1)+1;
KMx100 = size(Alt30(Alt30<100),1)+1;
KMx105 = size(Alt30(Alt30<105),1)+1;
KMx115 = size(Alt30(Alt30<115),1)+1;
KMx120 = size(Alt30(Alt30<120),1)+1;
KMx150 = size(Alt30(Alt30<150),1)+1;

% 去噪声
VF0 = VF_0 - mean(VF_0(KM120:KM150,:),1);
VFR = VF_R - mean(VF_R(KM120:KM150,:),1);
VFL = VF_L - mean(VF_L(KM120:KM150,:),1);
NF0 = NF_0 - mean(NF_0(KMx120:KMx150,:),1);
NFR = NF_R - mean(NF_R(KMx120:KMx150,:),1);
NFL = NF_L - mean(NF_L(KMx120:KMx150,:),1);
EF0 = EF_0 - mean(EF_0(KMx120:KMx150,:),1);
EFR = EF_R - mean(EF_R(KMx120:KMx150,:),1);
EFL = EF_L - mean(EF_L(KMx120:KMx150,:),1);

RunRead = '第二次合并去噪OK了'

% 归一化
FL = 1;
FR = size(VF0,2);
V0 = VF0(:,FL:FR)./VF0(KM35,FL:FR);
VR = VFR(:,FL:FR)./VFR(KM35,FL:FR);
VL = VFL(:,FL:FR)./VFL(KM35,FL:FR);
N0 = NF0(:,FL:FR)./NF0(KMx35,FL:FR);
NR = NFR(:,FL:FR)./NFR(KMx35,FL:FR);
NL = NFL(:,FL:FR)./NFL(KMx35,FL:FR);
E0 = EF0(:,FL:FR)./EF0(KMx35,FL:FR);
ER = EFR(:,FL:FR)./EFR(KMx35,FL:FR);
EL = EFL(:,FL:FR)./EFL(KMx35,FL:FR);

RaySV0 = sum(VF_0(KM35:KM40,:),1);
RaySVR = sum(VF_R(KM35:KM40,:),1);
RaySVL = sum(VF_L(KM35:KM40,:),1);
RaySN0 = sum(NF_0(KMx35:KMx40,:),1);
RaySNR = sum(NF_R(KMx35:KMx40,:),1);
RaySNL = sum(NF_L(KMx35:KMx40,:),1);
RaySE0 = sum(EF_0(KMx35:KMx40,:),1);
RaySER = sum(EF_R(KMx35:KMx40,:),1);
RaySEL = sum(EF_L(KMx35:KMx40,:),1);

Tlist1 = datestr(TimeX','yyyymmdd HH:MM');
Tlist = Tlist1(FL:FR,:);
RayV0 = RaySV0(:,FL:FR);
RayVR = RaySVR(:,FL:FR);
RayVL = RaySVL(:,FL:FR);
RayN0 = RaySN0(:,FL:FR);
RayNR = RaySNR(:,FL:FR);
RayNL = RaySNL(:,FL:FR);
RayE0 = RaySE0(:,FL:FR);
RayER = RaySER(:,FL:FR);
RayEL = RaySEL(:,FL:FR);

ErayVR = RayV0./RayVR;
ErayVL = RayV0./RayVL;
ErayNR = RayN0./RayNR;
ErayNL = RayN0./RayNL;
ErayER = RayE0./RayER;
ErayEL = RayE0./RayEL;

NoiseV0 = VF_0 - VF0;
NoiseVR = VF_R - VFR;
NoiseVL = VF_L - VFL;

RunRead = '第二次归一化OK了'


%% Ray温密度段
Alt_Res = 983;

for jpf = 1:size(VF0,2)
    pf = jpf;
    TlistRay = cellstr(Tlist);
    AltVX = Altitude;

    % 各个参考高度的行数
    i_40_mix = find(AltVX<40);
    i_70_mix = find(AltVX>70);
    i_120_mix = find(AltVX>120);
    i_140_mix = find(AltVX>140);
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
    G = (G_all*M_Earth)./((R_Earth+(AltVX.*1e3)).^2);

    % 获取NRLMSIS大气模型温度与密度,Rho_Z_0取40km，T_Z_0取70km
    Z_0_NRL = AltVX(i_40)*1e3;
    Z_top = AltVX(i_70)*1e3;
    [TemRay,DenRay] = atmosnrlmsise00([Z_0_NRL,Z_top],40.5,116.0,YYYYNum,DateNum,17);
    Rho_Z_0 = (sum(DenRay(1,:))./NA) * M;
    T_Z_top = TemRay(2,2);

    % 密度反演
    Ph = VF0(:,pf);              % 取20:00 光子数
    Ph = smoothdata(Ph, 'movmean', 5);
    Z = AltVX;        % 高度矩阵
    Z_0 = Z_0_NRL*1e-3;               % 密度参考高度
    N_Z = Ph;          % 30-60 km光子数
    N_B = mean(Ph(i_120:i_140,:));  % 光子数平均噪声水平
    % N_STD = std(Ph(126:147,:));
    % Ratio_STD = 1-(N_STD./N_B);.*Ratio_STD
    N_Z_0 = Ph(i_40,:);           % 密度参考高度处光子数
    Rho_ZNa = (Rho_Z_0 .* ((Z.^2)./(Z_0^2)) .* ((N_Z-N_B)./(N_Z_0-N_B)));

    % 密度相对误差
    DenErr = Rho_ZNa.*(1./sqrt(N_Z));

    % 曲线
    Curve = Rho_ZNa.*G;
    % plot(Altitude,Curve);

    % 温度反演
    % Rho_Z_top = (sum(DenRay(2,:))./NA) * M;
    Rho_Z_top = Rho_ZNa(i_70,:);
    T_ZNa = (1:(i_70-1))';
    for i = 1:(i_70-1)
        T_ZNa(i,:) = (Rho_Z_top/Rho_ZNa(i,:))*T_Z_top + ((M/R)/Rho_ZNa(i,:))*((sum(Curve(i:i_70,:)))*Alt_Res);
    end

    AltTZ = AltVX(1:i_70-1,:);
    
    fC = figure('name','NaRayTem');
    plot(T_ZNa,AltTZ,'-b','linewidth',2)
    ylim([25 65])
    xlabel('Temperature (K)');
    ylabel('Altitude (km)');
    ttstr = "Rayleigh Temperature of Na " + string(TlistRay(pf));
    title(ttstr);
    box on;
    grid on;

    figPath = "D:\File\20241010PPT\NaRay\"+ReadDate+"\";
    if ~exist(figPath, 'dir')
        mkdir(figPath);
    end
    figV = datestr(TimeX(:,jpf),'yyyymmdd')+"-NaRay-"+string(jpf);
    print(fC, '-dpng', '-r300', figPath+figV);
    close all;
end

%% 单条Ray温密度段
Alt_Res = 983;

for jpf = 488
    pf = jpf;
    TlistRay = cellstr(Tlist);
    AltVX = Alt30;

    % 各个参考高度的行数
    i_40_mix = find(AltVX<40);
    i_70_mix = find(AltVX>50);
    i_120_mix = find(AltVX>120);
    i_140_mix = find(AltVX>140);
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
    G = (G_all*M_Earth)./((R_Earth+(AltVX.*1e3)).^2);

    % 获取NRLMSIS大气模型温度与密度,Rho_Z_0取40km，T_Z_0取70km
    Z_0_NRL = AltVX(i_40)*1e3;
    Z_top = AltVX(i_70)*1e3;
    [TemRay,DenRay] = atmosnrlmsise00([Z_0_NRL,Z_top],40.5,116.0,YYYYNum,DateNum,17);
    Rho_Z_0 = (sum(DenRay(1,:))./NA) * M;
    T_Z_top = TemRay(2,2);

    % 密度反演
    Ph = EF0(:,pf);              % 取20:00 光子数
    Ph = smoothdata(Ph, 'movmean', 5);
    PhNa = Ph;
    Z = AltVX;        % 高度矩阵
    Z_0 = Z_0_NRL*1e-3;               % 密度参考高度
    N_Z = Ph;          % 30-60 km光子数
    N_B = mean(Ph(i_120:i_140,:));  % 光子数平均噪声水平
    % N_STD = std(Ph(126:147,:));
    % Ratio_STD = 1-(N_STD./N_B);.*Ratio_STD
    N_Z_0 = Ph(i_40,:);           % 密度参考高度处光子数
    Rho_ZNa = (Rho_Z_0 .* ((Z.^2)./(Z_0^2)) .* ((N_Z-N_B)./(N_Z_0-N_B)));

    % 密度相对误差
    DenErr = Rho_ZNa.*(1./sqrt(N_Z));

    % 曲线
    Curve = Rho_ZNa.*G;
    % plot(Altitude,Curve);

    % 温度反演
    % Rho_Z_top = (sum(DenRay(2,:))./NA) * M;
    Rho_Z_top = Rho_ZNa(i_70,:);
    T_ZNa = (1:(i_70-1))';
    for i = 1:(i_70-1)
        T_ZNa(i,:) = (Rho_Z_top/Rho_ZNa(i,:))*T_Z_top + ((M/R)/Rho_ZNa(i,:))*((sum(Curve(i:i_70,:)))*Alt_Res);
    end

    AltTZ = AltVX(1:i_70-1,:);
    
    fC = figure('name','NaRayTem');
    plot(T_ZNa,AltTZ,'-b','linewidth',2)
    hold on;
    BHalt = 30.55;
    plot([220 340],[BHalt BHalt],'--r','linewidth',2)
    ylim([25 55])
    xlabel('Temperature (K)');
    ylabel('Altitude (km)');
    ttstr = "Rayleigh Temperature of Na " + string(TlistRay(pf));
    title(ttstr);
    box on;
    grid on;
    hold on;
end


%% 在饱和线段
AltOg30 = height_num_origin*(sqrt(3)/2);
AltOg = height_num_origin;
pbhf = 488;
BHct = 105;
figure('name','饱和线','position',[50 50 800 800])
subplot(3,1,1)
plot(AltOg30,E_F_0(:,pbhf),'-b','linewidth',1.5);
hold on;
plot([20 115],[BHct BHct],'--r','linewidth',1);
xlim([20 115])
grid on;
title('f0')
xlabel('Altitude (km)')
ylabel('Counts')
% set(gca, 'YScale', 'log');

subplot(3,1,2)
plot(AltOg30,E_F_R(:,pbhf),'-b','linewidth',1.5);
hold on;
plot([20 115],[BHct BHct],'--r','linewidth',1);
xlim([20 115])
grid on;
title('f+')
xlabel('Altitude (km)')
ylabel('Counts')
% set(gca, 'YScale', 'log');

subplot(3,1,3)
plot(AltOg30,E_F_L(:,pbhf),'-b','linewidth',1.5);
hold on;
plot([20 115],[BHct BHct],'--r','linewidth',1);
xlim([20 115])
grid on;
title('f-')
xlabel('Altitude (km)')
ylabel('Counts')
% set(gca, 'YScale', 'log');



























