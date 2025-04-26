% Date--1
% Time--2
% Rank--3
% NaV---4
% NaM---5
% NaZ---6
% K-----7
% Fe----8
% Ni----9
% Ca----10
% Cap---11
% Ray---12

clc; clear all;

%%

RunDate = '20240119';
ReadDate = RunDate;
Year = str2double(char(RunDate(1:4)));
Month = str2double(char(RunDate(5:6)));
Day = str2double(char(RunDate(7:8)));
DayNum = Month*30-30+Day;

NajNum = 1; NaiNum = 1;NaDex = 1/1;
KjNum = 1; KiNum = 1;KDex = 1/2;
FejNum = 1; FeiNum = 1;FeDex = 1/2;
NijNum = 1; NiiNum = 1;NiDex = 1/1;
CajNum = 1; CaiNum = 1;CaDex = 1/1;
CapjNum = 1; CapiNum = 1;CapDex = 1/1;
RayjNum = 1; RayiNum = 1;RayDex = KDex;

%% Na

TS = 'Na Start '

fldstart = ['F:\RawData\ZWDATA\MOHEnew\Na\',ReadDate,'\'];
folder = [fldstart,'Na\'];
files = dir(fullfile(folder, '*.dat'));
num_files = length(files)-1;
num_files = 3*num_files;

FileNa = num_files;

if FileNa>=0
    V_F_0 = zeros(4096, num_files/3);
    V_F_R = zeros(4096, num_files/3);
    V_F_L = zeros(4096, num_files/3);
    N_F_0 = zeros(4096, num_files/3);
    N_F_R = zeros(4096, num_files/3);
    N_F_L = zeros(4096, num_files/3);
    E_F_0 = zeros(4096, num_files/3);
    E_F_R = zeros(4096, num_files/3);
    E_F_L = zeros(4096, num_files/3);
    for j = 1:(num_files/3)
        filename = fullfile(folder, files(j).name);
        Na_table = readtable(filename);
        Na_data = table2array(Na_table);
        V_F_0(:,j) = Na_data(1:4096,2);
        V_F_R(:,j) = Na_data(1:4096,3);
        V_F_L(:,j) = Na_data(1:4096,4);
        N_F_0(:,j) = Na_data(1:4096,5);
        N_F_R(:,j) = Na_data(1:4096,6);
        N_F_L(:,j) = Na_data(1:4096,7);
        E_F_0(:,j) = Na_data(1:4096,8);
        E_F_R(:,j) = Na_data(1:4096,9);
        E_F_L(:,j) = Na_data(1:4096,10);
    end
    DenPhV = zeros(4096, num_files);
    DenPhN = zeros(4096, num_files);
    DenPhE = zeros(4096, num_files);
    for jD = 1:(num_files/3)
        DenPhV(:,3*(jD-1)+1) = V_F_0(:,jD);
        DenPhV(:,3*(jD-1)+2) = V_F_R(:,jD);
        DenPhV(:,3*(jD-1)+3) = V_F_L(:,jD);
        DenPhN(:,3*(jD-1)+1) = N_F_0(:,jD);
        DenPhN(:,3*(jD-1)+2) = N_F_R(:,jD);
        DenPhN(:,3*(jD-1)+3) = N_F_L(:,jD);
        DenPhE(:,3*(jD-1)+1) = E_F_0(:,jD);
        DenPhE(:,3*(jD-1)+2) = E_F_R(:,jD);
        DenPhE(:,3*(jD-1)+3) = E_F_L(:,jD);
    end
    TimeValues = cell(1,num_files);
    for j = 1:(num_files/3)
        filename = fullfile(folder, files(j).name);
        fid = fopen(filename);
        TXT = textscan(fid, '%s', 'Delimiter', '\n');
        RowTime = TXT{1}{4};
        TimeData = strsplit(RowTime);
        fclose(fid);
        TimeValues(:,3*(j-1)+1) = cellstr(TimeData{4});
        TimeValues(:,3*(j-1)+2) = cellstr(TimeData{4});
        TimeValues(:,3*(j-1)+3) = cellstr(TimeData{4});
    end
    jNum = NajNum;
    iNum = NaiNum;
    TimeList = cell(1,floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
    end
    TimeXNa = datetime(TimeList, 'InputFormat', ':yyyy-MM-dd''T''HH:mm:ss.SSS');
    height_num_origin = Na_data(1:4096,1);
    
    % V
    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPhV(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumNaV = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumNaV(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeNaV = height_num_origin(1:floor(4096/iNum))*iNum;
    KM30NaV = size(AltitudeNaV(AltitudeNaV<30),1)+1;
    KM75NaV = size(AltitudeNaV(AltitudeNaV<75),1)+1;
    KM80NaV = size(AltitudeNaV(AltitudeNaV<80),1)+1;
    KM110NaV = size(AltitudeNaV(AltitudeNaV<110),1)+1;
    KM115NaV = size(AltitudeNaV(AltitudeNaV<115),1)+1;
    KM120NaV = size(AltitudeNaV(AltitudeNaV<120),1)+1;
    KM125NaV = size(AltitudeNaV(AltitudeNaV<125),1)+1;
    
    % N
    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPhN(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumNaM = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumNaM(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeNaM = height_num_origin(1:floor(4096/iNum))*iNum*(sqrt(3)/2);
    KM30NaM = size(AltitudeNaM(AltitudeNaM<30),1)+1;
    KM75NaM = size(AltitudeNaM(AltitudeNaM<75),1)+1;
    KM80NaM = size(AltitudeNaM(AltitudeNaM<80),1)+1;
    KM110NaM = size(AltitudeNaM(AltitudeNaM<110),1)+1;
    KM115NaM = size(AltitudeNaM(AltitudeNaM<115),1)+1;
    KM120NaM = size(AltitudeNaM(AltitudeNaM<120),1)+1;
    KM125NaM = size(AltitudeNaM(AltitudeNaM<125),1)+1;
    
    % E
    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPhE(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumNaZ = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumNaZ(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeNaZ = height_num_origin(1:floor(4096/iNum))*iNum*(sqrt(3)/2);
    KM30NaZ = size(AltitudeNaZ(AltitudeNaZ<30),1)+1;
    KM75NaZ = size(AltitudeNaZ(AltitudeNaZ<75),1)+1;
    KM80NaZ = size(AltitudeNaZ(AltitudeNaZ<80),1)+1;
    KM110NaZ = size(AltitudeNaZ(AltitudeNaZ<110),1)+1;
    KM115NaZ = size(AltitudeNaZ(AltitudeNaZ<115),1)+1;
    KM120NaZ = size(AltitudeNaZ(AltitudeNaZ<120),1)+1;
    KM125NaZ = size(AltitudeNaZ(AltitudeNaZ<125),1)+1;
end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Na ok '

%% Na 密度反演

TimeShuNa = TimeXNa';

Noise = mean(PhjSumNaV(KM120NaV:KM125NaV,:),1);

ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));
RayScatter = Scatter(5,1);
EffScatter = Scatter(5,2);

[TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,Year,DayNum,17);

Z = AltitudeNaV;                       % 高度矩阵
ZR = 30;                            % 参考高度
SigmaRay = RayScatter;              % 瑞利后向散射截面
SigmaEff = EffScatter;              % 有效后向散射截面
NZ = PhjSumNaV;                        % 光子数矩阵
NB = Noise;                         % 噪声矩阵
NZR = PhjSumNaV(KM30NaV,:);                % 参考高度处光子数
NumberRay = sum(DenRay)*1e-6;       % 参考高度处大气模型数密度
NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

DenNaV = NumberZ;

% 密度廓线

Prof = 4822;
Ns = NumberZ(:,Prof);
plot(Z,Ns,'-','linewidth',1.5)
grid on;
Profstr = string(Prof)+" | "+string(TimeShuNa(Prof));
title(Profstr)
xlim([70 125]);
% ylim([0 80]);
xlabel('Altitude (km)');
ylabel('Density (cm^{-3})');
legend('Na');

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Na Prof ok '

%% K

TS = 'K Start '

fldstart = ['F:\RawData\ZWDATA\MOHEnew\K532\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;

FileK = num_files;

if FileK>=0

    DenPh = zeros(4096, num_files);
    for j = 1:num_files
        filename = fullfile(folder, files(j).name);
        fid = fopen(filename);
        PhData = textscan(fid, '%f %f','HeaderLines',20);   % 7月22日及之前均为23，后应改为20
        fclose(fid);
        DenPh(:,j) = PhData{1,2}(1:4096);
    end
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
    jNum = KjNum;
    iNum = KiNum;
    TimeList = cell(1,floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
    end

    TimeXK = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
    height_num_origin = PhData{1,1}(1:4096);

    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumK = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumK(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeK = height_num_origin(1:floor(4096/iNum))*iNum;
    KM30K = size(AltitudeK(AltitudeK<30),1)+1;
    KM75K = size(AltitudeK(AltitudeK<75),1)+1;
    KM80K = size(AltitudeK(AltitudeK<80),1)+1;
    KM110K = size(AltitudeK(AltitudeK<110),1)+1;
    KM115K = size(AltitudeK(AltitudeK<115),1)+1;
    KM120K = size(AltitudeK(AltitudeK<120),1)+1;
    KM125K = size(AltitudeK(AltitudeK<125),1)+1;

end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'K ok '

%% K 密度反演

TimeShuK = TimeXK';

Noise = mean(PhjSumK(KM120K:KM125K,:),1);

ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));
RayScatter = Scatter(6,1);
EffScatter = Scatter(6,2);

[TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,Year,DayNum,17);

Z = AltitudeK;                       % 高度矩阵
ZR = 30;                            % 参考高度
SigmaRay = RayScatter;              % 瑞利后向散射截面
SigmaEff = EffScatter;              % 有效后向散射截面
NZ = PhjSumK;                        % 光子数矩阵
NB = Noise;                         % 噪声矩阵
NZR = PhjSumK(KM30K,:);                % 参考高度处光子数
NumberRay = sum(DenRay)*1e-6;       % 参考高度处大气模型数密度
NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

DenK = NumberZ;

% 密度廓线

Prof = 11;
Ns = NumberZ(:,Prof);
% Ns = NumberZ(:,Prof) + NumberZ(:,Prof+1);
plot(Z,Ns,'-b','linewidth',2)
grid on;
Profstr = string(Prof)+" | "+string(TimeShuK(Prof));
title(Profstr)
xlim([70 125]);
% ylim([0 20]);
xlabel('Altitude (km)');
ylabel('Density (cm^{-3})');
legend('K');

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'K Prof ok '

%% Fe

TS = 'Fe Start '

fldstart = ['F:\RawData\ZWDATA\MOHEnew\Fe\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;

FileFe = num_files;

if FileFe>=0

    DenPh = zeros(4096, num_files);
    for j = 1:num_files
        filename = fullfile(folder, files(j).name);
        fid = fopen(filename);
        PhData = textscan(fid, '%f %f','HeaderLines',20);   % 7月22日及之前均为23，后应改为20
        fclose(fid);
        DenPh(:,j) = PhData{1,2}(1:4096);
    end
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
    jNum = FejNum;
    iNum = FeiNum;
    TimeList = cell(1,floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
    end

    TimeXFe = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
    height_num_origin = PhData{1,1}(1:4096);

    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumFe = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumFe(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeFe = height_num_origin(1:floor(4096/iNum))*iNum;
    KM30Fe = size(AltitudeFe(AltitudeFe<30),1)+1;
    KM75Fe = size(AltitudeFe(AltitudeFe<75),1)+1;
    KM80Fe = size(AltitudeFe(AltitudeFe<80),1)+1;
    KM110Fe = size(AltitudeFe(AltitudeFe<110),1)+1;
    KM115Fe = size(AltitudeFe(AltitudeFe<115),1)+1;
    KM120Fe = size(AltitudeFe(AltitudeFe<120),1)+1;
    KM125Fe = size(AltitudeFe(AltitudeFe<125),1)+1;

end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Fe ok '

%% Fe 密度反演

TimeShuFe = TimeXFe';

Noise = mean(PhjSumFe(KM120Fe:KM125Fe,:),1);

ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));
RayScatter = Scatter(2,1);
EffScatter = Scatter(2,2);

[TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,Year,DayNum,17);

Z = AltitudeFe;                       % 高度矩阵
ZR = 30;                            % 参考高度
SigmaRay = RayScatter;              % 瑞利后向散射截面
SigmaEff = EffScatter;              % 有效后向散射截面
NZ = PhjSumFe;                        % 光子数矩阵
NB = Noise;                         % 噪声矩阵
NZR = PhjSumFe(KM30Fe,:);                % 参考高度处光子数
NumberRay = sum(DenRay)*1e-6;       % 参考高度处大气模型数密度
NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

DenFe = NumberZ;

% 密度廓线

Prof = 11;
Ns = NumberZ(:,Prof);
% Ns = NumberZ(:,Prof) + NumberZ(:,Prof+1);
plot(Z,Ns,'-b','linewidth',2)
grid on;
Profstr = string(Prof)+" | "+string(TimeShuFe(Prof));
title(Profstr)
xlim([70 125]);
% ylim([0 3e4]);
xlabel('Altitude (km)');
ylabel('Density (cm^{-3})');
legend('Fe');

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Fe Prof ok '

%% Ni

TS = 'Ni Start '

fldstart = ['F:\RawData\ZWDATA\MOHEnew\Ni\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;

FileNi = num_files;

if FileNi>=0

    DenPh = zeros(4096, num_files);
    for j = 1:num_files
        filename = fullfile(folder, files(j).name);
        fid = fopen(filename);
        PhData = textscan(fid, '%f %f','HeaderLines',20);   % 7月22日及之前均为23，后应改为20
        fclose(fid);
        DenPh(:,j) = PhData{1,2}(1:4096);
    end
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
    jNum = NijNum;
    iNum = NiiNum;
    TimeList = cell(1,floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
    end

    TimeXNi = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
    height_num_origin = PhData{1,1}(1:4096);

    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumNi = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumNi(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeNi = height_num_origin(1:floor(4096/iNum))*iNum;
    KM30Ni = size(AltitudeNi(AltitudeNi<30),1)+1;
    KM75Ni = size(AltitudeNi(AltitudeNi<75),1)+1;
    KM80Ni = size(AltitudeNi(AltitudeNi<80),1)+1;
    KM110Ni = size(AltitudeNi(AltitudeNi<110),1)+1;
    KM115Ni = size(AltitudeNi(AltitudeNi<115),1)+1;
    KM120Ni = size(AltitudeNi(AltitudeNi<120),1)+1;
    KM125Ni = size(AltitudeNi(AltitudeNi<125),1)+1;

end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Ni ok '

%% Ni 密度反演

TimeShuNi = TimeXNi';

Noise = mean(PhjSumNi(KM120Ni:KM125Ni,:),1);

ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));
RayScatter = Scatter(1,1);
EffScatter = Scatter(1,2);

[TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,Year,DayNum,17);

Z = AltitudeNi;                       % 高度矩阵
ZR = 30;                            % 参考高度
SigmaRay = RayScatter;              % 瑞利后向散射截面
SigmaEff = EffScatter;              % 有效后向散射截面
NZ = PhjSumNi;                        % 光子数矩阵
NB = Noise;                         % 噪声矩阵
NZR = PhjSumNi(KM30Ni,:);                % 参考高度处光子数
NumberRay = sum(DenRay)*1e-6;       % 参考高度处大气模型数密度
NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

DenNi= NumberZ;

% 密度廓线

Prof = 8;
Ns = NumberZ(:,Prof);
plot(Z,Ns,'-b','linewidth',1.5)
grid on;
Profstr = string(Prof)+" | "+string(TimeShuNi(Prof));
title(Profstr)
xlim([70 125]);
ylim([0 500]);
xlabel('Altitude (km)');
ylabel('Density (cm^{-3})');
legend('Ni');

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Ni Prof ok '


%% Ca

TS = 'Ca Start '

fldstart = ['F:\RawData\ZWDATA\MOHEnew\Ca\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;

FileCa = num_files;

if FileCa>=0

    DenPh1 = zeros(4096, num_files);
    folder = [fldstart,'CH1\'];
    files = dir(fullfile(folder, '*.txt'));
    for j = 1:num_files
        filename = fullfile(folder,files(j).name);
        fid = fopen(filename);
        PhData1 = textscan(fid, '%f %f','HeaderLines',20);
        fclose(fid);
        DenPh1(:,j) = PhData1{1,2}(1:4096,:);
    end

    DenPh2 = zeros(4096, num_files);
    folder = [fldstart,'CH2\'];
    files = dir(fullfile(folder, '*.txt'));
    for j = 1:num_files
        filename = fullfile(folder,files(j).name);
        fid = fopen(filename);
        PhData2 = textscan(fid, '%f %f','HeaderLines',20);
        fclose(fid);
        DenPh2(:,j) = PhData2{1,2}(1:4096,:);
    end

    DenPh3 = zeros(4096, num_files);
    folder = [fldstart,'CH3\'];
    files = dir(fullfile(folder, '*.txt'));
    for j = 1:num_files
        filename = fullfile(folder,files(j).name);
        fid = fopen(filename);
        PhData3 = textscan(fid, '%f %f','HeaderLines',20);
        fclose(fid);
        DenPh3(:,j) = PhData3{1,2}(1:4096,:);
    end

    DenPh4 = zeros(4096, num_files);
    folder = [fldstart,'CH4\'];
    files = dir(fullfile(folder, '*.txt'));
    for j = 1:num_files
        filename = fullfile(folder,files(j).name);
        fid = fopen(filename);
        PhData4 = textscan(fid, '%f %f','HeaderLines',20);
        fclose(fid);
        DenPh4(:,j) = PhData4{1,2}(1:4096,:);
    end

    DenPh = DenPh1+DenPh2+DenPh3+DenPh4;

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
    jNum = CajNum;
    iNum = CaiNum;
    TimeList = cell(1,floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
    end

    TimeXCa = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
    height_num_origin = PhData1{1,1}(1:4096);

    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumCa = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumCa(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeCa = height_num_origin(1:floor(4096/iNum))*iNum;
    KM30Ca = size(AltitudeCa(AltitudeCa<30),1)+1;
    KM75Ca = size(AltitudeCa(AltitudeCa<75),1)+1;
    KM80Ca = size(AltitudeCa(AltitudeCa<80),1)+1;
    KM110Ca = size(AltitudeCa(AltitudeCa<110),1)+1;
    KM115Ca = size(AltitudeCa(AltitudeCa<115),1)+1;
    KM120Ca = size(AltitudeCa(AltitudeCa<120),1)+1;
    KM125Ca = size(AltitudeCa(AltitudeCa<125),1)+1;

end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Ca ok '

%% Ca 密度反演

TimeShuCa = TimeXCa';

Noise = mean(PhjSumCa(KM120Ca:KM125Ca,:),1);

ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));
RayScatter = Scatter(4,1);
EffScatter = Scatter(4,2);

[TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,Year,DayNum,17);

Z = AltitudeCa;                       % 高度矩阵
ZR = 30;                            % 参考高度
SigmaRay = RayScatter;              % 瑞利后向散射截面
SigmaEff = EffScatter;              % 有效后向散射截面
NZ = PhjSumCa;                        % 光子数矩阵
NB = Noise;                         % 噪声矩阵
NZR = PhjSumCa(KM30Ca,:);                % 参考高度处光子数
NumberRay = sum(DenRay)*1e-6;       % 参考高度处大气模型数密度
NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

DenCa = NumberZ;

% 密度廓线

Prof = 5;
Ns = NumberZ(:,Prof);
Ns = NumberZ(:,Prof) + NumberZ(:,Prof+1);
plot(Z,Ns,'-','linewidth',1.5)
grid on;
Profstr = string(Prof)+" | "+string(TimeShuCa(Prof));
title(Profstr)
xlim([70 125]);
% ylim([0 80]);
xlabel('Altitude (km)');
ylabel('Density (cm^{-3})');
legend('Ca');

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Ca Prof ok '


%% Cap

TS = 'Cap Start '

fldstart = ['F:\RawData\ZWDATA\MOHEnew\Cap\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;

FileCap = num_files;

if FileCap>=0

    DenPh1 = zeros(4096, num_files);
    folder = [fldstart,'CH1\'];
    files = dir(fullfile(folder, '*.txt'));
    for j = 1:num_files
        filename = fullfile(folder,files(j).name);
        fid = fopen(filename);
        PhData1 = textscan(fid, '%f %f','HeaderLines',20);
        fclose(fid);
        DenPh1(:,j) = PhData1{1,2}(1:4096,:);
    end

    DenPh2 = zeros(4096, num_files);
    folder = [fldstart,'CH2\'];
    files = dir(fullfile(folder, '*.txt'));
    for j = 1:num_files
        filename = fullfile(folder,files(j).name);
        fid = fopen(filename);
        PhData2 = textscan(fid, '%f %f','HeaderLines',20);
        fclose(fid);
        DenPh2(:,j) = PhData2{1,2}(1:4096,:);
    end

    DenPh3 = zeros(4096, num_files);
    folder = [fldstart,'CH3\'];
    files = dir(fullfile(folder, '*.txt'));
    for j = 1:num_files
        filename = fullfile(folder,files(j).name);
        fid = fopen(filename);
        PhData3 = textscan(fid, '%f %f','HeaderLines',20);
        fclose(fid);
        DenPh3(:,j) = PhData3{1,2}(1:4096,:);
    end

    DenPh4 = zeros(4096, num_files);
    folder = [fldstart,'CH4\'];
    files = dir(fullfile(folder, '*.txt'));
    for j = 1:num_files
        filename = fullfile(folder,files(j).name);
        fid = fopen(filename);
        PhData4 = textscan(fid, '%f %f','HeaderLines',20);
        fclose(fid);
        DenPh4(:,j) = PhData4{1,2}(1:4096,:);
    end

    DenPh = DenPh1+DenPh2+DenPh3+DenPh4;

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
    jNum = CapjNum;
    iNum = CapiNum;
    TimeList = cell(1,floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
    end

    TimeXCap = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
    height_num_origin = PhData1{1,1}(1:4096);

    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumCap = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumCap(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeCap = height_num_origin(1:floor(4096/iNum))*iNum;
    KM30Cap = size(AltitudeCap(AltitudeCap<30),1)+1;
    KM75Cap = size(AltitudeCap(AltitudeCap<75),1)+1;
    KM80Cap = size(AltitudeCap(AltitudeCap<80),1)+1;
    KM110Cap = size(AltitudeCap(AltitudeCap<110),1)+1;
    KM115Cap = size(AltitudeCap(AltitudeCap<115),1)+1;
    KM120Cap = size(AltitudeCap(AltitudeCap<120),1)+1;
    KM125Cap = size(AltitudeCap(AltitudeCap<125),1)+1;

end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Cap ok '

%% Cap 密度反演

TimeShuCap = TimeXCap';

Noise = mean(PhjSumCap(KM120Cap:KM125Cap,:),1);

ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));
RayScatter = Scatter(3,1);
EffScatter = Scatter(3,2);

[TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,Year,DayNum,17);

Z = AltitudeCap;                       % 高度矩阵
ZR = 30;                            % 参考高度
SigmaRay = RayScatter;              % 瑞利后向散射截面
SigmaEff = EffScatter;              % 有效后向散射截面
NZ = PhjSumCap;                        % 光子数矩阵
NB = Noise;                         % 噪声矩阵
NZR = PhjSumCap(KM30Cap,:);                % 参考高度处光子数
NumberRay = sum(DenRay)*1e-6;       % 参考高度处大气模型数密度
NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

DenCap = NumberZ;

% 密度廓线

Prof = 4;
Ns = NumberZ(:,Prof);
% Ns = NumberZ(:,Prof) + NumberZ(:,Prof+1);
plot(Z,Ns,'-b','linewidth',1.5)
grid on;
Profstr = string(Prof)+" | "+string(TimeShuCap(Prof));
title(Profstr)
xlim([70 125]);
% ylim([0 700]);
xlabel('Altitude (km)');
ylabel('Density (cm^{-3})');
legend('Ca ion');

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Cap Prof ok '







