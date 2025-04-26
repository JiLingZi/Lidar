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

RunDate = '20231027';
ReadDate = RunDate;
Year = str2double(char(RunDate(1:4)));
Month = str2double(char(RunDate(5:6)));
Day = str2double(char(RunDate(7:8)));
DayNum = Month*30-30+Day;

CajNum = 1; CaiNum = 1;CaDex = 1/1;
CapjNum = 1; CapiNum = 1;CapDex = 1/1;


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

% 密度廓线

Prof = 305;
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

% 密度廓线

Prof = 219;
Ns = NumberZ(:,Prof);
% Ns = NumberZ(:,Prof) + NumberZ(:,Prof+1);
plot(Z,Ns,'-','linewidth',1.5)
grid on;
Profstr = string(Prof)+" | "+string(TimeShuCap(Prof));
title(Profstr)
xlim([70 125]);
% ylim([0 80]);
xlabel('Altitude (km)');
ylabel('Density (cm^{-3})');
legend('Ca ion');

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Cap Prof ok '







