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
CHN = 'Ca';
XlsName = char("F:\TrailIMG\" + CHN + ".xlsx");
SDate = RunDate(5:8);
SheetName = char('S'+string(SDate));
Table = readtable(XlsName,'sheet',SheetName);

CaList = table2array(Table(:,10));
CapList = table2array(Table(:,11));

CajNum = 1; CaiNum = 1;CaDex = 1/1;
CapjNum = 1; CapiNum = 1;CapDex = 1/1;

List = CaList;

PF = List(~isnan(List));
Ttk = 1/1440;
TtkUp = 5/1440;
TtkDp = 5/1440;

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Excel ok '


%% Ca

TS = 'Ca start '

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


%% Cap

TS = 'Cap start '

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



%% Plot Time

TimeSum = TimeXCa;

for jt = PF'

start_time = datenum(TimeSum(:,jt))-TtkUp;
end_time = datenum(TimeSum(:,jt))+TtkDp;

fC = figure('name','全队列','position',[10 50 1520 820]);


subplot(2,1,1);

if FileCa >= 0

    PhX = PhjSumCa(KM75Ca:KM115Ca,jt);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumCa(KM75Ca:KM115Ca,:).*((1/1440)./PhMax).*CaDex;
    Ca = PhMin;
    TimeD = datenum(TimeXCa);
    for jn = 1:size(PhMin,2)
        Ca(:,jn) = datenum(TimeXCa(:,jn)) + PhMin(:,jn);
    end
    plot(Ca,AltitudeCa(KM75Ca:KM115Ca,:),'linewidth',1.5)
    grid on;
    ylim([75 115])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits','keepticks');
    CaStr  = "Ca " + RunDate;
    title(CaStr)
    
end

subplot(2,1,2);

if FileCap >= 0

    PhX = PhjSumCap(KM75Cap:KM115Cap,:);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumCap(KM75Cap:KM115Cap,:).*((1/1440)./PhMax).*CapDex;
    Cap = PhMin;
    TimeD = datenum(TimeXCap);
    for jn = 1:size(PhMin,2)
        Cap(:,jn) = datenum(TimeXCap(:,jn)) + PhMin(:,jn);
    end
    plot(Cap,AltitudeCap(KM75Cap:KM115Cap,:),'linewidth',1.5)
    grid on;
    ylim([75 115])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits','keepticks');
    CapStr  = "Cap " + RunDate;
    title(CapStr)
    
end


    % 指定保存路径和文件名
    figPath = "F:\TrailIMG\SumCaCapTrailFind\"+ReadDate+"\";
    if ~exist(figPath, 'dir')
        mkdir(figPath);
    end
    figV = datestr(TimeSum(:,jt),'yyyymmdd-HH-MM-SS')+"-Ca-"+string(jt);
    print(fC, '-dpng', '-r300', figPath+figV);
    close all;

end

TS = 'Figures ok '


fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*1200*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*900*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.5;y = sin(2*pi*800*t);sound(y, fs);


