% Ca密度

clc;clear all;

DateList = {'20240114';'20240115';'20240117';'20240118';'20240119';...
            '20240120';'20240121';'20240124';'20240125';'20240128';...
            '20240129';'20240131';'20240311';'20240510';'20240511';...
            '20240512';'20240514';'20240515';'20240518';'20240519';...
            '20240524';'20240525';'20240526';'20240528';'20240529';...
            '20240530';'20240531';'20240601';'20240629';'20240707';};

for jd = 1:size(DateList,1)
% 读取光子数矩阵
ReadDate = char(DateList(jd));
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Ca\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1
apple(jd,:) = num_files;


% 为了统一文件数，减一
% 减一是因为那个破机器动不动少一个文件，导致不晓得哪个通道文件数不同
% 格式不要再改了，探测高度上限也不要再改了，否则实现自动化反演困难重重


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

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = '数据读取OK了'

%% 合并时间文件数
% jNum = 15;
iNum = 1;
jNum = 1;
TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

TimeX = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
TimeY = 1:size(TimeX,2);
hours = hour(TimeX);
minutes = minute(TimeX);
seconds = second(TimeX);
Time = (((hours(1,:)+8)*3600+minutes(1,:)*60+seconds(1,:))+((((hours(end,:)+8)*3600+minutes(end,:)*60+seconds(end,:))-((hours(1,:)+8)*3600+minutes(1,:)*60)+seconds(1,:)).*(TimeY./TimeY(end))))./3600;



% 获取原始高度矩阵
height_num_origin = PhData1{1,1};

% 合并行，3行合并，96米
PhiSum = zeros(floor(4096/iNum),num_files);
for i = 1:floor(4096/iNum)
    PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
end

% 合并列，15列合并，15分钟
PhjSum = zeros(floor(4096/iNum),floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    PhjSum(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
end

% 噪声120-150，计算光子数，高度
Altitude = height_num_origin(1:floor(4096/iNum))*iNum;
KM30 = size(Altitude(Altitude<30),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM125 = size(Altitude(Altitude<125),1)+1;


%% Ca队列检查

for jt = 1:ceil(num_files/24)
    
    fC = figure('name','Ca队列','position',[10 50 2400 1200]);
    subplot(4,6,1);
    if 24*(jt-1)+1 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+1))
        xlim([75 115])
        t1str = string(24*(jt-1)+1)+" | "+datestr(TimeX(:,24*(jt-1)+1),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,2);
    if 24*(jt-1)+2 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+2))
        xlim([75 115])
        t1str = string(24*(jt-1)+2)+" | "+datestr(TimeX(:,24*(jt-1)+2),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,3);
    if 24*(jt-1)+3 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+3))
        xlim([75 115])
        t1str = string(24*(jt-1)+3)+" | "+datestr(TimeX(:,24*(jt-1)+3),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,4);
    if 24*(jt-1)+4 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+4))
        xlim([75 115])
        t1str = string(24*(jt-1)+4)+" | "+datestr(TimeX(:,24*(jt-1)+4),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,5);
    if 24*(jt-1)+5 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+5))
        xlim([75 115])
        t1str = string(24*(jt-1)+5)+" | "+datestr(TimeX(:,24*(jt-1)+5),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,6);
    if 24*(jt-1)+6 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+6))
        xlim([75 115])
        t1str = string(24*(jt-1)+6)+" | "+datestr(TimeX(:,24*(jt-1)+6),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,7);
    if 24*(jt-1)+7 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+7))
        xlim([75 115])
        t1str = string(24*(jt-1)+7)+" | "+datestr(TimeX(:,24*(jt-1)+7),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,8);
    if 24*(jt-1)+8 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+8))
        xlim([75 115])
        t1str = string(24*(jt-1)+8)+" | "+datestr(TimeX(:,24*(jt-1)+8),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,9);
    if 24*(jt-1)+9 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+9))
        xlim([75 115])
        t1str = string(24*(jt-1)+9)+" | "+datestr(TimeX(:,24*(jt-1)+9),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,10);
    if 24*(jt-1)+10 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+10))
        xlim([75 115])
        t1str = string(24*(jt-1)+10)+" | "+datestr(TimeX(:,24*(jt-1)+10),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,11);
    if 24*(jt-1)+11 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+11))
        xlim([75 115])
        t1str = string(24*(jt-1)+11)+" | "+datestr(TimeX(:,24*(jt-1)+11),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,12);
    if 24*(jt-1)+12 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+12))
        xlim([75 115])
        t1str = string(24*(jt-1)+12)+" | "+datestr(TimeX(:,24*(jt-1)+12),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,13);
    if 24*(jt-1)+13 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+13))
        xlim([75 115])
        t1str = string(24*(jt-1)+13)+" | "+datestr(TimeX(:,24*(jt-1)+13),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,14);
    if 24*(jt-1)+14 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+14))
        xlim([75 115])
        t1str = string(24*(jt-1)+14)+" | "+datestr(TimeX(:,24*(jt-1)+14),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,15);
    if 24*(jt-1)+15 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+15))
        xlim([75 115])
        t1str = string(24*(jt-1)+15)+" | "+datestr(TimeX(:,24*(jt-1)+15),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,16);
    if 24*(jt-1)+16 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+16))
        xlim([75 115])
        t1str = string(24*(jt-1)+16)+" | "+datestr(TimeX(:,24*(jt-1)+16),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,17);
    if 24*(jt-1)+17 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+17))
        xlim([75 115])
        t1str = string(24*(jt-1)+17)+" | "+datestr(TimeX(:,24*(jt-1)+17),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,18);
    if 24*(jt-1)+18 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+18))
        xlim([75 115])
        t1str = string(24*(jt-1)+18)+" | "+datestr(TimeX(:,24*(jt-1)+18),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,19);
    if 24*(jt-1)+19 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+19))
        xlim([75 115])
        t1str = string(24*(jt-1)+19)+" | "+datestr(TimeX(:,24*(jt-1)+19),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,20);
    if 24*(jt-1)+20 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+20))
        xlim([75 115])
        t1str = string(24*(jt-1)+20)+" | "+datestr(TimeX(:,24*(jt-1)+20),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,21);
    if 24*(jt-1)+21 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+21))
        xlim([75 115])
        t1str = string(24*(jt-1)+21)+" | "+datestr(TimeX(:,24*(jt-1)+21),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,22);
    if 24*(jt-1)+22 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+22))
        xlim([75 115])
        t1str = string(24*(jt-1)+22)+" | "+datestr(TimeX(:,24*(jt-1)+22),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,23);
    if 24*(jt-1)+23 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+23))
        xlim([75 115])
        t1str = string(24*(jt-1)+23)+" | "+datestr(TimeX(:,24*(jt-1)+23),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,24);
    if 24*(jt-1)+24 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+24))
        xlim([75 115])
        t1str = string(24*(jt-1)+24)+" | "+datestr(TimeX(:,24*(jt-1)+24),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    
    % 指定保存路径和文件名
    figPath = "D:\KYBF\TimingTrailIMG\Ca\"+ReadDate+"\";
    if ~exist(figPath, 'dir')
        mkdir(figPath);
    end
    figV = datestr(TimeX(:,jt),'yyyymmdd')+"-Ca-"+string(jt);
    print(fC, '-dpng', '-r300', figPath+figV);
    close all;
    
end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*1200*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*900*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.5;y = sin(2*pi*800*t);sound(y, fs);

end

% CaIon密度

clc;clear all;

DateList = {'20240114';'20240115';'20240117';'20240118';'20240119';...
            '20240120';'20240121';'20240124';'20240125';'20240128';...
            '20240129';'20240131';'20240311';'20240510';'20240511';...
            '20240512';'20240515';'20240519';...
            '20240524';'20240525';'20240526';'20240528';'20240529';...
            '20240530';'20240531';'20240601';};

for jd = 1:size(DateList,1)
% 读取光子数矩阵
ReadDate = char(DateList(jd));
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Cap\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1
apple(jd,:) = num_files;


% 为了统一文件数，减一
% 减一是因为那个破机器动不动少一个文件，导致不晓得哪个通道文件数不同
% 格式不要再改了，探测高度上限也不要再改了，否则实现自动化反演困难重重

bin = 26432;


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

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);


%% 合并时间文件数
% jNum = 15;
iNum = 1;
jNum = 1;
TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

TimeX = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
TimeY = 1:size(TimeX,2);
hours = hour(TimeX);
minutes = minute(TimeX);
seconds = second(TimeX);
Time = (((hours(1,:)+8)*3600+minutes(1,:)*60+seconds(1,:))+((((hours(end,:)+8)*3600+minutes(end,:)*60+seconds(end,:))-((hours(1,:)+8)*3600+minutes(1,:)*60)+seconds(1,:)).*(TimeY./TimeY(end))))./3600;



% 获取原始高度矩阵
height_num_origin = PhData1{1,1};

% 合并行，3行合并，96米
PhiSum = zeros(floor(4096/iNum),num_files);
for i = 1:floor(4096/iNum)
    PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
end

% 合并列，15列合并，15分钟
PhjSum = zeros(floor(4096/iNum),floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    PhjSum(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
end

% 噪声120-150，计算光子数，高度
Altitude = height_num_origin(1:floor(4096/iNum))*iNum;
KM30 = size(Altitude(Altitude<30),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM125 = size(Altitude(Altitude<125),1)+1;
KM350 = size(Altitude(Altitude<350),1)+1;
KM400 = size(Altitude(Altitude<400),1)+1;

%% Cap队列检查

for jt = 1:ceil(num_files/24)
    
    fC = figure('name','Cap队列','position',[10 50 2400 1200]);
    subplot(4,6,1);
    if 24*(jt-1)+1 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+1))
        xlim([75 115])
        t1str = string(24*(jt-1)+1)+" | "+datestr(TimeX(:,24*(jt-1)+1),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,2);
    if 24*(jt-1)+2 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+2))
        xlim([75 115])
        t1str = string(24*(jt-1)+2)+" | "+datestr(TimeX(:,24*(jt-1)+2),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,3);
    if 24*(jt-1)+3 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+3))
        xlim([75 115])
        t1str = string(24*(jt-1)+3)+" | "+datestr(TimeX(:,24*(jt-1)+3),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,4);
    if 24*(jt-1)+4 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+4))
        xlim([75 115])
        t1str = string(24*(jt-1)+4)+" | "+datestr(TimeX(:,24*(jt-1)+4),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,5);
    if 24*(jt-1)+5 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+5))
        xlim([75 115])
        t1str = string(24*(jt-1)+5)+" | "+datestr(TimeX(:,24*(jt-1)+5),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,6);
    if 24*(jt-1)+6 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+6))
        xlim([75 115])
        t1str = string(24*(jt-1)+6)+" | "+datestr(TimeX(:,24*(jt-1)+6),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,7);
    if 24*(jt-1)+7 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+7))
        xlim([75 115])
        t1str = string(24*(jt-1)+7)+" | "+datestr(TimeX(:,24*(jt-1)+7),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,8);
    if 24*(jt-1)+8 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+8))
        xlim([75 115])
        t1str = string(24*(jt-1)+8)+" | "+datestr(TimeX(:,24*(jt-1)+8),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,9);
    if 24*(jt-1)+9 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+9))
        xlim([75 115])
        t1str = string(24*(jt-1)+9)+" | "+datestr(TimeX(:,24*(jt-1)+9),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,10);
    if 24*(jt-1)+10 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+10))
        xlim([75 115])
        t1str = string(24*(jt-1)+10)+" | "+datestr(TimeX(:,24*(jt-1)+10),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,11);
    if 24*(jt-1)+11 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+11))
        xlim([75 115])
        t1str = string(24*(jt-1)+11)+" | "+datestr(TimeX(:,24*(jt-1)+11),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,12);
    if 24*(jt-1)+12 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+12))
        xlim([75 115])
        t1str = string(24*(jt-1)+12)+" | "+datestr(TimeX(:,24*(jt-1)+12),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,13);
    if 24*(jt-1)+13 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+13))
        xlim([75 115])
        t1str = string(24*(jt-1)+13)+" | "+datestr(TimeX(:,24*(jt-1)+13),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,14);
    if 24*(jt-1)+14 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+14))
        xlim([75 115])
        t1str = string(24*(jt-1)+14)+" | "+datestr(TimeX(:,24*(jt-1)+14),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,15);
    if 24*(jt-1)+15 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+15))
        xlim([75 115])
        t1str = string(24*(jt-1)+15)+" | "+datestr(TimeX(:,24*(jt-1)+15),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,16);
    if 24*(jt-1)+16 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+16))
        xlim([75 115])
        t1str = string(24*(jt-1)+16)+" | "+datestr(TimeX(:,24*(jt-1)+16),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,17);
    if 24*(jt-1)+17 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+17))
        xlim([75 115])
        t1str = string(24*(jt-1)+17)+" | "+datestr(TimeX(:,24*(jt-1)+17),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,18);
    if 24*(jt-1)+18 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+18))
        xlim([75 115])
        t1str = string(24*(jt-1)+18)+" | "+datestr(TimeX(:,24*(jt-1)+18),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,19);
    if 24*(jt-1)+19 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+19))
        xlim([75 115])
        t1str = string(24*(jt-1)+19)+" | "+datestr(TimeX(:,24*(jt-1)+19),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,20);
    if 24*(jt-1)+20 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+20))
        xlim([75 115])
        t1str = string(24*(jt-1)+20)+" | "+datestr(TimeX(:,24*(jt-1)+20),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,21);
    if 24*(jt-1)+21 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+21))
        xlim([75 115])
        t1str = string(24*(jt-1)+21)+" | "+datestr(TimeX(:,24*(jt-1)+21),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,22);
    if 24*(jt-1)+22 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+22))
        xlim([75 115])
        t1str = string(24*(jt-1)+22)+" | "+datestr(TimeX(:,24*(jt-1)+22),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,23);
    if 24*(jt-1)+23 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+23))
        xlim([75 115])
        t1str = string(24*(jt-1)+23)+" | "+datestr(TimeX(:,24*(jt-1)+23),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    subplot(4,6,24);
    if 24*(jt-1)+24 <= num_files
        plot(Altitude,PhjSum(:,24*(jt-1)+24))
        xlim([75 115])
        t1str = string(24*(jt-1)+24)+" | "+datestr(TimeX(:,24*(jt-1)+24),'yyyymmdd|HH:MM:SS');
        title(t1str)
        grid on;
    end
    
    % 指定保存路径和文件名
    figPath = "D:\KYBF\TimingTrailIMG\Cap\"+ReadDate+"\";
    if ~exist(figPath, 'dir')
        mkdir(figPath);
    end
    figV = datestr(TimeX(:,jt),'yyyymmdd')+"-Cap-"+string(jt);
    print(fC, '-dpng', '-r300', figPath+figV);
    close all;
    
end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*1200*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*900*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.5;y = sin(2*pi*800*t);sound(y, fs);



end

