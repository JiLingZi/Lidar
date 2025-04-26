% Na密度

clc;clear all;

DateList = {'20240810';'20240811';'20240826';'20240828';'20240830';...
            '20240901';'20240902';'20240904';'20240906';'20240910';};

for jd = 1:size(DateList,1)
% 读取光子数矩阵
ReadDate = char(DateList(jd));
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Na\',ReadDate,'\'];
folder = [fldstart,'Na\'];
files = dir(fullfile(folder, '*.dat'));
num_files = length(files)-1;
num_files = 3*num_files

% 为了统一文件数，减一
% 减一是因为那个破机器动不动少一个文件，导致不晓得哪个通道文件数不同
% 格式不要再改了，探测高度上限也不要再改了，否则实现自动化反演困难重重

V_F_0 = zeros(8192, num_files/3);
V_F_R = zeros(8192, num_files/3);
V_F_L = zeros(8192, num_files/3);
N_F_0 = zeros(8192, num_files/3);
N_F_R = zeros(8192, num_files/3);
N_F_L = zeros(8192, num_files/3);
E_F_0 = zeros(8192, num_files/3);
E_F_R = zeros(8192, num_files/3);
E_F_L = zeros(8192, num_files/3);
for j = 1:(num_files/3)
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

% DenPh = zeros(8192,num_files*3);
% for j = 1:num_files
%     DenPh(:,(j-1)*3+1) = V_F_0(:,j);
%     DenPh(:,(j-1)*3+2) = V_F_R(:,j);
%     DenPh(:,(j-1)*3+3) = V_F_L(:,j);
% end

DenPhV = zeros(8192, num_files);
DenPhN = zeros(8192, num_files);
DenPhE = zeros(8192, num_files);

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

% DenPhV = V_F_0 + V_F_R + V_F_L;
% % DenPhV = V_F_0 ;
% DenPhN = N_F_0 + N_F_R + N_F_L;
% 
% DenPhE = E_F_0 + E_F_R + E_F_L;




% 读取时间序列
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

% 合并时间文件数
jNum = 1;
iNum = 1;
% jNum = 1;
TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
	TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

TimeX = datetime(TimeList, 'InputFormat', ':yyyy-MM-dd''T''HH:mm:ss.SSS');
TimeY = 1:size(TimeX,2);
hours = hour(TimeX);
minutes = minute(TimeX);
TimeLong = (((hours(1,:)+8)*3600+minutes(1,:)*60)+((((hours(end,:)+8)*3600+minutes(end,:)*60)-((hours(1,:)+8)*3600+minutes(1,:)*60)).*(TimeY./TimeY(end))))./3600;
Timexx = TimeLong;

TimeL= 63;
TimeR= 78;
xTimeL = Timexx(TimeL);
xTimeR = Timexx(TimeR);
TimeP = Timexx(TimeL:TimeR);
Time = Timexx;

RunRead = '数据读取OK了'

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);


%% V 获取原始高度矩阵
height_num_origin = Na_data(1:8192,1);

% 合并行，3行合并，96米
PhiSum = zeros(floor(8192/iNum),num_files);
for i = 1:floor(8192/iNum)
    PhiSum(i,:) = sum(DenPhV(iNum*(i-1)+1:iNum*i,:),1);
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

% % 归一化光子数检查
% CheckH = 30;% 光子数检查归一化高度，KM
% KMCheck = size(Altitude(Altitude<CheckH),1);
% CheckPh = PhjSum./PhjSum(KMCheck,:);
% figure('Name','V光子数检查','position',[10 100 1350 500])
% for j = 1:size(CheckPh,2)
%     plot(CheckPh(:,j)+j,Altitude,'-','linewidth',1.5);
%     hold on;
% end
% set(gca,'FontSize',12,'FontName','Times New Roman');
% title('Fe Density Profiles 20231102 @Mohe');
% xlabel('Local Time');
% ylabel('Altitude (km)');
% xlim([-1 11]);
% ylim([75 115]);
% % xlim([-1 41]);
% % ylim([75 115]);

%% V队列检查

for jt = 1:ceil(num_files/24)
    
    fC = figure('name','V队列','position',[10 50 2400 1200]);
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
    figPath = "D:\KYBF\TimingTrailIMG\Na\"+ReadDate+"\";
    if ~exist(figPath, 'dir')
        mkdir(figPath);
    end
    figV = datestr(TimeX(:,jt),'yyyymmdd')+"-V-"+string(jt);
    print(fC, '-dpng', '-r300', figPath+figV);
    close all;
    
end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*600*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);



%% N 合并行，3行合并，96米
PhiSum = zeros(floor(8192/iNum),num_files);
for i = 1:floor(8192/iNum)
    PhiSum(i,:) = sum(DenPhN(iNum*(i-1)+1:iNum*i,:),1);
end

% 合并列，15列合并，15分钟
PhjSum = zeros(floor(8192/iNum),floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    PhjSum(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
end

% 噪声120-150，计算光子数，高度
Altitude = height_num_origin(1:floor(8192/iNum))*iNum*(sqrt(3)/2);
KM30 = size(Altitude(Altitude<30),1)+1;
KM75 = size(Altitude(Altitude<75),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
KM115 = size(Altitude(Altitude<115),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM125 = size(Altitude(Altitude<125),1)+1;

% % 归一化光子数检查
% CheckH = 30;% 光子数检查归一化高度，KM
% KMCheck = size(Altitude(Altitude<CheckH),1);
% CheckPh = PhjSum./PhjSum(KMCheck,:);
% figure('Name','N光子数检查','position',[10 100 1350 500])
% for j = 1:size(CheckPh,2)
%     plot(CheckPh(:,j)+j,Altitude,'-','linewidth',1.5);
%     hold on;
% end
% set(gca,'FontSize',12,'FontName','Times New Roman');
% title('Fe Density Profiles 20231102 @Mohe');
% xlabel('Local Time');
% ylabel('Altitude (km)');
% xlim([-1 31]);
% ylim([75 115]);
% % xlim([-1 41]);
% % ylim([75 115]);


% N队列检查

for jt = 1:ceil(num_files/24)
    
    fC = figure('name','V队列','position',[10 50 2400 1200]);
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
    figPath = "D:\KYBF\TimingTrailIMG\Na\"+ReadDate+"\";
    if ~exist(figPath, 'dir')
        mkdir(figPath);
    end
    figV = datestr(TimeX(:,jt),'yyyymmdd')+"-N-"+string(jt);
    print(fC, '-dpng', '-r300', figPath+figV);
    close all;
    
end


fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*600*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);




%% E 合并行，3行合并，96米
PhiSum = zeros(floor(8192/iNum),num_files);
for i = 1:floor(8192/iNum)
    PhiSum(i,:) = sum(DenPhE(iNum*(i-1)+1:iNum*i,:),1);
end

% 合并列，15列合并，15分钟
PhjSum = zeros(floor(8192/iNum),floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    PhjSum(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
end

% 噪声120-150，计算光子数，高度
Altitude = height_num_origin(1:floor(8192/iNum))*iNum*(sqrt(3)/2);
KM30 = size(Altitude(Altitude<30),1)+1;
KM75 = size(Altitude(Altitude<75),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
KM115 = size(Altitude(Altitude<115),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM125 = size(Altitude(Altitude<125),1)+1;

% % 归一化光子数检查
% CheckH = 30;% 光子数检查归一化高度，KM
% KMCheck = size(Altitude(Altitude<CheckH),1);
% CheckPh = PhjSum./PhjSum(KMCheck,:);
% figure('Name','E光子数检查','position',[10 100 1350 500])
% for j = 1:size(CheckPh,2)
%     plot(CheckPh(:,j)+j,Altitude,'-','linewidth',1.5);
%     hold on;
% end
% set(gca,'FontSize',12,'FontName','Times New Roman');
% title('Fe Density Profiles 20231102 @Mohe');
% xlabel('Local Time');
% ylabel('Altitude (km)');
% xlim([-1 31]);
% ylim([75 115]);
% % xlim([-1 41]);
% % ylim([75 115]);



% E队列检查

for jt = 1:ceil(num_files/24)
    
    fC = figure('name','V队列','position',[10 50 2400 1200]);
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
    figPath = "D:\KYBF\TimingTrailIMG\Na\"+ReadDate+"\";
    if ~exist(figPath, 'dir')
        mkdir(figPath);
    end
    figV = datestr(TimeX(:,jt),'yyyymmdd')+"-E-"+string(jt);
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
%% 
% figure('Name','光子数Time检查','position',[300 300 800 250])
% TXcal = (Time(end)-Time(1))/(size(Time,2)-1);
% TXcheck = CheckPh;
% for j = 1:size(CheckPh,2)
%     TXcheck(:,j) = (CheckPh(:,j)*0.5+j-1).*TXcal+Time(1);
% end
% plot(TXcheck,Altitude,'-','linewidth',1.5);
% hold on;
% title('Na East 20231102')
% TX1 = double(16+(53/60))+8;
% TX2 = double(16+(60/60))+8;
% % xticks = Time;
% xlim([TX1 TX2]);
% ylim([80 100]);
% set(gca,'FontSize',12);
% xlabel('Time (LT)');
% ylabel('Altitude (km)');


fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*1200*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*900*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.5;y = sin(2*pi*800*t);sound(y, fs);


