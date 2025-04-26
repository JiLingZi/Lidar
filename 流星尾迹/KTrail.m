% K密度

clc;clear all;

DateList = {'20240805';'20240811';'20240826';'20240828';'20240830';...
            '20240901';'20240902';'20240904';'20240906';'20240910';};

for jd = 1:size(DateList,1)
% 读取光子数矩阵
ReadDate = char(DateList(jd));
fldstart = ['F:\RawData\ZWDATA\MOHEnew\K532\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1
apple(jd,:) = num_files;

% 为了统一文件数，减一
% 减一是因为那个破机器动不动少一个文件，导致不晓得哪个通道文件数不同
% 格式不要再改了，探测高度上限也不要再改了，否则实现自动化反演困难重重



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

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);


%% 合并时间文件数
% jNum = 30;
jNum = 1;
iNum = 1;
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
KM80 = size(Altitude(Altitude<80),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM125 = size(Altitude(Altitude<125),1)+1;

%% K队列检查

for jt = 1:ceil(num_files/24)
    
    fC = figure('name','K队列','position',[10 50 2400 1200]);
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
    figPath = "D:\KYBF\TimingTrailIMG\K\"+ReadDate+"\";
    if ~exist(figPath, 'dir')
        mkdir(figPath);
    end
    figV = datestr(TimeX(:,jt),'yyyymmdd')+"-K-"+string(jt);
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
% figure('name','检查','position',[600 50 350 800])
% pf = 210;
% plot(PhjSum(:,pf), Altitude, '-k', 'linewidth',1.5)
% set(gca, 'XScale', 'log');
% ylim([30 125])
% grid on
% ttstr =char(TimeList(pf));
% title({
% ['(a) K']
% [ttstr]
% })
% set(gca,'FontSize',12.5);
% % ylim([3000 7000])
% xlabel('Photon Counts')
% ylabel('Altitude (km)')

%% 归一化光子数检查
% CheckH = 60;% 光子数检查归一化高度，KM
% KMCheck = size(Altitude(Altitude<CheckH),1);
% CheckPhH = PhjSum./PhjSum(KMCheck,:);
% MeanCheck = 2000;
% CheckCof = MeanCheck./PhjSum(KMCheck,:);
% CheckPh = CheckPhH./CheckCof;
% figure('Name','光子数检查','position',[10 100 1350 500])
% for j = 1:size(CheckPh,2)
%     plot(CheckPh(:,j)+j,Altitude,'-','linewidth',1.5);
%     hold on;
% end
% xlim([1193 1198]);
% ylim([85 90]);
% % xlim([-1 41]);
% % ylim([75 115]);
% xInterval = 1;  % x轴刻度间隔
% 
% % 噪声
% Noise = mean(PhjSum(KM120:KM125,:),1);
% SumPh = PhjSum - Noise;
% 
% SNR = SumPh./Noise;
% 
% 
% 
% % 获取瑞利后向散射截面和有效后向散射截面
% ScatterTable = readtable('RayEffScatter.txt');
% Scatter = table2array(ScatterTable(:,2:3));
% RayScatter = Scatter(6,1);
% EffScatter = Scatter(6,2);
% 
% % 获取30 km 处大气模型温度与密度
% [TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,2023,201,17);
% 
% % 反演密度
% Z = Altitude;                       % 高度矩阵
% ZR = 30;                            % 参考高度
% SigmaRay = RayScatter;              % 瑞利后向散射截面
% SigmaEff = EffScatter;              % 有效后向散射截面
% NZ = PhjSum;                        % 光子数矩阵
% NB = Noise;                         % 噪声矩阵
% NZR = PhjSum(KM30,:);                % 参考高度处光子数
% NumberRay = sum(DenRay)*1e-6;       % 参考高度处大气模型数密度
% NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

%% 归一化密度检查
% CheckH = 30;% 光子数检查归一化高度，KM
% KMCheck = size(Altitude(Altitude<CheckH),1);
% CheckPh = NumberZ./NumberZ(KMCheck,:);
% figure('Name','密度检查','position',[10 100 1350 500])
% for j = 1:size(CheckPh,2)
%     plot(CheckPh(:,j)+j,Altitude,'-k','linewidth',1.5);
%     hold on;
% end
% set(gca,'FontSize',12,'FontName','Times New Roman');
% title('Fe Density Profiles 20231102 @Mohe');
% xlabel('Local Time');
% ylabel('Altitude (km)');
% xlim([-1 41]);
% % xlim([120 130]);
% ylim([75 115]);



%%
% FindTime = TimeList';
% Prof = 802;
% plot(Altitude,NumberZ(:,Prof),'-','linewidth',1)
% grid on;
% Profstr = string(Prof)+' | '+TimeList(:,Prof);
% title(Profstr)
% xlim([70 125]);
% % ylim([0 200]);
% xlabel('Altitude (km)');
% ylabel('Density (cm^{-3})');
% legend('K');

%% 尾迹变化图
% TimeTrail = TimeX(:,799:801)';
% DenTrail = [24 54 20];
% AltTrail = [98.27 98.3 98.27];
% Profx = 799;
% 
% yyaxis left
% plot(TimeTrail,DenTrail,'-b','linewidth',1.5)
% % ylim([0 1000])
% ylabel('Density of Trail (cm^{-3})')
% grid on
% yyaxis right
% plot(TimeTrail,AltTrail,'-r','linewidth',1.5)
% ylim([98 98.5])
% ylabel('Altitude of Trail (km)')
% TDstr = 'starttime: '+string(Profx)+' | '+TimeList(:,Profx);
% title(TDstr)
% xlabel('Time')
% datetick('x','HH:MM','keeplimits','keepticks');
% legend('Density','Altitude')

%%

% figure('Name','光子数细致检查','position',[600 100 350 500])
% for j = 1:size(CheckPh,2)
%     plot(CheckPh(:,j)+j,Altitude,'-','linewidth',1.5);
%     hold on;
% end
% title('K 20231106T13:31:47')
% xlim([370 377]);
% ylim([75 115]);
% xInterval = 1;  % x轴刻度间隔
% 
% figure('Name','光子数heng检查','position',[600 300 400 150])
% for j = 1:size(CheckPh,2)
%     plot(CheckPh(:,j)+j,Altitude,'-','linewidth',1.5);
%     hold on;
% end
% title('K 20231106T21:09:30')
% xlim([1193 1199]);
% ylim([80 100]);
% xInterval = 1;  % x轴刻度间隔

% figure('Name','光子数Time检查','position',[300 300 800 250])
% TXcal = (Time(end)-Time(1))/(size(Time,2)-1);
% TXcheck = CheckPh;
% for j = 1:size(CheckPh,2)
%     TXcheck(:,j) = (CheckPh(:,j)*3+j-1).*TXcal+Time(1);
% end
% plot(TXcheck,Altitude,'-','linewidth',1.5);
% hold on;
% title('K 20231103')
% TX1 = double(14+(33/60))+8;
% TX2 = double(14+(43/60))+8;
% % xticks = Time;
% xlim([TX1 TX2]);
% ylim([80 100]);
% set(gca,'FontSize',12,'FontName','Times New Roman');
% xlabel('Time (LT)');
% ylabel('Altitude (km)');





